#include <math.h>
static double T0=0.0;
void set_t0(double v){T0=v;}
static double D=0.08,D2=0.03,F0=0.12,G1=0.50,G3=0.08,S2=1.593*1.593,W=1.72,GAM=3.0,KK=8.0;
static void driven_rhs(double drive,const double*z,int pl,double*o){
  double x=z[0],xd=z[1],r=z[2],rd=z[3];
  double F=drive*(F0+(G1-G3*x*x)*x);
  o[0]=xd; o[1]=x-x*x*x-D*xd-D2*xd*fabs(xd)+F-GAM*x*r*r;
  if(pl){o[2]=0;o[3]=0;} else {o[2]=rd; o[3]=S2*r-KK*r*r*r-GAM*x*x*r-D*rd-D2*rd*fabs(rd);}
}
static void rhs(double t,const double*z,double A,int pl,double*o){
  driven_rhs(A*cos(W*t),z,pl,o);
}
static void step(double t,double*z,double h,double A,int pl){
  double k1[4],k2[4],k3[4],k4[4],y[4];int i;
  rhs(t,z,A,pl,k1); for(i=0;i<4;i++)y[i]=z[i]+h/2*k1[i];
  rhs(t+h/2,y,A,pl,k2); for(i=0;i<4;i++)y[i]=z[i]+h/2*k2[i];
  rhs(t+h/2,y,A,pl,k3); for(i=0;i<4;i++)y[i]=z[i]+h*k3[i];
  rhs(t+h,y,A,pl,k4); for(i=0;i<4;i++)z[i]+=h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
}
/* integrate n states (n x 4, row-major) for nper periods, steps per period; optional record of x at every period end into rec (n x nper) and xd into rec2 */
void integrate(double*z,int n,double A,int nper,int steps,int pl,double*recx,double*recv,double*recr){
  double T=2*M_PI/W,h=T/steps;
  #pragma omp parallel for
  for(int i=0;i<n;i++){
    double*zi=z+4*i,t=T0;
    for(int p=0;p<nper;p++){
      for(int s=0;s<steps;s++){step(t,zi,h,A,pl);t+=h;}
      if(recx){recx[(long)i*nper+p]=zi[0];recv[(long)i*nper+p]=zi[1];recr[(long)i*nper+p]=zi[2];}
    }
  }
}
/* Benettin largest Lyapunov exponent for full 4D (or planar 2D if pl) from z; returns exponent per unit time and also normal exponent via separate tangent in rho-direction */
double lyap(double*z0,double A,int ntrans,int nper,int steps,int pl,double d0){
  double T=2*M_PI/W,h=T/steps,t=0,z[4],w[4];int i;
  for(i=0;i<4;i++)z[i]=z0[i];
  for(int p=0;p<ntrans;p++)for(int s=0;s<steps;s++){step(t,z,h,A,pl);t+=h;}
  for(i=0;i<4;i++)w[i]=z[i]; w[0]+=d0; if(!pl) w[2]+=0; 
  double sum=0;
  for(int p=0;p<nper;p++){
    for(int s=0;s<steps;s++){step(t,z,h,A,pl);step(t,w,h,A,pl);t+=h;}
    double d=0;for(i=0;i<4;i++)d+=(w[i]-z[i])*(w[i]-z[i]); d=sqrt(d);
    sum+=log(d/d0); for(i=0;i<4;i++)w[i]=z[i]+(w[i]-z[i])*d0/d;
  }
  return sum/(nper*T);
}
/* normal (rho) Lyapunov exponent: growth of linearised rho-equation along planar trajectory */
double nlyap(double*z0,double A,int ntrans,int nper,int steps){
  double T=2*M_PI/W,h=T/steps,t=0,z[4]={z0[0],z0[1],0,0};
  for(int p=0;p<ntrans;p++)for(int s=0;s<steps;s++){step(t,z,h,A,1);t+=h;}
  double u=1e-6,ud=0,sum=0;
  for(int p=0;p<nper;p++){
    for(int s=0;s<steps;s++){
      /* RK4 on (x,xd,u,ud) with linear u-eq */
      double y[4]={z[0],z[1],u,ud},k1[4],k2[4],k3[4],k4[4],q[4];int i;
      #define LIN(tt,yy,o) { double F=A*(F0+(G1-G3*yy[0]*yy[0])*yy[0])*cos(W*(tt)); o[0]=yy[1]; o[1]=yy[0]-yy[0]*yy[0]*yy[0]-D*yy[1]-D2*fabs(yy[1])*yy[1]+F; o[2]=yy[3]; o[3]=(S2-GAM*yy[0]*yy[0])*yy[2]-D*yy[3]-D2*fabs(yy[3])*yy[3]; }
      LIN(t,y,k1); for(i=0;i<4;i++)q[i]=y[i]+h/2*k1[i];
      LIN(t+h/2,q,k2); for(i=0;i<4;i++)q[i]=y[i]+h/2*k2[i];
      LIN(t+h/2,q,k3); for(i=0;i<4;i++)q[i]=y[i]+h*k3[i];
      LIN(t+h,q,k4); for(i=0;i<4;i++)y[i]+=h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
      z[0]=y[0];z[1]=y[1];u=y[2];ud=y[3];t+=h;
    }
    double nrm=sqrt(u*u+ud*ud); sum+=log(nrm/1e-6); u*=1e-6/nrm; ud*=1e-6/nrm;
  }
  return sum/(nper*T);
}

/* First hitting of x >= 0 at an accepted RK4 step. Each start begins in x < 0.
   Precompute the time-periodic drive once per phase, shared by all starts. */
#include <stdint.h>
#include <stdlib.h>
void first_cross(const double *x, const double *v, int count, double A,
                 int planar, int steps, uint8_t *hit) {
    const double h = 2*M_PI/(W*steps);
    double drive[2*steps+1];
    for (int k=0; k<=2*steps; k++)
        drive[k] = A*cos(M_PI*k/steps);
    #pragma omp parallel for schedule(dynamic,32)
    for (int i=0; i<count; i++) {
        double z[4]={x[i],v[i],planar ? 0.0 : 1e-4,0.0};
        hit[i]=0;
        for (int t=0; t<20*steps && !hit[i]; t++) {
            int k=2*(t%steps);
            double a[4],b[4],c[4],d[4],y[4];
            driven_rhs(drive[k],z,planar,a);
            for (int j=0;j<4;j++) y[j]=z[j]+h*a[j]/2;
            driven_rhs(drive[k+1],y,planar,b);
            for (int j=0;j<4;j++) y[j]=z[j]+h*b[j]/2;
            driven_rhs(drive[k+1],y,planar,c);
            for (int j=0;j<4;j++) y[j]=z[j]+h*c[j];
            driven_rhs(drive[k+2],y,planar,d);
            for (int j=0;j<4;j++) z[j]+=h*(a[j]+2*b[j]+2*c[j]+d[j])/6;
            hit[i]=(z[0]>=0);
        }
    }
}
