"""Create the published PNGs solely from the archived CSV/NPZ outputs."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from common17 import DATA, FIG, T, read_rows

plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.spines.top':False,
                     'axes.spines.right':False,'savefig.facecolor':'white','axes.axisbelow':True})
BLUE='#235b83'; RED='#a23f35'
def save(fig,name):
    fig.savefig(FIG/name,dpi=235)
    plt.close(fig)

def basins():
    summaries=read_rows('basin_summary.csv')
    with np.load(DATA/'basin_labels.npz') as d:
        x,v,m=d['x'],d['v'],d['eligible']
        fig,axes=plt.subplots(2,4,figsize=(10,5.1),sharex=True,sharey=True,layout='constrained')
        curves,cax=plt.subplots(1,2,figsize=(8,4.6),sharey=True,layout='constrained')
        for i,(row,ax) in enumerate(zip(summaries,axes.flat)):
            labels=d[f'case{i}_s200']; hits=int(labels[m].sum())
            assert hits==int(row['escaped'])
            ax.pcolormesh(x,v,np.where(m,labels+1,0),cmap=ListedColormap(['#e5e7e8','white','#17242b']),vmin=0,vmax=2,shading='nearest',rasterized=True)
            ax.set_title(f"{row['model'].capitalize()}  $A={row['A']}$\n{100*hits/m.sum():.1f}% cross",fontsize=9)
            ax.label_outer(); ax.set_xlabel('$x(0)$') if i>=4 else None
            if i%4==0: ax.set_ylabel("$x'(0)$")
            offsets=np.array([1,2,4,8]); delta=offsets*(x[1]-x[0])
            frac=[]
            for k in offsets:
                pairs=m[:,k:] & m[:,:-k]
                frac.append(np.count_nonzero((labels[:,k:]!=labels[:,:-k])&pairs)/pairs.sum())
            alpha,intercept=np.polyfit(np.log(delta),np.log(frac),1)
            assert abs(alpha-float(row['alpha']))<1e-12
            panel=cax[int(row['model']=='planar')]
            line,=panel.loglog(delta,frac,'o',ms=4,label=f"A={row['A']}, α={alpha:.3f}")
            panel.loglog(delta,np.exp(intercept)*delta**alpha,color=line.get_color(),lw=1)
        for ax,title in zip(cax,['Coupled','Planar']):
            ax.set_title(title); ax.grid(alpha=.25,which='both'); ax.set_xlabel('Initial displacement offset, ε')
            ax.legend(fontsize=8)
        cax[0].set_ylabel('Fraction of eligible pairs with different outcomes')
        save(fig,'fig5_first_crossing.png');save(curves,'fig6_uncertainty.png')

def manifolds():
    points=read_rows('manifold_splitting_table3.csv')
    fig,axes=plt.subplots(1,2,figsize=(9.5,3.8),layout='constrained')
    for s in read_rows('table3.csv'):
        name=s['loop']; color=BLUE if name=='right' else RED
        r=[p for p in points if p['loop']==name and p['N_phase']=='32']
        amps=np.unique([float(p['A']) for p in r])
        section=[p for p in r if float(p['A'])==amps[-1]]
        x=np.array([float(p['phase_rad']) for p in section]+[2*np.pi])
        y=np.array([float(p['d']) for p in section]+[float(section[0]['d'])])
        axes[0].plot(x/(2*np.pi),y,'o-',color=color,ms=2.7,label=f"{name.capitalize()}, A={amps[-1]:g}")
        maxima=[max(float(p['d']) for p in r if float(p['A'])==a) for a in amps]
        axes[1].plot(amps/float(s['first_order']),maxima,'o-',color=color,ms=3,label=name.capitalize())
    for ax in axes:
        ax.axhline(0,lw=.7,color='black');ax.legend();ax.grid(alpha=.15)
    axes[0].set(xlabel='Forcing phase / 2π',ylabel='Signed splitting, d',title='Primary turning section')
    axes[1].set(xlabel='Forcing / first-order threshold',ylabel='Maximum sampled splitting',title='Direct onset on each loop')
    save(fig,'fig7_manifold_splitting.png')

def sections():
    rows=read_rows('table4.csv')
    fig,axes=plt.subplots(2,4,figsize=(10,5),layout='constrained')
    with np.load(DATA/'poincare_table4.npz') as data:
        for ax,r in zip(axes.flat,rows):
            x,v=data[f"{r['model']}_{float(r['A'])}"]
            ax.scatter(x.ravel(),v.ravel(),s=1.0,alpha=.5,color=BLUE if r['model']=='planar' else RED,linewidths=0)
            ax.set_title(f"{r['model'].capitalize()}  A={r['A']}\nλ = {float(r['median']):+.3f}",fontsize=9)
            ax.set_ylim(-1.5,1.5);ax.set_xlim(-1.7,1.7)
            if ax.get_subplotspec().is_first_col(): ax.set_ylabel("$x'(kT)$")
            if ax.get_subplotspec().is_last_row(): ax.set_xlabel('$x(kT)$')
    save(fig,'fig8_poincare.png')

def trajectories():
    rows=read_rows('table5.csv')
    fig,axes=plt.subplots(2,2,figsize=(8,6),layout='constrained')
    with np.load(DATA/'table5_traces.npz') as data:
        for col,A in enumerate((.3,.5)):
            i=next(i for i,r in enumerate(rows) if float(r['Knn'])==8 and float(r['A'])==A)
            # Velocity for the portrait is archived from the ODE solver.
            t,x,r=data[f'run{i}']; v=data[f'velocity{i}']
            axes[0,col].plot(x,v,lw=.5,color=BLUE)
            axes[0,col].set(xlabel='x',ylabel="x′",title=f'A = {A:.2f}')
            axes[1,col].plot(t/T,x,lw=.55,color=BLUE,label='x')
            axes[1,col].plot(t/T,r,lw=.55,color=RED,label='ρ')
            axes[1,col].set(xlabel='Forcing periods',ylabel='Displacement');axes[1,col].legend(ncol=2)
    save(fig,'fig1_trajectories.png')
    fig,ax=plt.subplots(figsize=(7.5,4.6),layout='constrained')
    rows=read_rows('figure4_sweep.csv')
    for model,style in [('coupled','o-'),('planar','s--')]:
        r=[p for p in rows if p['model']==model]
        ax.plot([float(p['A']) for p in r],[max(.5,int(p['crossings'])) for p in r],style,label=model.capitalize(),ms=4)
    for r in read_rows('table3.csv'):
        ax.axvline(float(r['first_order']),ls=':' if r['loop']=='right' else '--',color='.5',lw=.9,label=f"First-order {r['loop']}")
    ax.set(yscale='log',xlabel='Forcing amplitude, A',ylabel='Sign changes (zero displayed at 0.5)')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    save(fig,'fig4_crossing_counts.png')

if __name__=='__main__':
    trajectories(); basins(); manifolds(); sections()
