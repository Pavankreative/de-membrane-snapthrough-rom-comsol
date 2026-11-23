import ctypes, os, subprocess, numpy as np
_here = os.path.dirname(os.path.abspath(__file__)); _so = os.path.join(_here, 'librk.so')
if not os.path.exists(_so):   # compile on first use (needs gcc)
    subprocess.check_call(['gcc', '-O2', '-shared', '-fPIC', '-fopenmp', os.path.join(_here, 'rk.c'), '-o', _so, '-lm'])
lib=ctypes.CDLL(_so)
dp=ctypes.POINTER(ctypes.c_double)
lib.integrate.argtypes=[dp,ctypes.c_int,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_int,dp,dp,dp]
lib.lyap.restype=ctypes.c_double; lib.lyap.argtypes=[dp,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double]
lib.nlyap.restype=ctypes.c_double; lib.nlyap.argtypes=[dp,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_int]
def P_(a): return a.ctypes.data_as(dp)
def integ(z0,A,nper,steps=200,planar=False,record=False):
    z=np.ascontiguousarray(z0.T.copy()); n=z.shape[0]
    if record:
        rx=np.zeros((n,nper));rv=np.zeros((n,nper));rr=np.zeros((n,nper))
        lib.integrate(P_(z),n,A,nper,steps,int(planar),P_(rx),P_(rv),P_(rr)); return z.T,(rx,rv,rr)
    lib.integrate(P_(z),n,A,nper,steps,int(planar),None,None,None); return z.T,None
def lyap(z0,A,ntrans,nper,steps=200,planar=False,d0=1e-8):
    z=np.array(z0,float); return lib.lyap(P_(z),A,ntrans,nper,steps,int(planar),d0)
def nlyap(z0,A,ntrans,nper,steps=200):
    z=np.array(z0,float); return lib.nlyap(P_(z),A,ntrans,nper,steps)
