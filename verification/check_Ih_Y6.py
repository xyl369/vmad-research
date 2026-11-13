import numpy as np, warnings
warnings.filterwarnings("ignore")
from scipy.special import sph_harm
def Y(l,m,th,ph): return sph_harm(m,l,ph,th)
def f(th,ph,s):
    return np.sqrt(11/25)*Y(6,0,th,ph)+s*np.sqrt(7/25)*(Y(6,5,th,ph)-Y(6,-5,th,ph))
t0=np.arctan(2)
for off in [0,np.pi/5]:
  V=[(0,0),(np.pi,0)]+[(t0,2*np.pi*k/5+off) for k in range(5)]+[(np.pi-t0,2*np.pi*k/5+np.pi/5+off) for k in range(5)]
  for s in [1,-1,1j,-1j]:
    vals=np.array([f(th,ph,s) for th,ph in V])
    # invariance test under 3-fold rotation about a face normal: random points
    print("off",round(off,3),"s",s,"vertex spread",np.ptp(vals.real).round(6),np.abs(vals.imag).max().round(6))
