import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
import os as _os; sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from s1_schwinger_u1 import ground_energy
N, w, m, J = 20, 1.0, 0.5, 0.5
print(f"B(frac): U(1) Casimir scaling, FRACTIONAL (unscreened) charges. N={N}", flush=True)
E0,_,_ = ground_energy(N,w,m,J,{})
sig={}
for q in (0.3,0.6,0.9):
    rs,Vs=[],[]
    for r in range(1,7):
        iL=(N-r)//2; iR=iL+r
        Er,_,_ = ground_energy(N,w,m,J,{iL:+q,iR:-q})
        rs.append(r); Vs.append(Er-E0)
    rs=np.array(rs,float); Vs=np.array(Vs,float)
    A=np.vstack([rs,np.ones(len(rs))]).T
    (s,b),*_=np.linalg.lstsq(A,Vs,rcond=None)
    r2=1-np.sum((Vs-(s*rs+b))**2)/np.sum((Vs-Vs.mean())**2)
    sig[q]=s
    print(f"  q={q}: sigma={s:.5f} R^2={r2:.3f}", flush=True)
print("\n  Casimir scaling (expect sigma_q/sigma_0.3 = (q/0.3)^2 = 1,4,9):", flush=True)
for q in (0.3,0.6,0.9):
    print(f"    sigma_{q}/sigma_0.3 = {sig[q]/sig[0.3]:.3f}   ((q/0.3)^2 = {(q/0.3)**2:.2f})", flush=True)
print("DONE", flush=True)
