import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
import os as _os; sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from s1_schwinger_u1 import ground_energy

N, w, m, J = 20, 1.0, 0.5, 0.5
print(f"B: U(1) Casimir scaling via DMRG (reusing S1). N={N}, w={w}, m={m}, J={J}", flush=True)
print("test: does the (screened) string tension scale as sigma_q ~ q^2 (abelian Casimir)?", flush=True)
E0, _, _ = ground_energy(N, w, m, J, {})
sig = {}
for q in (1.0, 2.0, 3.0):
    rs, Vs = [], []
    for r in range(1, 7):
        iL = (N - r) // 2; iR = iL + r
        Er, _, _ = ground_energy(N, w, m, J, {iL: +q, iR: -q})
        rs.append(r); Vs.append(Er - E0)
    rs = np.array(rs, float); Vs = np.array(Vs, float)
    A = np.vstack([rs, np.ones(len(rs))]).T
    (s, b), *_ = np.linalg.lstsq(A, Vs, rcond=None)
    r2 = 1 - np.sum((Vs-(s*rs+b))**2)/np.sum((Vs-Vs.mean())**2)
    sig[q] = s
    print(f"  q={q:.0f}: sigma_q={s:.4f}  R^2={r2:.3f}   V(r)={[round(v,2) for v in Vs]}", flush=True)
print("\n  Casimir-scaling check (expect sigma_q/sigma_1 = q^2 = 1,4,9):", flush=True)
for q in (1.0,2.0,3.0):
    print(f"    sigma_{q:.0f}/sigma_1 = {sig[q]/sig[1.0]:.3f}   (q^2 = {q*q:.0f})", flush=True)
print("DONE", flush=True)
