import sys, numpy as np
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
import smg_dmrg.model_3450 as m
sp = m.Model3450Spec(unit_cells=16, t1=1.0, t2=0.5, g1=6.5, g2=6.5)
modes = m.build_3450_modes(sp)
print("n_modes:", len(modes), " sp.total_modes:", sp.total_modes)
# mode charge attribute
c0 = modes[0]
print("mode0 attrs:", [a for a in dir(c0) if not a.startswith('_')][:20])
print("mode0 charge:", getattr(c0,'charge',None), " index:", getattr(c0,'index',None))
ips = m.initial_product_state_3450(sp) if hasattr(m,'initial_product_state_3450') else m.initial_product_state(sp)
print("len initial state:", len(ips), " sample:", ips[:16])
# filling
def occ(s):
    return 1 if str(s).lower() in ('full','1','occupied','up') or s==1 else 0
nf = sum(occ(s) for s in ips)
print("filled:", nf, " filling fraction:", round(nf/len(ips),4))
# total charge assuming site order == modes order
q = np.zeros(2, dtype=int)
for i,s in enumerate(ips):
    if occ(s) and i < len(modes):
        ch = modes[i].charge
        q += np.array(ch, dtype=int)
print("total (Q1,Q2) of initial state (if order matches):", tuple(q.tolist()))
