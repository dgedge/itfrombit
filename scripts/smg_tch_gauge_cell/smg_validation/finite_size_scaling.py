import sys, time, numpy as np
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
from smg_dmrg.model_3450 import Model3450Spec, mass_correlation_term, centered_left_position
from smg_dmrg.tenpy_backend import run_ground_dmrg_3450
from smg_dmrg.analyze_correlations import classify
GS = (8.0, 10.0)        # safely above g_c~5.7, gapped phase (out of the critical fan)
LS = (8, 12, 16)
CHI = 768               # gapped mirror uses far less; cap driven only by the gapless light edge
print("START scaling", flush=True)
for g in GS:
    print(f"\n### g={g}  (g_c~5.7; gapped) ###", flush=True)
    for L in LS:
        t0=time.time()
        sp = Model3450Spec(unit_cells=L, t1=1.0, t2=0.5, g1=g, g2=g)
        res = run_ground_dmrg_3450(sp, chi_max=CHI, max_sweeps=30, svd_min=1e-12, conserve="q3450", combine=True)
        psi = res["psi"]; md = sp.edge_sites - 1
        rs=[]; mags=[]
        for d in range(1, md+1):
            left = centered_left_position(sp.edge_sites, d)
            mags.append(float(abs(psi.expectation_value_term(mass_correlation_term("B","dirac",("4","0"),left,d), autoJW=True))))
            rs.append(d)
        c = classify(rs, mags)
        print(f"L={L:>2} chi_used={max(psi.chi):>4} E={res['energy']:.3f} maxdist={md} ({time.time()-t0:.0f}s)", flush=True)
        print(f"   dirac_40(r): {[f'{m:.2e}' for m in mags]}", flush=True)
        print(f"   largest-r value = {mags[-1]:.2e}   [{c.get('form')}; {c.get('detail')}]", flush=True)
print("DONE scaling", flush=True)
