import sys, numpy as np
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
from smg_dmrg.model_3450 import Model3450Spec, mass_correlation_term, centered_left_position
from smg_dmrg.tenpy_backend import run_ground_dmrg_3450
from smg_dmrg.analyze_correlations import classify
L=12; CHI=384; MAXD=8
print("START largeg", flush=True)
for g in (6.5, 12.0, 20.0):
    sp = Model3450Spec(unit_cells=L, t1=1.0, t2=0.5, g1=g, g2=g)
    res = run_ground_dmrg_3450(sp, chi_max=CHI, max_sweeps=30, svd_min=1e-12, conserve="q3450", combine=True)
    psi = res["psi"]; md = min(MAXD, sp.edge_sites-1)
    rs=[]; mags=[]
    for d in range(1, md+1):
        left = centered_left_position(sp.edge_sites, d)
        v = abs(psi.expectation_value_term(mass_correlation_term("B","dirac",("4","0"),left,d), autoJW=True))
        rs.append(d); mags.append(float(v))
    c = classify(rs, mags)
    print(f"g={g:>5} E={res['energy']:.3f} chi={max(psi.chi)}", flush=True)
    print(f"   dirac_40(r): {[f'{m:.2e}' for m in mags]}", flush=True)
    print(f"   -> {c.get('form')}  {c.get('detail')}", flush=True)
print("DONE largeg", flush=True)
