import sys, time, numpy as np
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
from smg_dmrg.model_3450 import (Model3450Spec, mass_correlation_term, centered_left_position,
                                 build_3450_terms)
from smg_dmrg.tenpy_backend import build_tenpy_term_model, run_ground_dmrg_3450, build_3450_charged_lattice
from tenpy.algorithms import dmrg
L=16; CHI=512; MAXD=10
def model_for(g):
    sp = Model3450Spec(unit_cells=L, t1=1.0, t2=0.5, g1=g, g2=g)
    _, terms = build_3450_terms(sp)
    lat = build_3450_charged_lattice(sp)
    model = build_tenpy_term_model(sp.total_modes, terms, conserve=None, lattice=lat, site_index_mode="3450")
    return sp, model
def d40(psi, sp):
    md=min(MAXD, sp.edge_sites-1); out=[]
    for d in range(1,md+1):
        left=centered_left_position(sp.edge_sites,d)
        out.append(float(abs(psi.expectation_value_term(mass_correlation_term("B","dirac",("4","0"),left,d),autoJW=True))))
    return out
def run_from(psi, model, chi, sweeps=30):
    p={"mixer":True,"max_sweeps":sweeps,"max_trunc_err":None,"combine":True,"trunc_params":{"chi_max":chi,"svd_min":1e-12}}
    return dmrg.run(psi, model, p)["E"]
print("START", flush=True)
res20 = run_ground_dmrg_3450(Model3450Spec(unit_cells=L,t1=1,t2=0.5,g1=20,g2=20), chi_max=CHI, max_sweeps=24, svd_min=1e-12, conserve="q3450", combine=True)
sp20,_=model_for(20.0)
print(f"g=20 fresh:   E={res20['energy']:.4f}  d40tail={np.mean(d40(res20['psi'],sp20)[-3:]):.2e}", flush=True)
res65 = run_ground_dmrg_3450(Model3450Spec(unit_cells=L,t1=1,t2=0.5,g1=6.5,g2=6.5), chi_max=CHI, max_sweeps=30, svd_min=1e-12, conserve="q3450", combine=True)
sp65,_=model_for(6.5)
print(f"g=6.5 FRESH:  E={res65['energy']:.4f}  d40tail={np.mean(d40(res65['psi'],sp65)[-3:]):.2e}", flush=True)
sp65b, model65 = model_for(6.5)
psi_start = res20['psi'].copy()
E_ad = run_from(psi_start, model65, CHI, sweeps=30)
print(f"g=6.5 SEEDED: E={E_ad:.4f}  d40tail={np.mean(d40(psi_start,sp65b)[-3:]):.2e}", flush=True)
print("DONE", flush=True)
