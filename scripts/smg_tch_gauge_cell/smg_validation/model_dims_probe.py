import sys, time
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
from smg_dmrg.model_3450 import Model3450Spec
from smg_dmrg.tenpy_backend import (build_3450_terms, build_3450_charged_lattice,
                                    build_tenpy_term_model, charge_violations, is_3450_charge_conserve)

for L in (8, 16):
    sp = Model3450Spec(unit_cells=L, t1=1.0, t2=0.5, g1=6.5, g2=6.5)
    modes, terms = build_3450_terms(sp)
    lat = build_3450_charged_lattice(sp)
    model = build_tenpy_term_model(sp.total_modes, terms, conserve=None, lattice=lat, site_index_mode="3450")
    sites = model.lat.mps_sites()
    dims = [s.dim for s in sites]
    wmpo = model.H_MPO.chi  # list of MPO bond dims
    print(f"L={L}: total_sites={sp.total_modes} n_mps_sites={len(sites)} site_dims(min/max)={min(dims)}/{max(dims)}")
    print(f"       MPO bond dim W: min={min(wmpo)} max={max(wmpo)}  (n_terms={len(terms)})")
