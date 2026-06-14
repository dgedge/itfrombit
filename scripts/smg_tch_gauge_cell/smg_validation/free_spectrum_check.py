import sys, numpy as np
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
from smg_dmrg.model_3450 import Model3450Spec, free_spectrum_3450
for L in (8,16,24):
    sp = Model3450Spec(unit_cells=L, t1=1.0, t2=0.5, g1=0.0, g2=0.0)
    ev, herm, gap = free_spectrum_3450(sp)
    ev = np.sort(np.real(np.asarray(ev)))
    absev = np.sort(np.abs(ev))
    nz = absev[absev < 0.6]
    print(f"L={L:>2}: n_modes={len(ev)} hermOK={herm} gap={gap:.4f} gap*L={gap*L:.2f} #|E|<0.6={len(nz)}")
    print(f"      lowest |E|: {np.round(absev[:10],4)}")
