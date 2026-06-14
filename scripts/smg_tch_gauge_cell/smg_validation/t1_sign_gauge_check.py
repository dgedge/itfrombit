import sys, numpy as np
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
from smg_dmrg.model_3450 import Model3450Spec, free_spectrum_3450
sp_p = Model3450Spec(unit_cells=16, t1=1.0,  t2=0.5, g1=0.0, g2=0.0)
sp_m = Model3450Spec(unit_cells=16, t1=-1.0, t2=0.5, g1=0.0, g2=0.0)
ep = np.sort(np.real(np.asarray(free_spectrum_3450(sp_p)[0])))
em = np.sort(np.real(np.asarray(free_spectrum_3450(sp_m)[0])))
print("signed spectra identical (+t1 vs -t1)? ", np.allclose(ep, em))
print("is -t1 spectrum the NEGATION of +t1 (E->-E)? ", np.allclose(ep, -em[::-1]))
print("sum(+t1 negative-energy modes) =", round(float(ep[ep<0].sum()),4))
print("sum(-t1 negative-energy modes) =", round(float(em[em<0].sum()),4))
