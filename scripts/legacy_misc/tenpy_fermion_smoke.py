#!/usr/bin/env python3
"""
Verified fermion-DMRG smoke test for the chiral-SMG DMRG programme (see smg_dmrg_scope.md).
Confirms the TeNPy fermion + U(1)-charge-conservation pipeline works and matches exact theory.

Run:  ~/tenpy-env/bin/python python_code/tenpy_fermion_smoke.py
(TeNPy lives in ~/tenpy-env:  python3 -m venv ~/tenpy-env && ~/tenpy-env/bin/pip install physics-tenpy)

Model: spinless free-fermion hopping chain, particle number conserved (conserve='N'), with a chi
RAMP. The ground state is exactly solvable, so DMRG must match to ~1e-14. Verified (L=20, half
filling): DMRG E0 = exact E0 = -12.381490 (diff 1.8e-14), <N>=10, max chi=100.

KEY DMRG LESSONS this encodes (the difference between a light run and a heavy one):
  * ALWAYS conserve the charge (conserve='N'): dense tensors split into small symmetry blocks
    -> ~10-100x less time/memory at the SAME chi. conserve=None is the slow/heavy choice.
  * RAMP chi via chi_list; never jump straight to a huge chi_max.
  * The SMG phase is GAPPED -> modest chi (hundreds to low thousands) suffices. chi~6000 is a
    'something is wrong / sitting at criticality' number, not a default. Watch whether chi stays
    pinned at chi_max (raise it) or settles below it (converged).

This is the exact building block for the 3-4-5-0 kinetic sector (arXiv:2202.12355): swap the
free hopping for the multilayer Chern-ladder kinetic term, add the six-fermion edge-B interaction,
keep conserve + the chi ramp, and measure the SMG-vs-SSB observables (see smg_dmrg_scope.md sec 3).

Requires: numpy + physics-tenpy.
"""
import numpy as np
import tenpy
from tenpy.models.fermions_spinless import FermionChain
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

L = 20
M = FermionChain({'L': L, 'J': 1.0, 'V': 0.0, 'mu': 0.0, 'bc_MPS': 'finite', 'conserve': 'N'})
psi = MPS.from_product_state(M.lat.mps_sites(), ['full', 'empty'] * (L // 2), 'finite')
eng = dmrg.TwoSiteDMRGEngine(psi, M, {
    'trunc_params': {'chi_max': 100, 'svd_min': 1e-10},
    'chi_list': {0: 16, 5: 64, 10: 100},     # ramp chi over sweeps
    'max_sweeps': 20, 'mixer': True})
E, psi = eng.run()

# exact free-fermion ground state (OBC): fill the lowest L/2 single-particle levels
e_k = np.sort([-2 * 1.0 * np.cos(k * np.pi / (L + 1)) for k in range(1, L + 1)])
E_exact = e_k[:L // 2].sum()

print(f"\nFERMION SMOKE (conserve='N'):")
print(f"  DMRG E0  = {E:.6f}")
print(f"  exact E0 = {E_exact:.6f}   diff = {abs(E - E_exact):.2e}")
print(f"  <N>      = {psi.expectation_value('N').sum():.3f}  (expect {L // 2})")
print(f"  max chi  = {max(psi.chi)}")
assert abs(E - E_exact) < 1e-8, "DMRG did not match exact free-fermion energy"
print("  PASS: DMRG matches exact free-fermion theory; fermion + conserve pipeline works.")
