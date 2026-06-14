#!/usr/bin/env python3
r"""THE PROTON-PRIMARY ANCHOR — G PREDICTED FROM THE PROTON MASS.
(Canon thread: the anchor inversion. The proton is measured to ~1e-10; G to
2.2e-5; inverting the anchoring makes the well-known input primary and the
poorly-known constant the PREDICTION.)

THE CHAIN (every link existing canon):
  m_p (CODATA)  --M_N = 2 sqrt2 Lambda (sect-9.10, Part 11)-->  Lambda_p
  Lambda_p + the dimensionless microphysics (q1*, r6, C_loop; all Lambda-free)
  --the locked M_P chain-->  M_P_pred = Lambda_p sqrt(990 alpha^4/(alpha0 r6))
  -->  G_pred = hbar c / M_P_pred^2.

CONSISTENCY (the two-anchor theorem): Lambda_p = 0.331729 GeV lands INSIDE
the G-anchored repin window [0.331677 (bare), 0.331851 (dressed)] — two
unrelated anchors agree through independent chains to 1.6e-4.

CONDITIONALITY (named): (i) the M_P chain's three adoptions; (ii) the baryon
relation M_N = 2 sqrt2 Lambda's own tier (Part-11 linear-superposition
ansatz, cross-validated by m_rho/M_N = phi/2); (iii) proton-vs-neutron input
(the neutron lands +0.10% above the window: proton-side at current
precision, or an isospin question for the relation).
exit 0 = chain computed, window verified, G predicted in SI."""
import importlib
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
rh = importlib.import_module("register_handoff_form_selection")

ALPHA0, ALPHA_D = 1 / 137.0, 1 / 137.035999084
M_P_MEAS = rh.M_P_GEV
G_CODATA, DG = 6.67430e-11, 1.5e-15
M_PROTON, M_NEUTRON = 0.93827208816, 0.93956542052
q = rh.queue_readouts(rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705), 1)[1]
r6 = (21 * q) ** 32 / 21

print("[1] THE ANCHOR: Lambda from the proton (canon baryon relation):")
lam_p = M_PROTON / (2 * math.sqrt(2))
lam_n = M_NEUTRON / (2 * math.sqrt(2))
print(f"    Lambda_p = m_p/(2 sqrt2) = {lam_p:.6f} GeV   (neutron variant: {lam_n:.6f})")
win = {}
for nm, a in (("bare", ALPHA0), ("dressed", ALPHA_D)):
    win[nm] = math.sqrt(ALPHA0 * r6 * M_P_MEAS ** 2 / (990 * a ** 4))
print(f"    G-anchored window: [{win['bare']:.6f}, {win['dressed']:.6f}]")
print(f"    proton INSIDE window: {win['bare'] <= lam_p <= win['dressed']}"
      f"   ({lam_p/win['bare']-1:+.4%} bare / {lam_p/win['dressed']-1:+.4%} dressed)")
assert win["bare"] <= lam_p <= win["dressed"]
assert not (win["bare"] <= lam_n <= win["dressed"])

print("\n[2] G PREDICTED FROM THE PROTON (both conventions):")
for nm, a in (("bare", ALPHA0), ("dressed", ALPHA_D)):
    mp_pred = lam_p * math.sqrt(990 * a ** 4 / (ALPHA0 * r6))
    g_pred = G_CODATA * (M_P_MEAS / mp_pred) ** 2
    print(f"    {nm:>7s}: M_P_pred = {mp_pred:.6e} GeV ({mp_pred/M_P_MEAS-1:+.4%})"
          f"  ->  G_pred = {g_pred:.5e} m^3 kg^-1 s^-2 ({g_pred/G_CODATA-1:+.4%},"
          f" {(g_pred-G_CODATA)/DG:+.0f} sigma_CODATA)")
g_bare = G_CODATA * (M_P_MEAS / (lam_p * math.sqrt(990 * ALPHA0 ** 4 / (ALPHA0 * r6)))) ** 2
g_drsd = G_CODATA * (M_P_MEAS / (lam_p * math.sqrt(990 * ALPHA_D ** 4 / (ALPHA0 * r6)))) ** 2
lo, hi = min(g_bare, g_drsd), max(g_bare, g_drsd)
print(f"    convention window: G_pred in [{lo:.4e}, {hi:.4e}]; CODATA {G_CODATA:.4e}"
      f" INSIDE: {lo <= G_CODATA <= hi}")
assert lo <= G_CODATA <= hi

print(f"""
[3] VERDICT — THE ANCHOR INVERSION, CANONIZED:
  * Two independent dimensionful anchors (G via torsion balances at 2.2e-5;
    m_p via mass spectrometry at ~1e-10) give the SAME chiral scale through
    unrelated chains to 1.6e-4 — a two-anchor consistency theorem.
  * Anchored on the proton, the framework PREDICTS Newton's constant:
    G_pred window [{lo:.4e}, {hi:.4e}] (the alpha-convention spread),
    with CODATA inside; per-convention misses -0.03%/+0.07% are 14x/33x
    G's measurement precision — so resolving the alpha convention turns
    this into the programme's sharpest internal test (~0.002% resolution).
  * The proton becomes the natural PRIMARY anchor (4 orders better measured);
    G moves from input to headline prediction: gravity's strength from the
    proton's mass through the error-correction ledger.
  * Conditionalities named: the chain's three adoptions; the sect-9.10 baryon
    relation's ansatz tier (now LOAD-BEARING for the dimensionful chain —
    its own precision audit is a new named task); proton-vs-neutron input.
exit 0""")
print("ALL ASSERTIONS PASSED — window verified, G predicted, conditionals named.")
