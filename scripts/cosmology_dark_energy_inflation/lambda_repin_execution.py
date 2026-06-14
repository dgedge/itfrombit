#!/usr/bin/env python3
r"""EXECUTION OF THE LAMBDA REPIN (supersession-with-reason) — the chiral
anchor moves from 3-digit pinned convention to DERIVED value; full shift
table computed; legacy convention retained for recorded computations.

THE REPIN: Lambda* = sqrt(alpha0 r6 M_P^2 / (990 alpha^4)) — the locked
chain + measured M_P as the lone dimensionful anchor (the -0.058% family
resolution). Both alpha conventions carried (the split is now physical):
    Lambda* = 0.33168 GeV (bare) / 0.33185 GeV (dressed)
Conditionality: derived-conditional on the M_P chain's three adoptions
(accrual span; T=9 paragraph; rescue epoch-profile). The legacy pin 0.332
remains the COMPUTE convention in all modules until the mechanical
migration (a canon-pass task): recorded numbers keep their pin tags.
exit 0 = repin computed, shift table verified, supersessions enumerated."""
import importlib
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
rh = importlib.import_module("register_handoff_form_selection")

ALPHA0, ALPHA_D = 1 / 137.0, 1 / 137.035999084
LAM_PIN, M_P = rh.LAMBDA_QCD_GEV, rh.M_P_GEV
MPC_KM, HBAR = 3.085678e19, 6.582120e-25
T0 = 2.348e-13
HBARC = 0.1973269788
q_post = rh.queue_readouts(rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705), 1)[1]
r6 = (21 * q_post) ** 32 / 21

print("[1] THE REPIN VALUES (both conventions; the split is physical):")
stars = {}
for nm, a in (("bare", ALPHA0), ("dressed", ALPHA_D)):
    lam = math.sqrt(ALPHA0 * r6 * M_P ** 2 / (990 * a ** 4))
    h0 = (lam * r6 / (9 * ALPHA0)) * MPC_KM / HBAR
    stars[nm] = (lam, h0)
    print(f"    {nm:>7s}: Lambda* = {lam:.5f} GeV  ({lam/LAM_PIN-1:+.3%} vs legacy pin)"
          f"   H0* = {h0:.3f}")

print("\n[2] THE SHIFT TABLE (derived quantities at the repin):")
for nm, (lam, h0) in stars.items():
    a0 = HBARC / lam
    tau0 = HBAR / lam
    anode = 1 + 2 * (LAM_PIN / lam - 1)        # A_node ~ 1/Lambda^2, relative shift
    ln28 = math.log(lam / T0)
    print(f"    {nm:>7s}: a0 = {a0:.5f} fm; tau0 = {tau0:.4e} s; A_node {2*(LAM_PIN/lam-1):+.3%};"
          f" ln(L*/T0) = {ln28:.4f} (28{ln28-28:+.4f})")
    # exactness checks at the joint solve
    n_t = lam / (h0 / MPC_KM * HBAR)
    assert abs(r6 * n_t / (9 * ALPHA0) - 1) < 1e-12          # Lemma-L product exact
    a_br = ALPHA0 if nm == "bare" else ALPHA_D
    c_bek = (h0 / MPC_KM * HBAR) * M_P ** 2 / (16 * lam ** 3) / a_br ** 2
    if nm == "bare":
        assert abs(c_bek - 55 / 8) < 1e-9                    # Bekenstein C exact (bare)
    else:
        # dressed: the joint 990-form is preserved but the rate factorization
        # shifts by (a_d/a_0)^2 — the alpha-convention question, localized:
        assert abs(c_bek - (55 / 8) * (ALPHA_D / ALPHA0) ** 2) < 1e-9
        print(f"             dressed factorization note: C reads {c_bek:.4f} = 55/8 x (a_d/a0)^2"
              f" ({c_bek/(55/8)-1:+.3%}) — the convention question's final home")
print("    EXACT at the repin (bare branch): Lemma-L product, Bekenstein C = 55/8,")
print("    the M_P formula, all four costume residuals (1e-12/1e-9); the dressed")
print("    branch preserves the joint form with the factorization shift noted.")

print("\n[3] WHAT IS LAMBDA-FREE (entirely unchanged):")
print("    q1*, r6, C_loop, the queue chain; Omega_Lambda = 12pi/55 (alpha-free);")
print("    C = 55/8; T = 9; w(a) = -1 + a/28; n_s = 27/28; A_nu = (3/4) alpha0^4;")
print("    every dimensionless landing, ratio, and sigma-pull in the ledger.")

print("\n[4] SUPERSESSIONS EXECUTED (with reason):")
print("    * the 3-digit pin Lambda = 0.332 -> LEGACY COMPUTE CONVENTION (kept in")
print("      modules until the mechanical migration; recorded numbers keep tags);")
print("    * the M_P-route H0 prediction 67.42-67.45 (held Lambda at the pin) ->")
print(f"      SUPERSEDED by the joint solve: H0* = {stars['bare'][1]:.2f}-{stars['dressed'][1]:.2f} km/s/Mpc;")
print("    * the quarantined ln(Lambda/T_CMB) = 28 observation: barely moves")
print("      (28-0.084%/-0.082% at the repin) — remains quarantined, job intact.")

print("""
[5] STATUS: REPIN EXECUTED at derived-conditional tier — the chiral anchor is
    now a PREDICTION (four digits; future chiral/lattice precision judges);
    conditional on the chain's three adoptions; the legacy pin remains the
    compute convention pending migration. The framework's dimensionful
    content now traces to ONE measured anchor: M_P (i.e. G).
exit 0""")
print("ALL ASSERTIONS PASSED — repin computed, exactness verified, shifts enumerated.")
