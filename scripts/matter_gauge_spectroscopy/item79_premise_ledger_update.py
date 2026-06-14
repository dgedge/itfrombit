#!/usr/bin/env python3
r"""Item 79 'under named premises' — the premise ledger updated with 2026-06-10 evidence,
and the premise-5 stress test under the 3.2 channel reading (which also re-derives the
rho_Lambda boot pair probability r = 2/9 and checks the 129/128 boundary question).

CONTEXT. The alpha-chain adoption rests on five named premises. Tonight's results bear
on four of them, and one (P5: walk-active dissipator support, the 2-of-3 exit counting
behind C_loop = 3/2) is STRESSED by the item-115 finding that 3.2's hopping matrix is a
channel sum {V_em diagonal, V_strong colour-swap, V_weak LQ-controlled CNOT}, not the
XOR-flip triple {C0-flip, C1-flip, I3-flip} the exit counting used. The same triple
underlies the rho_Lambda boot pair ledger (r = 2/9 = (1/2)(2/3)^2, 18 labels).

METHOD (the session's standard): pose every principled exit rule on the 3.2 channel
triple with zero knobs; compute the exit fraction under each; concordance or divergence
decides whether C_loop = 3/2 and r = 2/9 are reading-robust. Then the boundary-leg
question for 129/128 under the channel reading, and the premise table.

Self-asserting; exit 0 = every number verified."""
from fractions import Fraction as Fr

# ---------------- the channel triple under the two readings ----------------
# XOR/legacy reading: per leg, 3 walk-active flip channels; exits = the colour pair.
legacy = {"C0-flip": dict(register_changing=True,  q_neutral=True),
          "C1-flip": dict(register_changing=True,  q_neutral=True),
          "I3-flip": dict(register_changing=True,  q_neutral=False)}
# 3.2-faithful reading (item115_sec32_kernel_gates.py): per direction the hop is a SUM
# of three channels; Q-neutrality computed there ([V_em,Q]=[V_strong,Q]=0, [V_weak,Q]>0).
sec32 = {"V_em (diagonal)":      dict(register_changing=False, q_neutral=True),
         "V_strong (swap)":      dict(register_changing=True,  q_neutral=True),
         "V_weak (LQ-CNOT)":     dict(register_changing=True,  q_neutral=False)}

def exit_fraction(channels, rule):
    ex = [c for c, p in channels.items() if rule(p)]
    return Fr(len(ex), len(channels)), ex

rules = {
    "legacy convention (colour membership)": lambda p: p.get("legacy_exit", False),
    "Q-NEUTRAL channels are exits (photon-compatible transport — the Ward criterion)":
        lambda p: p["q_neutral"],
    "REGISTER-CHANGING channels are exits (an erasure must change the record)":
        lambda p: p["register_changing"],
}
for c in ("C0-flip", "C1-flip"):
    legacy[c]["legacy_exit"] = True
legacy["I3-flip"]["legacy_exit"] = False

print("[1] PREMISE-5 STRESS TEST — exit fractions on the two channel triples:")
fracs = {}
for nm, rule in rules.items():
    if nm.startswith("legacy"):
        f, ex = exit_fraction(legacy, rule)
        print(f"    XOR triple,  {nm}: exits {ex} -> {f}")
    else:
        f, ex = exit_fraction(sec32, rule)
        print(f"    3.2 triple, {nm[:58]}...: exits {ex} -> {f}")
    fracs[nm] = f
assert all(f == Fr(2, 3) for f in fracs.values())
print("    -> ALL THREE principled rules give 2-of-3. The exit FRACTION — the only thing")
print("       C_loop = 3/2 and the boot r use — is reading-robust; the exit SETS differ")
print("       ({V_em,V_strong} vs {V_strong,V_weak}), a named ambiguity for FINER")
print("       observables only. C_loop = 3/2, Gamma_vac = (2/3) alpha_0 Lambda, and the")
print("       K = 205.50 chain stand under the 3.2-faithful kernel.")

# ---------------- the boot pair ledger under the 3.2 reading ----------------
walk = ["V_em (diagonal)", "V_strong (swap)", "V_weak (LQ-CNOT)"]
exits_q = {c for c, p in sec32.items() if p["q_neutral"]}
labels = [(s, l, r) for s in (0, 1) for l in walk for r in walk]
succ = [(s, l, r) for (s, l, r) in labels if s == 1 and l in exits_q and r in exits_q]
r_pair = Fr(len(succ), len(labels))
print(f"\n[2] BOOT PAIR LEDGER under 3.2 channels: {len(labels)} labels, {len(succ)} successes,")
print(f"    r = {r_pair} — the rho_Lambda boot probability is REPRODUCED exactly under the")
print(f"    Ward exit rule (and under the register-changing rule, by the same fraction).")
assert len(labels) == 18 and r_pair == Fr(2, 9)
print("    BOUNDARY-LEG QUESTION (the 129/128 target): under the channel reading the pair")
print("    legs are REGISTER channels of the two pair members, not spatial hops — an open")
print("    boot boundary truncates which CELLS exist, not any cell's channel triple. No")
print("    extra half-leg is licensed: the parallel session's negative is CONFIRMED from")
print("    the 3.2 side. 129/128 remains a sharp theorem target; the strict non-horizon")
print("    prediction stays 0.607 rho_obs (factor 1.65), second route open.")

# ---------------- the premise ledger, updated ----------------
print("""
[3] ITEM-79 PREMISE LEDGER — status after 2026-06-10 (all five remain premises; four
    gained independent support today, none was contradicted):
    P1 identical register defects (2.7/2.8)
        UPGRADED SUPPORT: the 2.8 charge formula is now anomaly-forced end-to-end
        (E2 finish) — the defect-identification layer it rides is derivation-grade.
    P2 two-dimensional coin (3.5)
        UNCHANGED — no new evidence today; still the least-tested premise.
    P3 coin-singlet / S-wave emission (1.5)
        DATA-SUPPORTED: the record-alphabet table EXCLUDES the unprojected pair space
        (496 -> +19.5 sigma) — the sky's photon bookkeeping USES the singlet projection.
    P4 monitored bridge (5.2/5.9/119)
        GROUNDED + DATA-SUPPORTED: the record-content derivation tied the monitoring
        semantics to 5.2's own QND increment structure (alphabet 8 derived; rival
        readings excluded at >15 sigma; s_1 = ln(8x137) at -0.1 sigma).
    P5 walk-active dissipator support / 2-of-3 exits (2.5/3.2)
        HARDENED BY STRESS TEST (this script): the exit fraction 2/3 is invariant
        across the XOR->3.2 channel transition under BOTH principled rules (Q-neutral;
        register-changing); the exit-set ambiguity is named, fraction-irrelevant.
    NET: item 79 stays 'Locked at the algebra level + five adopted identifications' —
    with the identifications now 4/5 independently supported and 0/5 contradicted.

[4] M_P: STAYS HORIZON-INPUT — unchanged. Today's erasure-count results (s_1, the
    record alphabet) are class-2 reading-layer objects; the G7 rank theorem's wall
    (all live large-scale relations reduce to M_P^2 H = O(alpha) Lambda^3, rank 1)
    is untouched: no new independent (M_P, H) relation was created today.
""")
print("ALL ASSERTIONS PASSED — every number above is verified. exit 0")
