#!/usr/bin/env python3
"""Item 123: boot-sequence residual-fault theorem attempt for rho_Lambda.

Question
--------
Can the horizon-consuming Dirac input in the cosmological-constant row be
replaced by a code-native residual-fault theorem?

Candidate theorem
-----------------
For a distance-4 [8,4,4] boot code whose leading logical failure at each
concatenation level is a malignant pair of lower-level faults,

    q_{ell+1} = q_ell^2 / p_th,
    q_0       = p_th r,

so

    q_ell = p_th r^(2^ell).

The candidate physical subthreshold ratio is the vacuum pair-fault ledger

    r = (coin-singlet pair factor) * (walk-active exiting fraction)^2
      = (1/2) * (2/3)^2
      = 2/9.

The code-native depth candidate is ell=6: equivalently C(4,2) stabilizer-pair
depth, the six coordinate faces of Q3, and the nearest integer to the exact
depth required by the observed rho_Lambda once r=2/9 is inserted.

This script verifies the algebra and applies the ANCHOR 16.4 T1-T9 gates.  The
result is deliberately not promoted to a parameter-free cosmological-constant
derivation: the boot-pair Kraus map, independence/covariance, and exact
coefficient are still open.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction


ALPHA0 = 1.0 / 137.0
LAMBDA_QCD_GEV = 0.332
P_TH = 0.110

# Same observational constants used by horizon_closure_scan.py.
M_P_GEV = 1.220890e19
HBAR_GEV_S = 6.582120e-25
H0_KM_S_MPC = 67.36
MPC_KM = 3.085678e19
OMEGA_L = 0.6847

BIT_NAMES = ["G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W"]
WALK_ACTIVE = [3, 4, 5]  # C0, C1, I3, as in item119_jump_channel.py.


class Status(Enum):
    PASS = "PASS"
    CONDITIONAL = "CONDITIONAL"
    OPEN = "OPEN"
    FAIL = "FAIL"


@dataclass(frozen=True)
class Gate:
    code: str
    name: str
    status: Status
    finding: str


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def rm13_words() -> set[tuple[int, ...]]:
    """RM(1,3), the classical [8,4,4] code used throughout canon."""

    def bits(x: int) -> tuple[int, int, int]:
        return ((x >> 2) & 1, (x >> 1) & 1, x & 1)

    points = range(8)
    gens = [tuple(1 for _ in points)]
    gens += [tuple(bits(x)[k] for x in points) for k in range(3)]
    words = set()
    for coeffs in itertools.product((0, 1), repeat=4):
        words.add(tuple(sum(c * g[i] for c, g in zip(coeffs, gens)) % 2 for i in points))
    return words


def code_distance(words: set[tuple[int, ...]]) -> int:
    weights = [sum(word) for word in words if any(word)]
    return min(weights)


def valid_cell_state(c: tuple[int, ...]) -> bool:
    G0, G1, LQ, C0, C1, I3, chi, W = c
    return (not (G0 and G1)) and W == chi and ((LQ == 0) == (C0 == 0 and C1 == 0))


def flip(index: int, bit: int) -> int:
    return index ^ (1 << (7 - bit))


def derive_walk_active_exit_fraction() -> Fraction:
    all_states = list(itertools.product((0, 1), repeat=8))
    in_p = [valid_cell_state(c) for c in all_states]
    vacuum = all_states.index((0,) * 8)
    exits = [bit for bit in WALK_ACTIVE if not in_p[flip(vacuum, bit)]]
    print("  walk-active bits:", ", ".join(BIT_NAMES[bit] for bit in WALK_ACTIVE))
    print("  vacuum-exiting bits:", ", ".join(BIT_NAMES[bit] for bit in exits))
    check(len(exits) == 2 and len(WALK_ACTIVE) == 3, "vacuum cell has 2 exiting channels of 3 walk-active channels")
    return Fraction(len(exits), len(WALK_ACTIVE))


def q_concatenated(p_th: float, ratio: float, levels: int) -> float:
    return p_th * ratio ** (2**levels)


def q_by_recurrence(p_th: float, ratio: float, levels: int) -> float:
    q = p_th * ratio
    for _ in range(levels):
        q = q * q / p_th
    return q


def observed_rho_lambda() -> float:
    h = H0_KM_S_MPC / MPC_KM * HBAR_GEV_S
    return 3.0 * OMEGA_L * h * h * M_P_GEV * M_P_GEV / (8.0 * math.pi)


def print_gate_table(gates: list[Gate]) -> None:
    for gate in gates:
        print(f"  {gate.code:>2s}  {gate.status.value:11s}  {gate.name}")
        print(f"      {gate.finding}")


def main() -> None:
    print("ITEM 123 BOOT-SEQUENCE RESIDUAL-FAULT THEOREM ATTEMPT")

    print("\n[1] Distance-4 concatenation algebra")
    words = rm13_words()
    weights = sorted(sum(word) for word in words)
    d = code_distance(words)
    print(f"  [8,4,4] words={len(words)}; weights={weights}; distance={d}")
    check(len(words) == 16 and d == 4, "RM(1,3) supplies a distance-4 finite code")
    check((d - 1) // 2 == 1, "distance 4 corrects one lower-level fault; leading logical failure is a pair")
    for ell in range(0, 7):
        q_closed = q_concatenated(P_TH, 2.0 / 9.0, ell)
        q_recur = q_by_recurrence(P_TH, 2.0 / 9.0, ell)
        check(abs(q_closed / q_recur - 1.0) < 1e-12, f"recurrence solution q_{ell}=p_th*r^(2^{ell})")

    print("\n[2] Candidate physical subthreshold ratio")
    exit_fraction = derive_walk_active_exit_fraction()
    singlet_pair_factor = Fraction(1, 2)
    r_pair = singlet_pair_factor * exit_fraction * exit_fraction
    print(f"  r_pair = (1/2) * ({exit_fraction})^2 = {r_pair} = {float(r_pair):.12f}")
    check(r_pair == Fraction(2, 9), "pair-fault subthreshold ratio candidate is exactly 2/9")

    print("\n[3] Depth candidate and rho_Lambda")
    rho_obs = observed_rho_lambda()
    q_obs = rho_obs / (ALPHA0 * LAMBDA_QCD_GEV**4)
    ell_code = math.comb(4, 2)
    ell_q3_faces = 6
    ell_exact = math.log2(math.log(P_TH / q_obs) / math.log(1.0 / float(r_pair)))
    q_pred = q_concatenated(P_TH, float(r_pair), ell_code)
    rho_pred = q_pred * ALPHA0 * LAMBDA_QCD_GEV**4
    exact_r_for_ell6 = (q_obs / P_TH) ** (1.0 / (2**ell_code))
    pth_needed = q_obs / (float(r_pair) ** (2**ell_code))
    print(f"  depth from C(4,2) stabilizer-pair count = {ell_code}")
    print(f"  depth from Q3 coordinate faces           = {ell_q3_faces}")
    print(f"  exact depth required for r=2/9           = {ell_exact:.6f}")
    print(f"  q_obs = rho_obs/(alpha0 Lambda^4)        = {q_obs:.12e} = 2^-{math.log2(1/q_obs):.6f}")
    print(f"  q_pred = p_th*(2/9)^64                  = {q_pred:.12e} = 2^-{math.log2(1/q_pred):.6f}")
    print(f"  rho_pred                                = {rho_pred:.12e} GeV^4")
    print(f"  rho_obs                                 = {rho_obs:.12e} GeV^4")
    print(f"  rho_pred/rho_obs                        = {rho_pred/rho_obs:.6f}")
    print(f"  exact r at ell=6                        = {exact_r_for_ell6:.12f} ({(exact_r_for_ell6/float(r_pair)-1)*100:.3f}% above 2/9)")
    print(f"  p_th needed at r=2/9, ell=6             = {pth_needed:.12f} ({pth_needed/P_TH:.3f} x canonical p_th)")
    check(0.5 < rho_pred / rho_obs < 2.0, "the boot residual-fault route hits rho_Lambda to within a factor of two")
    check(abs(ell_exact - 6.0) < 0.01, "the observed target selects depth six once r=2/9 is granted")
    check(abs(exact_r_for_ell6 / float(r_pair) - 1.0) < 0.01, "the exact physical ratio is within 1% of 2/9")

    print("\n[4] Nearby-ledger nulls")
    candidates = {
        "2/9 pair probability": 2.0 / 9.0,
        "3/14 colour channel": 3.0 / 14.0,
        "1/4 rule class": 1.0 / 4.0,
        "1/3 single channel": 1.0 / 3.0,
        "1/2 bisection": 1.0 / 2.0,
        "sqrt(2/9) amplitude-as-probability": math.sqrt(2.0 / 9.0),
    }
    rows = []
    for name, ratio in candidates.items():
        for ell in range(4, 9):
            q = q_concatenated(P_TH, ratio, ell)
            rows.append((abs(math.log(q / q_obs)), name, ratio, ell, q / q_obs))
    rows.sort()
    for err, name, ratio, ell, rel in rows[:8]:
        print(f"  {name:36s} ell={ell} ratio={ratio:.9f} q/q_obs={rel:.3e} |ln|={err:.3f}")
    best = rows[0]
    check(best[1] in {"2/9 pair probability", "sqrt(2/9) amplitude-as-probability"}, "best nearby ledger is the 2/9 relation or its amplitude duplicate")
    print("  amplitude-as-probability is the same numerical relation with ell shifted by one;")
    print("  the probability form r=2/9, ell=6 is the category-clean reading.")

    print("\n[5] T1-T9 thermodynamic gates")
    gates = [
        Gate(
            "T1",
            "event algebra / jump support",
            Status.CONDITIONAL,
            "[8,4,4] distance and the single-cell walk-active exit algebra are explicit; the boot-level malignant-pair Kraus map is not yet built.",
        ),
        Gate(
            "T2",
            "Landauer equality versus bound",
            Status.PASS,
            "No Landauer equality sets the small number; q is a logical residual-fault probability and alpha0*Lambda^4 is the adopted erasure-density scale.",
        ),
        Gate(
            "T3",
            "mean residual current",
            Status.CONDITIONAL,
            "q_6 follows from the concatenation recurrence once r=2/9 and ell=6 are granted; the actual boot service current has not been simulated from the microscopic ledger.",
        ),
        Gate(
            "T4",
            "covariance / Fano factor",
            Status.OPEN,
            "The recurrence assumes independent malignant pairs and Poisson-like level faults; correlations would exponentiate through 2^ell and must be derived.",
        ),
        Gate(
            "T5",
            "correlation volume",
            Status.OPEN,
            "One residual-fault probability per cell-tick is mapped to a homogeneous vacuum density; the cross-cell correlation volume is still an assumption.",
        ),
        Gate(
            "T6",
            "observable map",
            Status.CONDITIONAL,
            "rho_Lambda=q alpha0 Lambda^4 gives a scalar vacuum-energy density, but the covariant stress-tensor map T_mu_nu=-rho g_mu_nu remains to be derived.",
        ),
        Gate(
            "T7",
            "scale accounting",
            Status.PASS,
            "The prediction uses Lambda_QCD, alpha0, p_th, r=2/9, and ell=6; it does not insert H0 or R_dS. Exact-target comparisons are audit-only.",
        ),
        Gate(
            "T8",
            "alternatives / nulls",
            Status.PASS,
            "Nearby ledgers are enumerated; the result is a factor-1.65 scale theorem, not an exact coefficient claim.",
        ),
        Gate(
            "T9",
            "promotion statement",
            Status.CONDITIONAL,
            "Promote only to conditional residual-fault scale theorem. Full derivation needs T1/T3/T4/T5/T6 closure and the factor-1.65 coefficient.",
        ),
    ]
    print_gate_table(gates)
    blockers = [g for g in gates if g.status is Status.OPEN]
    conditionals = [g for g in gates if g.status is Status.CONDITIONAL]
    print(f"  blockers={len(blockers)}; conditionals={len(conditionals)}")
    check(len(blockers) == 2, "two gates remain fully open")
    check(len(conditionals) == 4, "four gates are conditional")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  A real conditional theorem is available:")
    print("      q_ell = p_th r^(2^ell),  r = 2/9,  ell = 6")
    print(f"      q = {q_pred:.3e},  rho_Lambda = {rho_pred:.3e} GeV^4 = {rho_pred/rho_obs:.3f} rho_obs.")
    print("  This replaces the old horizon-input expression with a code-native residual-fault")
    print("  scale relation, but it does NOT yet derive the cosmological constant exactly.")
    print("  The unresolved content is precise: build the boot-level Kraus/current ledger,")
    print("  derive independence/covariance and correlation volume, derive the covariant")
    print("  vacuum stress tensor, and explain the remaining factor 1.65.")
    print("=" * 108)
    print("exit 0 -- conditional boot residual-fault scale theorem derived; exact rho_Lambda remains open.")


if __name__ == "__main__":
    main()
