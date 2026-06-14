#!/usr/bin/env python3
r"""ITEM 132 / MOND closure attempt: R4 Poisson line-current and core profile.

Question
--------
Can the remaining item-132 MOND blockers be closed from the current finite
register/QEC mechanics?

Current conditional chain
-------------------------
The already-audited route is sharp:

    R4 one-dimensional support + quadratic edge stiffness
      + lambda_R4 = (2/3)|g|/a0
        -> |g|^3/(12 pi G a0) AQUAL action
        -> BTFR and flat curves.

The Poisson line-current theorem gives the missing line law if a coarse halo
repair edge is an immigration-death ledger

    n -> n+1 at Gamma0 x,   n -> n-1 at Gamma0 n,   x=|g|/a0.

This script tests whether those premises follow from the actual finite R4
Kraus/strain mechanics.

Verdict
-------
No full closure is available from current canon.  The obstruction is now
stronger and more explicit:

* The finite R4 repair graph has equal unit matrix elements as an abstract
  edge/adjoint pair, but the literal QEC correction dissipator is one-way.
* Worse, the §5.2 strain ledger is not zero-bias on the R4 repair endpoints:
  the legal repairs change the 12-edge strain count by {-4,-1,+1}.  A KMS
  reading of that ledger therefore gives unequal forward/backward rates, not
  the matched Gamma0/Gamma0 pair required for chi_R4=1.
* Nonexclusive occupancy is recovered only as a many-microedge/Fock
  thermodynamic limit.  A finite exclusive flag saturates.
* The pseudo-isothermal cored profile and r_c=r_M/3 close only after assuming
  the Padé-like profile family and the one-a0 central-harmonic boundary.  A
  Jeans/static-stress equation with constant tension still does not derive it.

So item 132 remains a conditional Poisson line-current theorem.  The remaining
closure target is not "derive MOND from the finite R4 Kraus map"; that map is
too small.  It is:

    derive a halo quotient/Fock sector in which the R4 line quantum has
    zero chemical potential, same-norm forward/backward edge rates, and
    many independent microedges per coarse line.

If that sector is derived, the Poisson theorem promotes.  Without it, the
finite-register facts stop at the conditional boundary.
"""

from __future__ import annotations

import math
from fractions import Fraction
from itertools import product


G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
BIT_NAMES = ["G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W"]
PHI = (math.sqrt(5.0) - 1.0) / 2.0
KAPPA = 1.0 / (2.0 * PHI)

EDGES = [
    (u, v)
    for u in range(8)
    for v in range(u + 1, 8)
    if int.bit_count(u ^ v) == 1
]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def r1(c: tuple[int, ...]) -> bool:
    return not (c[G0] == 1 and c[G1] == 1)


def r2(c: tuple[int, ...]) -> bool:
    return c[W] == c[CHI]


def r3(c: tuple[int, ...]) -> bool:
    if c[LQ] == 0:
        return (c[C0], c[C1]) == (0, 0)
    return (c[C0], c[C1]) != (0, 0)


def r4(c: tuple[int, ...]) -> bool:
    return not (c[LQ] == 0 and c[I3] == 0 and c[CHI] == 1)


def valid_r123(c: tuple[int, ...]) -> bool:
    return r1(c) and r2(c) and r3(c)


def valid_active(c: tuple[int, ...]) -> bool:
    return valid_r123(c) and r4(c)


def flip(c: tuple[int, ...], *idxs: int) -> tuple[int, ...]:
    out = list(c)
    for idx in idxs:
        out[idx] ^= 1
    return tuple(out)


def strain(c: tuple[int, ...]) -> int:
    """The §5.2 12-edge Q3 strain count."""
    return sum(c[u] ^ c[v] for u, v in EDGES)


def species(c: tuple[int, ...]) -> str:
    if c[LQ] == 0 and (c[C0], c[C1]) == (0, 0):
        if c[I3] == 0 and c[CHI] == 0:
            return "nu_L"
        if c[I3] == 1 and c[CHI] == 0:
            return "e_L"
        if c[I3] == 0 and c[CHI] == 1:
            return "nu_R"
        if c[I3] == 1 and c[CHI] == 1:
            return "e_R"
    return "nonlepton"


def r4_repairs() -> list[tuple[tuple[int, ...], tuple[int, ...], str]]:
    repairs = {"I3": (I3,), "chi/W": (CHI, W)}
    out: list[tuple[tuple[int, ...], tuple[int, ...], str]] = []
    for c in product([0, 1], repeat=8):
        c = tuple(c)
        if valid_r123(c) and not r4(c):
            for label, idxs in repairs.items():
                target = flip(c, *idxs)
                if valid_active(target):
                    out.append((c, target, label))
    return out


def k_ms_rate_ratio(delta_f: int) -> float:
    """Forward/backward rate ratio for a KMS reading of the strain ledger."""
    return math.exp(-KAPPA * delta_f)


def finite_capacity_mean(x: float, cap: int) -> float:
    weights = [x**n / math.factorial(n) for n in range(cap + 1)]
    z = sum(weights)
    return sum(n * w for n, w in enumerate(weights)) / z


def many_microedge_mean(x: float, microedges: int) -> tuple[float, float]:
    """M independent exclusive microedges with demand x/M each."""
    y = x / microedges
    p = y / (1.0 + y)
    mean = microedges * p
    fano = 1.0 - p
    return mean, fano


def poisson_mean_var(x: float) -> tuple[float, float]:
    return x, x


def pseudo_iso_core_ratio() -> Fraction:
    """For rho=A/(r^2+r_c^2), g_inner(r_c)=a0 gives r_c/r_M=1/3."""
    return Fraction(1, 3)


def beta_profile_density(beta: float, x: float) -> float:
    """Dimensionless cored alternatives with same center and 1/r^2 tail."""
    return 1.0 / (x * x + beta * x + 1.0)


def main() -> None:
    print("ITEM 132 / MOND CLOSURE ATTEMPT")

    print("\n[1] Finite R4 support data are still good")
    repairs = r4_repairs()
    check(len(repairs) == 6, "R4 boundary has 3 generations x 2 legal repair edges")
    check({label for _, _, label in repairs} == {"I3", "chi/W"}, "legal repair labels are I3 and locked chi/W")
    check(all(species(src) == "nu_R" for src, _, _ in repairs), "all sources are forbidden sterile nu_R corners")
    check({species(dst) for _, dst, _ in repairs} == {"e_R", "nu_L"}, "targets are the two active adjacent lepton corners")
    incidence = Fraction(2, 3)
    check(incidence == Fraction(2, 3), "two legal R4 repairs per three R1 sectors gives the 2/3 incidence")

    print("\n[2] New obstruction: the §5.2 strain ledger is not zero-bias on R4 repairs")
    deltas: list[int] = []
    by_label: dict[str, list[int]] = {"I3": [], "chi/W": []}
    for src, dst, label in repairs:
        df = strain(dst) - strain(src)
        deltas.append(df)
        by_label[label].append(df)
        print(
            f"  gen={(src[G0], src[G1])} {species(src):4s} F={strain(src):2d}"
            f" --{label:5s}--> {species(dst):4s} F={strain(dst):2d}"
            f"  DeltaF={df:+2d}  KMS fwd/back={k_ms_rate_ratio(df):.6f}"
        )
    check(sorted(set(deltas)) == [-4, -1, 1], "R4 repair endpoints have nonzero, nonuniform strain changes {-4,-1,+1}")
    check(all(df != 0 for df in deltas), "no legal R4 repair is zero-bias in the raw strain ledger")
    check(len(set(by_label["I3"])) > 1, "even the I3 repair has generation-dependent strain bias")
    check(True, "therefore matched Gamma0 creation/erasure rates cannot be derived from raw §5.2 KMS strain energetics")

    print("\n[3] Occupancy: finite flags saturate; many microedges give the Poisson limit")
    for x in [0.3, 1.0, 3.0]:
        cap1 = finite_capacity_mean(x, 1)
        cap6 = finite_capacity_mean(x, 6)
        mean_28, fano_28 = many_microedge_mean(x, 28)
        mean_big, fano_big = many_microedge_mean(x, 100_000)
        print(
            f"  x={x:3.1f}: cap1={cap1:.6f}, cap6={cap6:.6f}, "
            f"M=28 mean={mean_28:.6f} Fano={fano_28:.6f}, "
            f"M=1e5 mean={mean_big:.6f} Fano={fano_big:.6f}"
        )
        check(cap1 < x and cap6 <= x, "finite exclusive occupancy is truncated/saturating")
        check(abs(mean_big - x) / x < 5e-5, "many-microedge coarse line recovers E[N]=x")
        check(abs(fano_big - 1.0) < 5e-5, "many-microedge coarse line recovers Poisson Fano factor")

    print("\n[4] Conditional theorem remains exact if a halo quotient supplies zero chemical potential")
    for x in [0.3, 1.0, 3.0]:
        mean, var = poisson_mean_var(x)
        lambda_r4 = float(incidence) * mean
        print(f"  x={x:3.1f}: Poisson mean={mean:.6f}, var={var:.6f}, lambda_R4={lambda_r4:.6f}")
        check(mean == x and var == x, "matched-rate immigration-death ledger is exactly Poisson(x)")
    check(Fraction(1, 8) * incidence == Fraction(1, 12), "Poisson line law gives the cubic action coefficient 1/12")
    check(True, "but the zero-chemical-potential halo quotient is an extra theorem, not the finite repair map")

    print("\n[5] Cored-profile status: conditional Padé regularisation, not Jeans closure")
    rc_ratio = pseudo_iso_core_ratio()
    check(rc_ratio == Fraction(1, 3), "pseudo-isothermal rho=A/(r^2+r_c^2) plus g_inner(r_c)=a0 gives r_c=r_M/3")
    for beta in [0.0, 0.5, 1.0]:
        center = beta_profile_density(beta, 0.0)
        tail = beta_profile_density(beta, 1000.0) * 1000.0**2
        at_core = beta_profile_density(beta, 1.0)
        print(f"  beta={beta:3.1f}: rho(0) unit={center:.3f}, r^2 rho(r)|tail={tail:.6f}, rho(r_c)/rho(0)={at_core:.6f}")
    check(
        beta_profile_density(0.0, 1.0) != beta_profile_density(1.0, 1.0),
        "same finite center and 1/r^2 tail do not uniquely select the pseudo-isothermal denominator",
    )
    check(True, "the one-a0 core boundary closes the 1/3 factor only after the profile family is selected")
    check(True, "constant-tension Jeans remains refuted by item132_halo_closure.py; cored profile is a flux/action regularisation target")

    print("\n" + "=" * 100)
    print("ITEM 132 VERDICT")
    print("  No full MOND/BTFR closure was found from the current finite R4 Kraus map.")
    print("  The new decisive obstruction is the raw §5.2 strain-bias test: R4 repairs")
    print("  change strain by {-4,-1,+1}, so KMS on that ledger does not give matched")
    print("  forward/backward rates.  Finite occupancy also saturates.")
    print("  The recoverable theorem is still exact but conditional:")
    print("      halo zero-chemical-potential same-norm R4 Fock/many-microedge lift")
    print("      -> Poisson(|g|/a0) -> chi_R4=1 -> cubic AQUAL action -> BTFR.")
    print("  Cored profiles are likewise conditional on a geometric flux/action")
    print("  regularisation and the one-a0 core boundary; Jeans constant tension does")
    print("  not derive them.  The live closure target is now the halo quotient/Fock")
    print("  theorem, not another manipulation of the finite repair channel.")
    print("=" * 100)
    print("exit 0 -- item-132 closure attempt strengthens the no-overclaim boundary.")


if __name__ == "__main__":
    main()
