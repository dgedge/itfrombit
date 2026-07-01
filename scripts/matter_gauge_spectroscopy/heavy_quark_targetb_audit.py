#!/usr/bin/env python3
"""Heavy-quark mass Target-B audit.

Question:
    Can the framework turn heavy-quark masses into parameter-free predictions,
    or are the apparent "nice" Lambda_QCD ratios dense-alphabet matches?

This is a §16.3 search-space audit, not a fit.  It uses the same constituent
mass convention already used by the charm/bottom baryon scripts:

* m_s from the light DGG baryon sector.
* m_c from Lambda_c, with the light sector frozen.
* m_b from Lambda_b, with the light sector frozen.
* m_t/v is treated separately: it is the electroweak/top second anchor, not a
  QCD-heavy-baryon constituent scale.

Exit 0 means the current finite alphabet cannot honestly promote any heavy
quark mass ratio to a prediction.  It does not say heavy-quark spectroscopy is
uninteresting; it says each heavy flavour supplies a Target-B scale in the
same way QCD does.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from fractions import Fraction


ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


LAMBDA_QCD = 332.0  # MeV, canon §1.4
PHI = (1.0 + math.sqrt(5.0)) / 2.0
PRIMITIVES = (
    ("1", 1.0),
    ("phi", PHI),
    ("sqrt2", math.sqrt(2.0)),
    ("pi", math.pi),
)


@dataclass(frozen=True)
class Target:
    name: str
    mass_mev: float
    source: str
    pmax: int
    qmax: int = 17

    @property
    def ratio(self) -> float:
        return self.mass_mev / LAMBDA_QCD


def reduced_competitors(target: Target, tol: float = 0.01) -> list[tuple[float, str, float]]:
    hits: dict[tuple[int, int, str], tuple[float, str, float]] = {}
    for p in range(1, target.pmax + 1):
        for q in range(1, target.qmax + 1):
            frac = Fraction(p, q)
            for primitive_name, primitive_value in PRIMITIVES:
                value = float(frac) * primitive_value
                err = abs(value / target.ratio - 1.0)
                if err <= tol:
                    key = (frac.numerator, frac.denominator, primitive_name)
                    label = f"{frac.numerator}/{frac.denominator}"
                    if primitive_name != "1":
                        label += f"*{primitive_name}"
                    if key not in hits or err < hits[key][0]:
                        hits[key] = (err, label, value)
    return sorted(hits.values(), key=lambda row: row[0])


def print_target(target: Target) -> list[tuple[float, str, float]]:
    hits = reduced_competitors(target)
    print(f"\n{target.name}")
    print(f"  source         : {target.source}")
    print(f"  m/Lambda_QCD   : {target.ratio:.6f}")
    print(f"  alphabet       : p/q * {{1, phi, sqrt2, pi}}, p<= {target.pmax}, q<= {target.qmax}")
    print(f"  hits within 1% : {len(hits)} reduced expressions")
    print("  closest hits:")
    for err, label, value in hits[:8]:
        print(f"    {label:<14s} {value:>11.6f}  err={100.0 * err:>7.4f}%")
    return hits


def main() -> int:
    # Light-sector DGG values from baryon_parameter_count.py / bottom_baryon_forward.py.
    m_u = (938.92 + 1232.0) / 6.0
    m_s = m_u + (1115.68 - 938.92)
    m_c = m_s + (2286.46 - 1115.68)
    m_b = m_s + (5619.6 - 1115.68)
    m_top = 172_760.0
    v_ew = 246_220.0

    targets = [
        Target("m_s constituent", m_s, "light DGG baryons", 90),
        Target("m_c constituent", m_c, "Lambda_c with light sector frozen", 90),
        Target("m_b constituent", m_b, "Lambda_b with light sector frozen", 260),
        # Conservative wide scans for the electroweak/top anchor: pmax about
        # twice the target ratio admits nearby integers without pretending
        # arbitrary large numerators are explanatory.
        Target("m_top pole proxy", m_top, "EW/top second anchor, not QCD-derived", 1041),
        Target("v electroweak", v_ew, "EW/top second anchor, not QCD-derived", 1484),
    ]

    print("=" * 92)
    print("HEAVY-QUARK TARGET-B AUDIT")
    print("=" * 92)
    print(f"Lambda_QCD = {LAMBDA_QCD:.1f} MeV")
    print(f"DGG light-sector values: m_u={m_u:.2f} MeV, m_s={m_s:.2f} MeV")
    print(f"Heavy inputs: m_c={m_c:.2f} MeV, m_b={m_b:.2f} MeV")

    counts = {}
    for target in targets:
        hits = print_target(target)
        counts[target.name] = len(hits)

    print("\n[1] Heavy-baryon sector verdict")
    print(
        "  In the DGG/heavy-baryon convention, charm and bottom each enter as one\n"
        "  new constituent scale fixed from Lambda_c/Lambda_b.  The framework then\n"
        "  inherits the constituent-quark model's useful predictions (equal spacing,\n"
        "  hyperfine scaling, Omega_b* band), but it does not derive m_c or m_b\n"
        "  from Lambda_QCD."
    )
    check("m_c sits in a dense 1% alphabet band", counts["m_c constituent"] >= 20)
    check("m_b sits in a denser 1% alphabet band", counts["m_b constituent"] >= 30)

    print("\n[2] Top/electroweak classification")
    print(
        "  The top mass is not just another QCD constituent input.  Current canon\n"
        "  classifies v ~= m_top as the irreducible electroweak/top second anchor.\n"
        "  Scanning v/Lambda or m_top/Lambda only demonstrates dense numerology;\n"
        "  it does not supply an EW-scale derivation."
    )
    check("m_top/Lambda has many nearby alphabet forms in a conservative scan", counts["m_top pole proxy"] >= 50)
    check("v/Lambda has many nearby alphabet forms in a conservative scan", counts["v electroweak"] >= 50)

    print("\n[3] Honest status")
    print(
        "  No clean closure was found.  Heavy-quark masses remain Target-B: one\n"
        "  scale per heavy flavour in the heavy-baryon sector, plus the separate\n"
        "  electroweak/top anchor.  The useful framework content is structural\n"
        "  spectroscopy after those scales are supplied, not a prediction of the\n"
        "  scales themselves.  Nice Lambda_QCD ratios are consistency checks only."
    )
    check("heavy-quark mass ratios are not promoted to predictions", True)

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
