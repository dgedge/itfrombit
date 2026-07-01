#!/usr/bin/env python3
r"""Black-hole flux coefficient: the all-contact (Moore) alphabet is DATA-SELECTED -> a bounded local theorem.

B6 (bh_all_contact_severing_theorem_audit.py) showed the V_cell/V_Sch ISOMETRY does not FORCE the
all-contact (full Moore) service alphabet: O_h cubic symmetry + connectedness leave three nearest-
contact alphabets {F+E, E+C, F+E+C}, only the last giving 10/27. This script asks the OTHER question
-- not "does the isometry force it?" but "does the FLUX NORMALIZATION select it?" -- and finds it does,
uniquely.

Target: the Stefan-Hawking coefficient P = M_P^4/(15360 pi M^2), which is the naive
Stefan coefficient for a two-helicity massless bosonic channel.  The beta-one
freeze shell + escape-cone reduction (bh_flux_coefficient_gate.py) turns this into a
required local attempt rate Gamma_req/Lambda. For each symmetry-allowed alphabet the
attempt-rate prefactor is the outward-emission fraction (emit+latch)/(slots+latch);
P/P_SB scales linearly with it.

Result: among the O_h-invariant connected alphabets x {latch, no-latch}, ONLY all-contact
Moore + the V_cell latch (10/27) reproduces the two-helicity Stefan target (0.29%);
every alternative misses by >= 6.8%. So 10/27 is DATA-SELECTED -- the
predicted-structure + data-selection pattern used elsewhere -- and is a BOUNDED LOCAL
coefficient (the 27-slot Moore neighbourhood + latch, no global horizon graph).

Honest residual: the target is a two-helicity Stefan-channel normalization, not a full
species-summed Hawking luminosity.  The real flux also needs the QEC species/polarization
emission ledger plus the already-computed spin/partial-wave greybody transfer
(bh_greybody_transfer.py).  Bulk-service-universality (a horizon cell runs the same
isotropic bulk Landauer-Moore service as any lattice cell, the radial freezing selecting
only WHICH serviced slots emit) is the physical ground the isometry does not force.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bh_flux_coefficient_gate as gate          # noqa: E402  (tested gate: required rate, alpha0)
import bh_all_contact_severing_theorem_audit as audit  # noqa: E402  (Moore orbits, connectedness, prefactor)

GAMMA_REQ = gate.required_gamma_for_stefan(1.0, outward_hemisphere=True)   # /Lambda, beta-1 target
A0 = gate.ALPHA0


def p_over_psb(emit: int, total: int) -> float:
    """P/P_SB for attempt-rate prefactor (emit/total): linear in the prefactor."""
    return (emit / total) * A0 / GAMMA_REQ


def main() -> None:
    print("BLACK-HOLE FLUX COEFFICIENT — all-contact alphabet DATA-SELECTION (bounded local theorem)")
    print("=" * 100)

    counts = Counter(audit.orbit_name(v) for v in audit.moore_shell())
    assert counts == {"F": 6, "E": 12, "C": 8}, "Moore shell O_h orbits face/edge/corner = 6/12/8"
    print(f"  required local attempt rate (beta-1 shell, two-helicity Stefan target): Gamma_req/Lambda = {GAMMA_REQ:.12f}")
    print(f"  Moore shell O_h orbits: F=6  E=12  C=8")

    connected = [names for names in audit.all_orbit_subsets() if audit.connected(audit.orbit_union(names))]
    labels = sorted("".join(sorted(n)) for n in connected)
    print(f"  O_h-invariant CONNECTED nearest-contact alphabets: {labels}  (cubic symmetry + connectedness)")

    rows = []
    for names in connected:
        for latch in (True, False):
            emit, total = audit.emitted_prefactor(names, latch=latch)
            rows.append(("".join(sorted(names)), latch, emit, total, p_over_psb(emit, total)))
    rows.sort(key=lambda x: abs(x[4] - 1.0))

    print(f"\n  {'alphabet':<8s} {'latch':<8s} {'emit/total':>12s} {'frac':>7s} {'P/P_SB':>8s}")
    print("  " + "-" * 52)
    for label, latch, emit, total, r in rows:
        tag = "  <- DATA-SELECTED" if (label == "CEF" and latch) else ""
        print(f"  {label:<8s} {'+latch' if latch else 'no-latch':<8s} {f'{emit}/{total}':>12s} "
              f"{emit/total:>7.4f} {r:>8.4f}{tag}")

    winner = next(r for r in rows if r[0] == "CEF" and r[1] is True)
    others = [r for r in rows if not (r[0] == "CEF" and r[1] is True)]
    margin = min(abs(r[4] - 1.0) for r in others)

    assert abs(winner[2] / winner[3] - 10.0 / 27.0) < 1e-12, "all-contact Moore + latch = 10/27"
    assert abs(winner[4] - 1.0) < 0.01, "10/27 reproduces two-helicity Stefan target to <1%"
    assert margin > 0.05, f"every one of the other 5 symmetry-allowed candidates misses P_SB by >5% (closest {margin*100:.1f}%)"

    print(f"\n  UNIQUE: 10/27 lands at {abs(winner[4]-1)*100:.2f}% from the two-helicity Stefan target; next-closest alternative is {margin*100:.1f}% off.")
    print(f"""
{"=" * 100}
VERDICT (bounded local theorem, exit 0):  the all-contact Moore alphabet is DATA-SELECTED; 10/27 is a
bounded local flux coefficient.

  B6 stands at the ISOMETRY level: V_cell/V_Sch do not force the contact alphabet. But the FLUX
  NORMALIZATION does -- among the O_h-invariant connected nearest-contact alphabets the cubic symmetry
  allows, {{F+E: 6/19, E+C: 3/7, F+E+C: 10/27}} x {{latch, no-latch}}, ONLY all-contact Moore WITH the
  V_cell latch (10/27) reproduces the two-helicity Stefan-channel coefficient ({abs(winner[4]-1)*100:.2f}%);
  every one of the other five candidates misses by >= {margin*100:.1f}%. 10/27 thus moves from
  "isometry-unforced premise" (B6) to "DATA-SELECTED coefficient" -- the predicted-structure +
  data-selection pattern the framework uses elsewhere (cf. the sterile-release billing power).

  BOUNDED LOCAL: the coefficient is fixed entirely by the local 26-slot Moore neighbourhood + the one
  V_cell latch -- a bounded rational (9+1)/(26+1), with NO global horizon graph. (Contrast: fast
  scrambling is provably NON-local, the O(1)-gap expander of B5. Flux = local stencil; scrambling =
  global graph -- the two live BH maps are cleanly separated by locality.)

  PHYSICAL GROUND (the premise the isometry does not force): bulk-service-universality -- a horizon
  cell runs the SAME isotropic bulk Landauer-Moore service (item 120, the 26-slot shell) as any lattice
  cell; the radial freezing selects only WHICH serviced slots EMIT (the outward hemisphere, 9 of 26 +
  latch), not the service stencil. This separates the two horizon numbers consistently: the Bekenstein
  ENTROPY rides on radial PAIR severing (55/8 alpha0^2, derived); the Hawking FLUX rides on the
  isotropic SERVICE stencil (10/27 on alpha0) -- different structures, different alpha0-orders, no conflict.

  HONEST RESIDUAL: (a) target is the two-helicity Stefan-channel normalization, not the full
  species-summed Hawking luminosity; the real flux also needs the QEC species/polarization ledger
  and the spin/partial-wave greybody transfer (bh_greybody_transfer.py). (b) bulk-service-universality
  is argued from substrate uniformity, not forced by V_cell; the data-selection is independent
  evidence FOR it, but it remains the named premise.
{"=" * 100}""")
    print(f"exit 0 -- BH flux: all-contact Moore alphabet DATA-SELECTED (10/27 -> P/P_SB {winner[4]:.4f}; "
          f"5 symmetry-allowed alternatives miss >= {margin*100:.0f}%); bounded local; residual = species ledger + bulk-service premise.")


if __name__ == "__main__":
    main()
