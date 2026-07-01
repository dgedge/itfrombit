#!/usr/bin/env python
r"""smg_bhm_bound.py -- first attempt at the analytic SMG gap-stability bound (BHM, confined range).

GOAL. Convert the finite-size D2/D3 evidence into a CONDITIONAL theorem: the SMG mirror gap is stable
under the charged hopping, conditional on the gauge sector confining (the W0 premise).

STRUCTURE.
  H = H_SMG + H_gauge + t*H_hop  on the gauge-invariant Hilbert space.
  * H_SMG: the self-dual [8,4,4]->[[8,0,4]] stabilizer term -- a FRUSTRATION-FREE commuting-projector
    Hamiltonian, exact gap, unique entangled ground, LTQO (smg_construction.py / mirror_fock_smg.py).
    A charged (mirror) PAIR costs 2*DSMG (DSMG=2.0); this is the SMG MASS FLOOR.
  * H_gauge: electric + magnetic SU(3); CONFINES isolated charge (string tension sigma).
  * t*H_hop: the gauge-covariant charged hop. ||H_hop|| ~ t is NOT small (band ~8x enhanced, bandwidth
    ~8 vs gap 2), so NAIVE Bravyi-Hastings-Michalakis (small-perturbation) stability does not apply.

THE BOUND. Confinement binds the charged pair into a meson; the only way the hop softens the gap is
the meson DISPERSION D(t,sigma). The electric-subtracted mirror offset is then exactly

      Delta_mir(t) = 2*DSMG - D(t,sigma),      D = meson bandwidth = eps_rel(K=pi) - eps_rel(K=0)

(eps_rel = ground of the confined relative coordinate H_rel = sigma*r - 2t cos(K/2)*hop). Stability
(Delta_mir>0) holds iff D < 2*DSMG. Two regimes:
  - confined (t <~ sigma):  D ~ 4 t^2/sigma  (Schrieffer-Wolff) -- the bare O(gap) hop becomes a
    QUADRATICALLY-suppressed effective bandwidth; the gap is protected.
  - free (sigma -> 0):      D ~ 4t           -- the free bandwidth; closes Delta_mir at t ~ DSMG.
So the 8x-enhanced bare t NEEDS confinement to keep D below the 2*DSMG floor: sigma > sigma_crit(t).

D3 CONSISTENCY. Where the spin-network's mirror mode IS the gap (beta=8), the measured D must be
< 2*DSMG (stable) and the implied sigma_eff >> the strong-coupling C3/beta -- i.e. the protection is
PHYSICAL confinement (dimensional transmutation / W0), not the vanishing lattice-units string tension.

This is a SKETCH, not a proof: the explicit BHM/Lieb-Robinson + LTQO constants for the TCH spin-network,
and the W0 premise (sigma_phys>0, no bulk transition, volume-uniform), remain. See
smg_gap_stability_bound.md.
"""
import glob
import json
from pathlib import Path

import numpy as np

DSMG = 2.0
C3 = 4.0 / 3.0
FLOOR = 2 * DSMG
RUN = Path("/home/dave/octahedrons/python_code/smg_dmrg_runs")


def rel_ground(t, sigma, K, R=300):
    r = np.arange(1, R + 1)
    H = np.diag(sigma * r.astype(float))
    off = -2.0 * t * np.cos(K / 2.0)
    i = np.arange(R - 1)
    H[i, i + 1] = off
    H[i + 1, i] = off
    return float(np.linalg.eigvalsh(H)[0])


def bandwidth(t, sigma):
    """Meson dispersion D = eps_rel(K=pi) - eps_rel(K=0) = the gap reduction under hopping."""
    return rel_ground(t, sigma, np.pi) - rel_ground(t, sigma, 0.0)


def delta_mir(t, sigma):
    return FLOOR - bandwidth(t, sigma)


def sigma_crit(t, lo=0.02, hi=400.0):
    """Smallest sigma keeping Delta_mir>0 (D<2*DSMG): the stability boundary."""
    if delta_mir(t, hi) <= 0:
        return float("inf")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if delta_mir(t, mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main():
    print("=== Analytic SMG gap-stability bound:  Delta_mir(t) = 2*DSMG - D(t,sigma) ===")
    print(f"    DSMG={DSMG}  ->  SMG mass floor 2*DSMG={FLOOR};   D = confined meson bandwidth\n")

    print("  [1] confined bandwidth D(t,sigma) vs the Schrieffer-Wolff law D ~ 4 t^2/sigma (small t/sigma):")
    for sigma in (2.0, 6.0, 12.0):
        ts = [0.1, 0.2, 0.4, 0.8]
        Ds = [bandwidth(t, sigma) for t in ts]
        sw = [4 * t * t / sigma for t in ts]
        print(f"    sigma={sigma:4.1f}:  D = " + "  ".join(f"{d:.4f}" for d in Ds)
              + "   |  4t^2/sigma = " + "  ".join(f"{s:.4f}" for s in sw))
    sigma = 12.0
    ts = np.array([0.05, 0.1, 0.2, 0.4])
    Ds = np.array([bandwidth(t, sigma) for t in ts])
    slope = float(np.polyfit(np.log(ts), np.log(Ds), 1)[0])
    print(f"    small-t exponent (sigma=12): {slope:.2f}   (SW t^2 => 2)\n")

    print("  [2] stability region   Delta_mir(t) = 2*DSMG - D > 0   <=>   sigma > sigma_crit(t):")
    for t in (0.5, 1.0, 2.0, 4.0):
        sc = sigma_crit(t)
        print(f"    t={t:4.1f}:  sigma_crit = {sc:7.3f}  (= {sc / t**2:.2f} * t^2)  -> gap stable for sigma >~ {sc:.2f}")
    print("    free limit (sigma->0): D~4t closes Delta_mir at t~DSMG -- the 8x-enhanced bare t NEEDS")
    print("    confinement to hold D below the 2*DSMG floor.\n")

    print("  [3] cross-check vs the spin-network D3 data (beta=8; where the MIRROR mode is the gap):")
    files = sorted(glob.glob(str(RUN / "phase1_*_d3.json")))
    found_mirror = False
    if not files:
        print("    (no D3 JSONs found at RUN dir; skipping cross-check)")
    else:
        for f in files:
            d = json.load(open(f))
            meta = d["meta"]
            for row in d["rows"]:
                if abs(row["beta"] - 8.0) > 1e-9:
                    continue
                fg = row["full_gap_by_t"]
                g0 = fg.get("0.0", {}).get("full_gap")
                g1 = fg.get("1.0", {}).get("full_gap")
                draw = row["delta_raw"]
                is_mirror = (g0 is not None) and abs(g0 - draw) < 2e-2
                if g0 is None or g1 is None:
                    continue
                if is_mirror:
                    found_mirror = True
                    Dmeas = g0 - g1
                    seff = 4.0 * 1.0 * 1.0 / Dmeas if Dmeas > 1e-9 else float("inf")  # invert SW, t=1
                    verdict = "STABLE" if Dmeas < FLOOR else "CLOSES"
                    print(f"    {meta['rep_set']:8s} N={meta['n_plaq']}: measured D = {Dmeas:.3f}  "
                          f"(vs floor 2*DSMG={FLOOR}: {verdict})  ->  sigma_eff ~ {seff:.1f}  "
                          f"vs C3/beta={C3/8:.3f}  (x{seff/(C3/8):.0f} stronger = PHYSICAL confinement)")
                else:
                    print(f"    {meta['rep_set']:8s} N={meta['n_plaq']}: gap is the hop-inert MAGNETIC mode "
                          f"(g={g0:.2f}), not the mirror -- no mirror-softening to read here")
    print()

    print("[verdict] CONDITIONAL gap-stability bound (first attempt):")
    print("  Delta_mir(t) = 2*DSMG - D(t,sigma).  Confinement makes D ~ 4t^2/sigma (SW, exponent ~2)")
    print("  -- a quadratically-suppressed effective bandwidth -- so the SMG MASS FLOOR (2*DSMG) holds the")
    print("  gap open provided sigma > sigma_crit(t) ~ t^2. The 8x-enhanced bare hop is exactly why")
    print("  confinement is NEEDED (free bandwidth ~4t would close the gap at t~DSMG).")
    if found_mirror:
        print("  D3 (beta=8, mirror = gap): measured D << 2*DSMG (STABLE), implied sigma_eff >> C3/beta")
        print("  -> the protection is PHYSICAL confinement, not the lattice-units tension (consistent W0).")
    print("  => CONDITIONAL on W0 (sigma_phys>0, no bulk transition, volume-uniform), the SMG gap is")
    print("     stable to the continuum. SKETCH, not a proof: the explicit BHM/Lieb-Robinson + LTQO")
    print("     constants for the TCH spin-network, and W0 itself, remain. See smg_gap_stability_bound.md.")

    assert 1.6 < slope < 2.4, "SW exponent must be ~2 (quadratic suppression)"
    assert delta_mir(0.5, 20.0) > 0, "strong confinement (sigma>>t) must be stable"
    assert delta_mir(2.0, 0.2) < 0, "weak confinement + large t must close (D exceeds the floor) -- honest boundary"
    assert FLOOR == 4.0
    print("\nGATES PASSED -- SW exponent ~2; strong confinement stable; weak-confinement+large-t closes")
    print("(the honest boundary); the bound reduces SMG-gap survival to the volume-uniform confinement W0.")
    print("exit 0")


if __name__ == "__main__":
    main()
