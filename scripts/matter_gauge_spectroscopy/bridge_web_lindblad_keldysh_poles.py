#!/usr/bin/env python3
"""Bridge-plus-web Lindblad/Keldysh pole audit for dressed alpha.

This is the direct "do not introduce gamma_mon" test.

We integrate out the Wilson photon web in the zero-temperature, weak-coupling
Davies/Keldysh limit using the same pinned ingredients as the bridge/web escape
calculation:

  * Peierls current portal on the two bridge links;
  * eta = -1, the pi-flux branch;
  * E_ref = pair-band bottom;
  * g = 1/3, the Grover gauge-matter portal.

The resulting bridge-block Lindblad poles are then read directly.  If the web
bath itself supplies the missing monitored-bridge relaxation scale, its slow
pole should be O(2.17e-1 Lambda), the gap required by the dressed-alpha response
with the computed Gamma_esc.  If the slow pole is zero or O(alpha), the scalar
gamma_mon has not been derived; a separate substrate monitor remains necessary.
"""
from __future__ import annotations

import numpy as np

import dressed_alpha_bridge_web_open_system as bw


T_M = 1.0 / 3.0
GROUP_TOL = 1e-10
RATE_TOL = 1e-14


def eigen_groups(evals: np.ndarray, tol: float = GROUP_TOL) -> list[np.ndarray]:
    """Group numerically degenerate eigenvalues."""
    groups: list[np.ndarray] = []
    used = np.zeros(len(evals), dtype=bool)
    for k, value in enumerate(evals):
        if used[k]:
            continue
        group = np.where(np.abs(evals - value) < tol)[0]
        used[group] = True
        groups.append(group)
    return groups


def weighted_quantile(values: np.ndarray, weights: np.ndarray, p: float) -> float:
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    return float(np.interp(p, cdf, v))


def grouped_web_poles(h2: np.ndarray, idx: dict[tuple[int, int], int], s_eta):
    """Return the physical bright poles after integrating out the web.

    In a degenerate bridge eigenspace the Peierls current creates one bright
    vector and the orthogonal vectors are dark.  Grouping the eigenspaces avoids
    assigning basis-dependent rates inside accidental degeneracies.
    """
    evals_unit, evecs = np.linalg.eigh(h2)
    evals = T_M * evals_unit
    e_ref = float(evals.min())
    omega = evals - e_ref
    current = bw.current_portal(idx, bw.ETA_PIN)
    amps = evecs.T @ current

    groups = eigen_groups(evals)
    rows = []
    dark_dim = 0
    for group in groups:
        weight = float(np.sum(amps[group] ** 2))
        group_omega = float(np.mean(omega[group]))
        group_rate = (
            2.0
            * np.pi
            * bw.ALPHA0
            * bw.G_PORTAL
            * bw.G_PORTAL
            * weight
            * s_eta(group_omega, bw.ETA_PIN)
        )
        if weight > RATE_TOL:
            dark_dim += len(group) - 1
        else:
            dark_dim += len(group)
        rows.append(
            {
                "size": len(group),
                "omega": group_omega,
                "weight": weight,
                "rate": float(group_rate),
            }
        )

    rates = np.array([r["rate"] for r in rows])
    weights = np.array([r["weight"] for r in rows])
    positive = rates[rates > RATE_TOL]
    return rows, rates, weights, dark_dim


def main() -> None:
    h2, _pairs, idx, _bas = bw.build_pair_system()
    s_eta, omega_max = bw.photon_form_factor()
    c_j, gap_unit = bw.response_coefficient(h2, idx)
    rows, rates, weights, dark_dim = grouped_web_poles(h2, idx, s_eta)

    gamma_escape = float(rates.sum())
    target_gap = c_j * gamma_escape / bw.DELTA_TARGET
    positive = rates[rates > RATE_TOL]
    slow_bright_population_pole = float(positive.min()) if len(positive) else 0.0
    slow_bright_coherence_pole = 0.5 * slow_bright_population_pole
    max_bright_pole = float(positive.max()) if len(positive) else 0.0
    q05 = weighted_quantile(rates, weights, 0.05)
    q50 = weighted_quantile(rates, weights, 0.50)
    q90 = weighted_quantile(rates, weights, 0.90)

    print("BRIDGE + WILSON WEB LINDBLAD/KELDYSH POLES")
    print("=" * 86)
    print("Pinned inputs:")
    print(f"  t_m/t_web                 = {T_M:.6f}")
    print(f"  eta                       = {bw.ETA_PIN:+d}")
    print(f"  g_portal                  = {bw.G_PORTAL:.6f}")
    print(f"  photon omega_max          = {omega_max:.6f} Lambda")
    print(f"  response c_J              = {c_j:.6f}")
    print(f"  unit population gap       = {gap_unit:.6f} Lambda")
    print()
    print("Integrated-web output:")
    print(f"  bright eigenspaces        = {int(np.sum(rates > RATE_TOL))}")
    print(f"  bridge dark dimension     = {dark_dim} of {len(idx)}")
    print(f"  Gamma_esc = sum bright    = {gamma_escape:.6e} Lambda")
    print(f"  required response gap     = {target_gap:.6e} Lambda")
    print()
    print("Actual bridge-block Liouvillian decay poles:")
    print("  full bridge slow pole     = 0 exactly (dark bridge subspace)")
    print(f"  slow bright population    = {slow_bright_population_pole:.6e} Lambda")
    print(f"  slow bright coherence     = {slow_bright_coherence_pole:.6e} Lambda")
    print(f"  strongest bright pole     = {max_bright_pole:.6e} Lambda")
    print(f"  current-weighted q05      = {q05:.6e} Lambda")
    print(f"  current-weighted median   = {q50:.6e} Lambda")
    print(f"  current-weighted q90      = {q90:.6e} Lambda")
    print()
    print("Ratios to the required response gap:")
    print(f"  Gamma_esc / gap_required  = {gamma_escape / target_gap:.6e}")
    print(f"  strongest / gap_required  = {max_bright_pole / target_gap:.6e}")
    print(f"  median / gap_required     = {q50 / target_gap:.6e}")
    print()
    print("Dominant Keldysh/FGR bright poles:")
    print(f"  {'omega':>10} {'weight':>12} {'rate':>13} {'rate/target':>13}")
    for row in sorted(rows, key=lambda r: r["rate"], reverse=True)[:8]:
        print(
            f"  {row['omega']:>10.6f} {row['weight']:>12.6f}"
            f" {row['rate']:>13.6e} {row['rate'] / target_gap:>13.6e}"
        )

    print()
    print("VERDICT")
    print("  The Wilson web escape bath does not generate the missing bridge relaxation")
    print("  gap.  After integrating out the web, the exact slow Liouvillian pole of")
    print("  the bridge block is zero because the Peierls-current escape is low rank")
    print("  and leaves a dark bridge subspace.  Even the brightest web-induced poles")
    print("  are O(alpha) escape widths, hundreds of times below the response gap")
    print("  needed to reproduce the dressed-alpha shift.  Thus the full web-only")
    print("  pole calculation does not derive gamma_mon; it shows that a separate")
    print("  substrate monitoring/traffic mechanism is still required.")

    assert abs(gamma_escape - bw.gamma_escape(h2, idx, s_eta, t_m=T_M, eta=bw.ETA_PIN, g_portal=bw.G_PORTAL)[0]) < 1e-15
    assert dark_dim > 0
    assert max_bright_pole < 0.01 * target_gap
    assert gamma_escape < 0.01 * target_gap
    print("ALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
