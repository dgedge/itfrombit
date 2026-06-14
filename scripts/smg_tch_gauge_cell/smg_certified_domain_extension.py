#!/usr/bin/env python3
r"""SMG certified-domain extension beyond the max-degree KP bound.

The existing escalation theorem uses the deliberately conservative uniform
Kotecky-Preiss sufficient condition on the 4.8.8 plaquette-adjacency graph:

    z(beta) <= 1 / (e (Delta + 1)),  Delta = 8.

This script asks the cheap sharper questions: how far do progressively sharper
single-plaquette polymer criteria go on the actual two-type 4.8.8 face graph,
and do any of them reach beta=1?

Face adjacency facts on the 4.8.8 tiling:
  * a square shares edges with 4 octagons;
  * an octagon shares edges with 4 squares and 4 octagons.

Including self-incompatibility, the type-weighted KP condition is

    z (exp(a_s) + 4 exp(a_o))       <= a_s,
    z (4 exp(a_s) + 5 exp(a_o))     <= a_o.

Maximising over positive type weights gives the improved KP activity.  A second
pass applies the standard independent-neighbourhood/Fernandez-Procacci-style
criterion, using the exact local incompatibility subgraphs: C4 around a square
and alternating C8 around an octagon.

The magnetic activity remains the one already certified in
smg_escalation_theorem.py:

    z(beta) = beta^2 w / (8 C_3),  C_3 = 4/3, w = phi.

Exit 0 = the sharper bounds are computed and the residual beta-domain is tiered.
"""

from __future__ import annotations

import math

from scipy.optimize import differential_evolution, minimize


C3 = 4.0 / 3.0
W_PLAQUETTE = (1.0 + math.sqrt(5.0)) / 2.0
Z_UNIFORM = 1.0 / (9.0 * math.e)


def beta_from_z(z_star: float) -> float:
    return math.sqrt(8.0 * C3 * z_star / W_PLAQUETTE)


def z_limits(a_s: float, a_o: float) -> tuple[float, float]:
    z_s = a_s / (math.exp(a_s) + 4.0 * math.exp(a_o))
    z_o = a_o / (4.0 * math.exp(a_s) + 5.0 * math.exp(a_o))
    return z_s, z_o


def kp_objective(x: tuple[float, float]) -> float:
    a_s, a_o = x
    if a_s <= 0.0 or a_o <= 0.0:
        return 1.0
    return -min(z_limits(a_s, a_o))


def independent_cycle_polynomial(weights: list[float]) -> float:
    """Independent-set polynomial of a cycle with vertex weights."""
    n = len(weights)
    total = 0.0
    for mask in range(1 << n):
        prod = 1.0
        ok = True
        for i, weight in enumerate(weights):
            if (mask >> i) & 1:
                if (mask >> ((i + 1) % n)) & 1:
                    ok = False
                    break
                prod *= weight
        if ok:
            total += prod
    return total


def fp_neighbour_polynomials(mu_s: float, mu_o: float) -> tuple[float, float]:
    """Closed-neighbourhood independent polynomials for the 4.8.8 face graph.

    Square: central S plus four O neighbours whose induced graph is C4.
    Octagon: central O plus four S and four O neighbours whose induced graph is
    alternating C8.
    """
    square = mu_s + independent_cycle_polynomial([mu_o] * 4)
    octagon = mu_o + independent_cycle_polynomial([mu_o, mu_s] * 4)
    return square, octagon


def fp_objective(x: tuple[float, float]) -> float:
    mu_s, mu_o = x
    if mu_s <= 0.0 or mu_o <= 0.0:
        return 1.0
    psi_s, psi_o = fp_neighbour_polynomials(mu_s, mu_o)
    return -min(mu_s / psi_s, mu_o / psi_o)


def dobrushin_objective(x: tuple[float, float]) -> float:
    a_s, a_o = x
    if a_s <= 0.0 or a_o <= 0.0:
        return 1.0
    z_s = (math.exp(a_s) - 1.0) * math.exp(-4.0 * a_o)
    z_o = (math.exp(a_o) - 1.0) * math.exp(-4.0 * a_s - 4.0 * a_o)
    return -min(z_s, z_o)


def optimise(objective, bounds=(1e-8, 8.0), seed=20260612):
    de = differential_evolution(
        objective,
        bounds=[bounds, bounds],
        seed=seed,
        polish=False,
        tol=1e-13,
    )
    return minimize(
        objective,
        de.x,
        method="Nelder-Mead",
        options={"xatol": 1e-14, "fatol": 1e-16, "maxiter": 20000},
    )


def local_face_graph_signatures() -> tuple[tuple[int, list[int]], tuple[int, list[int]]]:
    """Build a small 4.8.8 face graph patch and read the local neighbour graphs."""
    adj: dict[tuple[str, int, int], set[tuple[str, int, int]]] = {}

    def add_node(node):
        adj.setdefault(node, set())

    def add_edge(a, b):
        add_node(a)
        add_node(b)
        adj[a].add(b)
        adj[b].add(a)

    r = 3
    for i in range(-r, r + 1):
        for j in range(-r, r + 1):
            o = ("O", i, j)
            add_node(o)
            add_edge(o, ("O", i + 1, j))
            add_edge(o, ("O", i, j + 1))
            s = ("S", i, j)
            for oi, oj in ((i - 1, j - 1), (i, j - 1), (i - 1, j), (i, j)):
                add_edge(s, ("O", oi, oj))

    def induced_signature(center):
        ns = sorted(adj[center])
        edge_count = 0
        deg = {n: 0 for n in ns}
        for i, a in enumerate(ns):
            for b in ns[i + 1:]:
                if b in adj[a]:
                    edge_count += 1
                    deg[a] += 1
                    deg[b] += 1
        return edge_count, sorted(deg.values())

    return induced_signature(("S", 0, 0)), induced_signature(("O", 0, 0))


def main() -> int:
    sig_s, sig_o = local_face_graph_signatures()
    assert sig_s == (4, [2, 2, 2, 2])
    assert sig_o == (8, [2, 2, 2, 2, 2, 2, 2, 2])

    res = optimise(kp_objective, bounds=(1e-6, 4.0))
    a_s, a_o = map(float, res.x)
    z_s, z_o = z_limits(a_s, a_o)
    z_kp = min(z_s, z_o)
    beta_uniform = beta_from_z(Z_UNIFORM)
    beta_kp = beta_from_z(z_kp)

    dob = optimise(dobrushin_objective, bounds=(1e-8, 4.0), seed=20260613)
    a_ds, a_do = map(float, dob.x)
    z_dob_s = (math.exp(a_ds) - 1.0) * math.exp(-4.0 * a_do)
    z_dob_o = (math.exp(a_do) - 1.0) * math.exp(-4.0 * a_ds - 4.0 * a_do)
    z_dob = min(z_dob_s, z_dob_o)
    beta_dob = beta_from_z(z_dob)

    fp = optimise(fp_objective, bounds=(1e-8, 20.0), seed=20260614)
    mu_s, mu_o = map(float, fp.x)
    psi_s, psi_o = fp_neighbour_polynomials(mu_s, mu_o)
    z_fp = min(mu_s / psi_s, mu_o / psi_o)
    beta_fp = beta_from_z(z_fp)
    z_beta_1 = 1.0 * W_PLAQUETTE / (8.0 * C3)

    print("[0] SMG certified-domain extension: two-type KP on 4.8.8")
    print(f"    uniform max-degree z* = 1/(9e) = {Z_UNIFORM:.12f}")
    print(f"    uniform beta_0        = {beta_uniform:.12f}")
    print(f"    local neighbour graph around S: edges={sig_s[0]}, degrees={sig_s[1]}  (C4)")
    print(f"    local neighbour graph around O: edges={sig_o[0]}, degrees={sig_o[1]}  (C8)")
    print()
    print("[1] two-type 4.8.8 KP optimisation")
    print("    square incompatibilities : self + 4 octagons")
    print("    octagon incompatibilities: self + 4 squares + 4 octagons")
    print(f"    optimal weights: a_square={a_s:.12f}, a_octagon={a_o:.12f}")
    print(f"    active constraints: z_square={z_s:.12f}, z_octagon={z_o:.12f}")
    print(f"    two-type z*       = {z_kp:.12f}")
    print(f"    improvement       = {z_kp / Z_UNIFORM:.6f} in activity")
    print(f"    certified beta_*  = {beta_kp:.12f}")
    print(f"    beta improvement  = {beta_kp / beta_uniform:.6f}")
    print()
    print("[2] sharper single-plaquette criteria")
    print(f"    Dobrushin two-type: z*={z_dob:.12f}, beta_*={beta_dob:.12f}")
    print(f"      weights a_square={a_ds:.12f}, a_octagon={a_do:.12f}")
    print("    FP/independent-neighbourhood two-type criterion:")
    print(f"      mu_square={mu_s:.12f}, mu_octagon={mu_o:.12f}")
    print(f"      Psi_square={psi_s:.12f}, Psi_octagon={psi_o:.12f}")
    print(f"      z*={z_fp:.12f}, beta_*={beta_fp:.12f}")
    print()
    print("[3] domain ledger")
    for beta in (0.5, 0.55, beta_kp, 0.65, beta_fp, 1.0):
        z = beta * beta * W_PLAQUETTE / (8.0 * C3)
        verdict = "FP-certified" if z <= z_fp * (1.0 + 1e-12) else "outside certificate"
        print(f"    beta={beta:.6f}: z={z:.12f}  {verdict}")
    print(f"    beta=1 needs z={z_beta_1:.12f}, i.e. {z_beta_1 / z_fp:.3f}x the best certified activity.")

    assert res.success, res.message
    assert abs(z_s - z_o) < 1e-10
    assert 0.519 < beta_uniform < 0.520
    assert 0.563 < beta_kp < 0.564
    assert beta_kp > beta_uniform
    assert dob.success, dob.message
    assert abs(z_dob_s - z_dob_o) < 1e-10
    assert 0.630 < beta_dob < 0.631
    assert fp.success, fp.message
    assert abs(mu_s / psi_s - mu_o / psi_o) < 1e-10
    assert 0.661 < beta_fp < 0.662
    assert beta_fp > beta_dob > beta_kp
    assert (0.65 * 0.65 * W_PLAQUETTE / (8.0 * C3)) < z_fp
    assert (1.0 * 1.0 * W_PLAQUETTE / (8.0 * C3)) > z_fp

    print()
    print("[VERDICT]")
    print("  The certified strong-coupling domain extends from beta<0.519 to beta<0.661")
    print("  by using the exact two-type local incompatibility structure and the")
    print("  independent-neighbourhood polymer criterion. This certifies beta=0.65,")
    print("  but beta=1 still needs a 2.29x activity improvement. Therefore the")
    print("  weak-coupling frontier is not closed by cheap polymer sharpening; it needs")
    print("  either a genuinely stronger contour/RG reorganisation or a continuum")
    print("  decoupling argument. No finite-row mirror-gap theorem is promoted beyond")
    print("  the beta<0.661 certified domain.")
    print("exit 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
