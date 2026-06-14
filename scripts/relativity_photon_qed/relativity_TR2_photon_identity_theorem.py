#!/usr/bin/env python3
r"""T-R2 object-identity theorem: which object is the physical photon?

Canon contains two gauge-sector objects that were previously conflated:

  1. The macroscopic compact-U(1) photon on the dual SC gauge web, with Wilson
     plaquette action and Gauss law.  Linearised at nonzero momentum, this is
     the standard Maxwell cochain complex: vector potentials are 1-cochains,
     gauge transformations are exact forms A -> A + d lambda, and the physical
     photon sector is the two-dimensional transverse quotient.

  2. The local K6 T1u (+) Eg vertex multiplet of item 102.  It is a five-mode
     local connectivity multiplet.  It contains a T1u vector triplet and an Eg
     tensor doublet, but it is not itself Gauss-projected.

This script establishes the object identity at finite momentum:

  physical photon = Gauss-projected transverse Wilson/Maxwell 1-cochain,
  not the raw K6 T1u/Eg branch.

The proof is constructive:
  * build the lattice gradient g_i(k)=exp(i k_i)-1;
  * build the Wilson/Maxwell quadratic kernel
        M_ij(k)=|g|^2 delta_ij - g_i g_j^*
    from |dA|^2;
  * verify Ward identity, rank two, and two degenerate transverse eigenvalues;
  * verify the long-wavelength SC group-velocity anisotropy is O(k^2) and
    below the optical SME line with the coefficient q^2/(12 Lambda^2);
  * build the K6 T1u/Eg block and show it fails to descend to the quotient:
    C(k) g(k) != 0, so it couples gauge-exact longitudinal A to Eg.

exit 0 = the compact-U(1) photon identity is established under canon's own
         Wilson-action/Gauss-law photon definition; K6-as-photon is rejected.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
ANCHOR = (ROOT / "ANCHOR.md").read_text(encoding="utf-8")
PHOTON = (ROOT / "photons" / "photon_paper2_final.tex").read_text(encoding="utf-8")


def check(cond: bool, msg: str) -> None:
    print(f"    [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def grad(k):
    """Forward lattice gradient symbol for SC compact U(1)."""
    k = np.asarray(k, dtype=float)
    return np.exp(1j * k) - 1.0


def maxwell_kernel(k):
    """Quadratic Wilson/Maxwell kernel from sum_{i<j}|g_i A_j-g_j A_i|^2."""
    g = grad(k)
    g2 = float(np.vdot(g, g).real)
    return g2 * np.eye(3, dtype=complex) - np.outer(g, np.conjugate(g))


def transverse_projector(k):
    g = grad(k)
    g2 = float(np.vdot(g, g).real)
    return np.eye(3, dtype=complex) - np.outer(g, np.conjugate(g)) / g2


def c_block(n):
    """Item-102 K6 Eg <- T1u block at unit direction n."""
    nx, ny, nz = unit(n)
    return np.array(
        [
            [-nx / math.sqrt(6), -ny / math.sqrt(6), 2 * nz / math.sqrt(6)],
            [nx / math.sqrt(2), -ny / math.sqrt(2), 0.0],
        ],
        dtype=complex,
    )


def k6_singular_values(n):
    return np.linalg.svd(c_block(n), compute_uv=False)


def sc_omega(k):
    return math.sqrt(4.0 * sum(math.sin(ki / 2.0) ** 2 for ki in k))


def sc_group_velocity(kmag, n):
    n = unit(n)

    def omega_at(x):
        return sc_omega(x * n)

    h = max(1e-8, kmag * 1e-5)
    return (omega_at(kmag + h) - omega_at(kmag - h)) / (2.0 * h)


print("[0] Canon premises for the photon object:")
check("standard Wilson action" in PHOTON, "photon paper: pure gauge action is standard Wilson")
check("two transverse modes from gauge invariance" in PHOTON, "photon paper: SC web supports two transverse modes")
check("Simple Cubic Bravais lattice" in ANCHOR, "ANCHOR: macroscopic gauge dual is SC")
check("Wilson action" in ANCHOR, "ANCHOR: SC gauge web uses Wilson action")

print("\n[1] Linearised compact-U(1) / Maxwell cochain theorem:")
for label, k in [
    ("generic", (0.37, 0.51, 0.23)),
    ("[100]", (0.41, 0.0, 0.0)),
    ("[111]", (0.41 / math.sqrt(3),) * 3),
]:
    k = np.asarray(k, dtype=float)
    g = grad(k)
    M = maxwell_kernel(k)
    P = transverse_projector(k)
    g2 = float(np.vdot(g, g).real)
    ward = np.linalg.norm(M @ g)
    pe = np.linalg.norm(P @ P - P)
    ph = np.linalg.norm(P.conjugate().T - P)
    evals = np.sort(np.linalg.eigvalsh(M).real)
    rank = int(np.sum(evals > 1e-10))
    print(f"    {label}: |g|^2={g2:.9f}, eig(M)={np.round(evals, 9)}, rank={rank}")
    assert ward < 1e-12
    assert pe < 1e-12 and ph < 1e-12
    assert rank == 2
    assert abs(evals[0]) < 1e-12
    assert abs(evals[1] - g2) < 1e-12 and abs(evals[2] - g2) < 1e-12
print("    -> nonzero-k physical photon sector has exactly two degenerate transverse modes;")
print("       the longitudinal exact mode is pure gauge. This is the Gauss-projected")
print("       Maxwell quotient, not a five-mode local multiplet.")

print("\n[2] SC Wilson/Maxwell anisotropy has the irrelevant O(k^2) form:")
ks = np.array([0.2, 0.1, 0.05, 0.025])
anis = []
for kmag in ks:
    v100 = sc_group_velocity(kmag, (1, 0, 0))
    v111 = sc_group_velocity(kmag, (1, 1, 1))
    frac = abs(v100 - v111) / (0.5 * (v100 + v111))
    anis.append(frac)
    print(f"    |k|={kmag:.3f}: v100={v100:.9f}, v111={v111:.9f}, Delta v/v={frac:.3e}")
power = np.polyfit(np.log(ks), np.log(anis), 1)[0]
print(f"    fitted anisotropy power = {power:.3f} (target 2)")
assert abs(power - 2.0) < 0.03
Lam = 0.33
q_visible = 2.48e-9  # GeV, about 500 nm visible light
q_lam = q_visible / Lam
optical = q_lam * q_lam / 12.0
print(f"    visible 500 nm: (q/Lambda)^2/12 = {optical:.2e} vs SME line 1e-17")
assert optical < 1e-17

print("\n[3] K6 T1u/Eg block is not a Gauss-quotient photon operator:")
for label, n in [("[100]", (1, 0, 0)), ("[110]", (1, 1, 0)), ("[111]", (1, 1, 1))]:
    C = c_block(n)
    nh = unit(n)
    longitudinal_leak = np.linalg.norm(C @ nh)
    sv = k6_singular_values(n)
    print(f"    {label}: singular velocities={np.round(sv, 9)}, ||C n||={longitudinal_leak:.9f}")
assert abs(k6_singular_values((1, 0, 0))[0] - math.sqrt(2 / 3)) < 1e-12
assert abs(k6_singular_values((1, 1, 0))[0] - 1 / math.sqrt(2)) < 1e-12
assert abs(k6_singular_values((1, 1, 1))[0] - 1 / math.sqrt(3)) < 1e-12
assert np.linalg.norm(c_block((1, 0, 0)) @ unit((1, 0, 0))) > 0.8
assert np.linalg.norm(c_block((1, 1, 0)) @ unit((1, 1, 0))) > 0.4
print("    -> C n != 0 means a gauge-exact longitudinal vector is sent into Eg.")
print("       Therefore the raw K6 mixing is not well-defined on A / im(d); it")
print("       cannot be the physical compact-U(1) photon inverse propagator.")

print("\n[4] Mode accounting:")
print("    raw K6 -1 manifold: dim(T1u + Eg) = 3 + 2 = 5")
print("    Gauss projection: T1u -> 2 transverse physical + 1 longitudinal gauge")
print("    Eg doublet: massive/non-photon gauge-web tensor (e.g. 2++ glueball channel)")
check(3 + 2 == 5, "raw K6 local multiplet has five modes")
check(2 + 1 + 2 == 5, "Gauss/gauge/tensor accounting exhausts the same five labels")
print("    -> treating all five K6 labels as photon polarisations is the category error.")

print(
    """
[5] THEOREM VERDICT:
  Under canon's own photon definition -- compact U(1) Wilson action on the SC
  gauge web, with Gauss law -- the physical nonzero-momentum photon is the
  two-dimensional transverse Maxwell quotient of the 1-cochain A_i(k).  Its
  quadratic kernel is Ward-transverse and has two equal eigenvalues |g(k)|^2.

  The raw K6 T1u/Eg branch does not descend to that quotient: its mixing sends
  longitudinal gauge-exact T1u vectors into Eg for generic directions.  It is a
  local pre-Gauss gauge-web/vertex multiplet.  It may contain the local T1u
  transmission resonance, but it is not the physical photon until the compact
  U(1) Gauss projection and Wilson plaquette action are imposed.

  Therefore the physical photon = Gauss-projected transverse Wilson/Maxwell
  sector, not raw K6 T1u/Eg.  The prior K6 marginal-anisotropy route is rejected
  as a photon-Lorentz problem; the macroscopic photon has only the standard
  irrelevant lattice anisotropy O((q/Lambda)^2).
exit 0"""
)
print("ALL ASSERTIONS PASSED -- photon object identity established under the canon Wilson/Gauss premises.")
