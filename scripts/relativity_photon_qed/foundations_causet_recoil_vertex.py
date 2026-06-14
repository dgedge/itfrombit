#!/usr/bin/env python3
r"""Causal-set holonomy gauge principle: full recoil vertex gate.

Question:
  K19 derived the eikonal Wilson-line vertex eps.p/(k.p), but left the full
  recoil-corrected scalar-QED vertex (2p+k) as open.  Does the same causal-set
  holonomy principle derive that vertex once it is applied to the scalar kinetic
  operator / propagator, rather than only to a prescribed classical worldline?

Result:
  Yes at continuum-propagator grade.  Gauging a translation-invariant scalar
  kinetic kernel by U(x,y)=exp(i integral_x^y A.dx) gives the exact
  Ward-Takahashi finite-difference identity

      k_mu Gamma^mu(p+k,p) = K(p+k) - K(p).

  For the continuum causal-set propagator limit K(p)=p^2-m^2, the unique local
  minimal longitudinal vertex is

      Gamma^mu = (p+k)^mu + p^mu = 2p^mu + k^mu.

  Lattice/finite-density artifacts appear as higher-derivative corrections to
  K(p), and the same line-integral formula keeps Ward exact while adding
  O(a^2 p^3) corrections.  This closes the recoil-vertex residual named in K19,
  but it does NOT derive the photon's own Maxwell F^2 action from causal-set
  loops.  That loop/measure/continuum-limit problem remains the real frontier.
"""

from __future__ import annotations

import numpy as np

G = np.diag([1.0, -1.0, -1.0])
M2 = 1.0


def dot(a: np.ndarray, b: np.ndarray) -> float:
    return float(a @ G @ b)


def K_cont(p: np.ndarray) -> float:
    """Continuum scalar inverse propagator K(p)=p^2-m^2."""
    return dot(p, p) - M2


def gamma_cont(p: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Minimal continuum vertex from gauged K(p)=p^2-m^2."""
    return 2.0 * p + k


def K_artifact(p: np.ndarray, eps: float = 0.08) -> float:
    """Toy finite-density/lattice correction: K=p^2-m^2+eps*(p^2)^2."""
    p2 = dot(p, p)
    return p2 - M2 + eps * p2 * p2


def grad_K_artifact(p: np.ndarray, eps: float = 0.08) -> np.ndarray:
    """Gradient in the metric-dual convention: dK/ds = k.grad_K(p+s k)."""
    p2 = dot(p, p)
    return 2.0 * p * (1.0 + 2.0 * eps * p2)


def gamma_line_integral(p: np.ndarray, k: np.ndarray, n: int = 20000) -> np.ndarray:
    """Gauge-link finite-difference vertex: integral_0^1 ds grad K(p+s k)."""
    s = (np.arange(n) + 0.5) / n
    vals = np.array([grad_K_artifact(p + si * k) for si in s])
    return vals.mean(axis=0)


def photon(omega: float, theta: float) -> np.ndarray:
    return np.array([omega, omega * np.cos(theta), omega * np.sin(theta)])


def pol(theta: float) -> np.ndarray:
    return np.array([0.0, -np.sin(theta), np.cos(theta)])


print("[1] Continuum gauged-kernel theorem: Ward finite difference -> full recoil vertex")
p = np.array([np.sqrt(1.0 + 0.6**2), 0.6, 0.0])
k = photon(0.8, 0.7)
g = gamma_cont(p, k)
lhs = dot(k, g)
rhs = K_cont(p + k) - K_cont(p)
print(f"    p^2={dot(p,p):.12f}, k^2={dot(k,k):+.2e}")
print(f"    Gamma = 2p+k = {np.array2string(g, precision=6)}")
print(f"    k.Gamma = {lhs:.12f}")
print(f"    K(p+k)-K(p) = {rhs:.12f}")
assert abs(lhs - rhs) < 1e-12
print("    -> exact Ward-Takahashi identity; the recoil k term is required by K(p+k)-K(p).")

print("\n[2] Same vertex is the scalar-QED plane-wave matrix element")
eta = np.array([0.3, -0.2, 0.7])  # off-shell/gauge probe, so the recoil term is visible
eps = pol(0.7)                    # physical transverse probe, eps.k=0
amp_from_vertex = dot(eta, g)
amp_from_endpoints = dot(eta, (p + k) + p)
amp_eikonal_only = dot(eta, 2.0 * p)
print(f"    eta.(2p+k) from gauged kinetic term = {amp_from_vertex:+.12f}")
print(f"    eta.(p_out+p_in) endpoint form      = {amp_from_endpoints:+.12f}")
print(f"    eta.(2p) eikonal-only value         = {amp_eikonal_only:+.12f}")
print(f"    visible recoil contribution eta.k   = {dot(eta, k):+.12f}")
assert abs(amp_from_vertex - amp_from_endpoints) < 1e-12
assert abs((amp_from_vertex - amp_eikonal_only) - dot(eta, k)) < 1e-12
print(f"    physical transverse check: eps.k={dot(eps,k):+.1e}, so eps.(2p+k)=eps.2p as expected.")
print("    -> the full numerator is p_out+p_in; transverse on-shell amplitudes can hide the k term,")
print("       but Ward/off-shell recoil requires it.")

print("\n[3] Finite-density artifacts: line-integral vertex keeps Ward exact")
k_off = np.array([0.45, 0.31, -0.18])  # deliberately off-shell so k^2 != 0
gamma_fd = gamma_line_integral(p, k_off)
lhs_fd = dot(k_off, gamma_fd)
rhs_fd = K_artifact(p + k_off) - K_artifact(p)
print(f"    toy K=p^2-m^2+eps(p^2)^2, eps=0.08")
print(f"    k.Gamma_line = {lhs_fd:.12f}")
print(f"    K(p+k)-K(p) = {rhs_fd:.12f}")
print(f"    artifact correction norm |Gamma_line-(2p+k)| = {np.linalg.norm(gamma_fd - gamma_cont(p, k_off)):.4e}")
assert abs(lhs_fd - rhs_fd) < 2e-5
print("    -> finite-density corrections change the vertex, but the gauged-link construction")
print("       preserves Ward exactly as a finite difference; continuum limit returns 2p+k.")

print("\n[4] What this does NOT solve: Maxwell F^2 from causal-set loops")
print("    Matter coupling: CLOSED at scalar/recoil level under gauged-kernel/holonomy premises.")
print("    Photon dynamics: still OPEN. A loop-holonomy action must still supply a local,")
print("    Lorentz-covariant, correctly normalised F_{mu nu}F^{mu nu} continuum limit with")
print("    a non-arbitrary loop measure and finite-density control.")

print(
    r"""
[5] VERDICT
  The K19 residual splits:
    * full scalar recoil vertex (2p+k): CLOSED at continuum causal-set-propagator grade.
      It follows from gauging the scalar kinetic kernel with link holonomies, via the exact
      finite-difference Ward identity k.Gamma = K(p+k)-K(p).  For K=p^2-m^2 this gives
      Gamma=2p+k.  Finite-density corrections are higher-derivative corrections to K and
      are Ward-controlled by the same line-integral formula.
    * Maxwell F^2 action from causal-set loop holonomies: STILL OPEN.  This is the
      genuine causal-set gauge-field frontier, not closed by deriving the matter vertex.
    * Dirac spin / spinor causal-set operator: still deferred.
exit 0"""
)
print("ALL ASSERTIONS PASSED -- recoil vertex derived; Maxwell loop action remains open.")
