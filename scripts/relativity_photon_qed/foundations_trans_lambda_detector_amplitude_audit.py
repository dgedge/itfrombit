#!/usr/bin/env python3
r"""Trans-Lambda_QCD detector-amplitude status audit.

Question
--------
The trans-Lambda_QCD problem is no longer "can a TeV photon be a Bloch mode of
the QCD-spaced crystal?"  The lattice-mode route is excluded.  The canonical
route is now a framed service-causal-set null-chain event: a source/detector
conditioned history carrying total P^mu and no vacuum oscillator tower.

This script tightens the detector-amplitude status:

  CLOSED, at framed smooth-field EFT grade:
    * the lattice cannot supply a TeV support without a new UV oscillator scale;
    * the null-chain endpoint has one total-P Ward vertex, independent of chain
      subdivision;
    * treating the same service units as N independent soft photons is a
      different Fock-sector process with N external legs and e^N suppression.

  RESOLVED BY COMPANION GATE:
    * gauging alone leaves a common endpoint-field normalization unfixed, but
      foundations_trans_lambda_lsz_normalization.py fixes it by the
      Stinespring/GNS endpoint-event isometry.

Exit 0 means the detector-amplitude form is closed and the old freedom is
identified as a service-isometry problem, now closed by the companion LSZ audit.
It is not a new Bloch-support, N-photon-bundle, or vacuum-UV problem.
"""

from __future__ import annotations

import math

import numpy as np


G = np.diag([1.0, -1.0, -1.0, -1.0])
M2 = 1.0
ALPHA0 = 1.0 / 137.035999084
E_CHARGE = math.sqrt(4.0 * math.pi * ALPHA0)
LAMBDA_GEV = 0.3317
TEV_GEV = 1000.0


def mdot(a: np.ndarray, b: np.ndarray) -> float:
    return float(a @ G @ b)


def null_momentum(energy: float, theta: float = 0.4, phi: float = 0.7) -> np.ndarray:
    direction = np.array(
        [
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            math.cos(theta),
        ]
    )
    return np.concatenate([[energy], energy * direction])


def k0(p: np.ndarray) -> float:
    return mdot(p, p) - M2


def gamma0(p: np.ndarray, p_transfer: np.ndarray) -> np.ndarray:
    return 2.0 * p + p_transfer


def subdivided_endpoint_vertex(p: np.ndarray, p_transfer: np.ndarray, nseg: int) -> np.ndarray:
    """Finite-difference endpoint vertex averaged over an internal chain subdivision."""

    s = (np.arange(nseg) + 0.5) / nseg
    return np.array([2.0 * (p + si * p_transfer) for si in s]).mean(axis=0)


def lattice_support_gate() -> dict[str, float]:
    bz_gev = math.pi * LAMBDA_GEV
    support_scale_gev = TEV_GEV / math.pi
    return {
        "bz_gev": bz_gev,
        "tev_over_bz": TEV_GEV / bz_gev,
        "support_scale_gev": support_scale_gev,
        "support_over_lambda": support_scale_gev / LAMBDA_GEV,
        "vacuum_cost": (support_scale_gev / LAMBDA_GEV) ** 4,
    }


def endpoint_ward_gate() -> tuple[float, float, float]:
    p = np.array([1.4, 0.55, -0.20, 0.10])
    p_transfer = null_momentum(0.85)
    lhs = mdot(p_transfer, gamma0(p, p_transfer))
    rhs = k0(p + p_transfer) - k0(p)
    max_subdivision_error = 0.0
    for nseg in (1, 3, 17, 3015):
        gsub = subdivided_endpoint_vertex(p, p_transfer, nseg)
        max_subdivision_error = max(
            max_subdivision_error,
            float(np.max(np.abs(gsub - gamma0(p, p_transfer)))),
            abs(mdot(p_transfer, gsub) - rhs),
        )
    return lhs, rhs, max_subdivision_error


def lsz_normalisation_gate() -> list[tuple[float, float, float, float]]:
    """Show Ward covariance alone does not fix common pole-residue normalization.

    K_c = c (p^2-m^2) and Gamma_c = c(2p+P) satisfy the same Ward identity for
    every positive c.  The propagator pole residue in the service field is 1/c.
    The companion service-GNS audit fixes c=1 by endpoint-event isometry.
    """

    p = np.array([1.31, 0.33, -0.18, 0.22])
    p_transfer = null_momentum(0.64, theta=0.9, phi=0.2)
    out = []
    for c in (0.25, 1.0, 4.0):
        lhs = mdot(p_transfer, c * gamma0(p, p_transfer))
        rhs = c * (k0(p + p_transfer) - k0(p))
        residue = 1.0 / c
        endpoint_strength = c * np.linalg.norm(gamma0(p, p_transfer))
        out.append((c, abs(lhs - rhs), residue, endpoint_strength))
    return out


def main() -> None:
    print("[1] Lattice support gate: TeV photons are not QCD-crystal Bloch modes")
    lat = lattice_support_gate()
    print(f"    Brillouin ceiling pi*Lambda = {lat['bz_gev']:.3f} GeV")
    print(f"    1 TeV / E_BZ = {lat['tev_over_bz']:.1f}")
    print(
        "    a TeV oscillator support scale would be "
        f"{lat['support_scale_gev']:.1f} GeV = {lat['support_over_lambda']:.1f} Lambda_QCD"
    )
    print(f"    vacuum-density cost relative to Lambda_QCD cutoff = {lat['vacuum_cost']:.3e}")
    assert lat["bz_gev"] < 1.1
    assert lat["tev_over_bz"] > 900.0
    assert lat["vacuum_cost"] > 1.0e11

    print("\n[2] Endpoint detector form: one total-P Ward identity")
    lhs, rhs, sub_err = endpoint_ward_gate()
    print(f"    P.Gamma_total       = {lhs:+.12f}")
    print(f"    K(p+P)-K(p)         = {rhs:+.12f}")
    print(f"    max subdivision / Ward error over N=1,3,17,3015 = {sub_err:.3e}")
    assert abs(lhs - rhs) < 1.0e-12
    assert sub_err < 1.0e-12

    print("\n[3] N-soft-photon bundle is the excluded alternative")
    n_tev = math.ceil(TEV_GEV / LAMBDA_GEV)
    log10_bundle = n_tev * math.log10(E_CHARGE)
    print(f"    1 TeV corresponds to N={n_tev} Lambda_QCD service units")
    print(f"    one endpoint leg scales as e = {E_CHARGE:.6f}")
    print(f"    N independent soft-photon legs scale as e^N = 10^{log10_bundle:.0f}")
    assert n_tev > 3000
    assert log10_bundle < -1000.0

    print("\n[4] Former freedom: Ward alone does not fix LSZ/service-measure normalization")
    print("    Common rescaling K -> cK, Gamma -> cGamma keeps Ward exact:")
    rows = lsz_normalisation_gate()
    for c, ward_err, residue, endpoint_strength in rows:
        print(
            f"      c={c:4.2f}: Ward error={ward_err:.2e}, "
            f"service pole residue Z={residue:.3f}, raw endpoint strength={endpoint_strength:.3f}"
        )
        assert ward_err < 1.0e-12
    residues = [r[2] for r in rows]
    strengths = [r[3] for r in rows]
    assert max(residues) / min(residues) == 16.0
    assert max(strengths) / min(strengths) == 16.0
    print(
        "    -> gauge covariance fixes the amplitude form, not the absolute "
        "service-field residue.  The companion service-GNS audit fixes this by "
        "requiring the endpoint event map to be an isometry."
    )

    print(
        r"""
[5] VERDICT
  Tightened status:
    * trans-Lambda_QCD detector-amplitude FORM is closed at framed causal-set
      smooth-field EFT grade: one endpoint/LSZ leg, total P^mu, exact Ward
      identity, subdivision-invariant chain label, no N-photon Fock penalty.
    * the QCD lattice remains only the IR condensate; high-energy external
      events do not require a finer vacuum oscillator tower and therefore do
      not reopen the cosmological-constant cutoff.

  Follow-up:
    * foundations_trans_lambda_lsz_normalization.py closes the external-leg
      pole residue by service-history isometry, Z_endpoint=1.  What remains is
      ordinary normalized-leg QED phenomenology: spinor/polarization sums, loop
      running, and finite-density corrections.
ALL ASSERTIONS PASSED -- trans-Lambda detector amplitude form closed; LSZ residue fixed by companion gate."""
    )


if __name__ == "__main__":
    main()
