#!/usr/bin/env python3
r"""Normalized-leg QED phenomenology for trans-Lambda_QCD null-chain events.

Question
--------
After the service-GNS audit fixes the endpoint LSZ residue,

    Z_endpoint = 1,

what remains of "high-energy QED phenomenology" around the trans-Lambda null
chain?  Is there another microscopic normalization gate hidden in spin sums,
polarization sums, loop running, or finite-density propagation?

Result
------
At the scoped, framed-causal-set EFT level: no.  Once the null-chain endpoint is
one normalized external leg carrying total P^mu, the surrounding machinery is
ordinary QED:

  * physical photon polarization sums reduce conserved detector currents to the
    two transverse polarizations, independent of axial-gauge reference vector;
  * Dirac spinor sums / Ward contractions work with the already-supplied
    E_{1/2} spin frame;
  * loop running is the standard charged-field vacuum-polarization running of
    alpha(Q), not a new oscillator support or vacuum zero-point tower;
  * finite-density corrections are medium self-energies / opacity, e.g. plasma
    frequency and gamma-gamma pair-production thresholds, not failures of the
    vacuum external-leg representation.

Scope
-----
This is a structural interface audit, not a precision phenomenology package.
It does not compute every astrophysical gamma-ray transfer function, higher-loop
SM correction, or detector-specific cross section.  It says those calculations
can be done using standard normalized-leg QED once the endpoint leg is supplied.
"""

from __future__ import annotations

import math

import numpy as np


ETA = np.diag([1.0, -1.0, -1.0, -1.0])
I2 = np.eye(2, dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
I4 = np.eye(4, dtype=complex)
GAMMA = [
    np.block([[I2, Z2], [Z2, -I2]]),
    np.block([[Z2, SX], [-SX, Z2]]),
    np.block([[Z2, SY], [-SY, Z2]]),
    np.block([[Z2, SZ], [-SZ, Z2]]),
]

ALPHA0_INV = 137.035999084
ALPHA0 = 1.0 / ALPHA0_INV
LAMBDA_GEV = 0.3317
TEV_GEV = 1000.0
ME_EV = 510_998.95


def mdot(a: np.ndarray, b: np.ndarray) -> float:
    return float(a @ ETA @ b)


def lower(v: np.ndarray) -> np.ndarray:
    return ETA @ v


def slash(p: np.ndarray) -> np.ndarray:
    return sum(GAMMA[mu] * p[mu] for mu in range(4))


def axial_photon_sum(p: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Physical-polarization completeness in axial gauge.

    Sum eps^mu eps^nu =
      -eta^mu nu + (p^mu n^nu+n^mu p^nu)/(p.n) - n^2 p^mu p^nu/(p.n)^2.
    """

    pn = mdot(p, n)
    n2 = mdot(n, n)
    return -ETA + (np.outer(p, n) + np.outer(n, p)) / pn - n2 * np.outer(p, p) / (pn * pn)


def polarization_gate() -> tuple[float, float, float]:
    """A conserved detector current sees only two transverse photon states."""

    p = np.array([TEV_GEV, 0.0, 0.0, TEV_GEV])
    n = np.array([1.0, 0.0, 0.0, -1.0])
    eps1 = np.array([0.0, 1.0, 0.0, 0.0])
    eps2 = np.array([0.0, 0.0, 1.0, 0.0])
    physical_sum = np.outer(eps1, eps1) + np.outer(eps2, eps2)
    axial_sum = axial_photon_sum(p, n)

    # P.J=0 but J may include a pure-gauge longitudinal component along P.
    j = np.array([0.73, 0.21, -0.46, 0.73])
    assert abs(mdot(p, j)) < 1.0e-9
    jl = lower(j)
    physical = float(jl @ physical_sum @ jl)
    axial = float(jl @ axial_sum @ jl)
    target = -mdot(j, j)
    return physical, axial, target


def on_shell(spatial: list[float], mass: float) -> np.ndarray:
    q = np.array(spatial, dtype=float)
    return np.concatenate([[math.sqrt(mass * mass + float(q @ q))], q])


def spinor_ward_gate() -> tuple[float, float]:
    """Spin-summed Dirac current obeys the Ward contraction.

    For p^2=q^2=m^2 and k=q-p, individual matrix elements satisfy
    k_mu ubar(q) gamma^mu u(p)=0.  The spin-summed version is checked by the
    trace Tr[(qslash+m) kslash (pslash+m) kslash].
    """

    m = 0.7
    p = on_shell([0.30, -0.20, 0.40], m)
    q = on_shell([-0.10, 0.50, 0.20], m)
    k = q - p
    ward_trace = np.trace((slash(q) + m * I4) @ slash(k) @ (slash(p) + m * I4) @ slash(k))
    p_shell = abs(mdot(p, p) - m * m)
    q_shell = abs(mdot(q, q) - m * m)
    return float(abs(ward_trace)), max(p_shell, q_shell)


def running_alpha_one_loop(q_gev: float) -> tuple[float, float, float]:
    """One-loop threshold model for alpha(Q), for orientation not precision."""

    fermions = [
        ("e", 0.00051099895, -1.0, 1),
        ("mu", 0.1056583755, -1.0, 1),
        ("tau", 1.77686, -1.0, 1),
        ("u", 0.00216, 2.0 / 3.0, 3),
        ("d", 0.00467, -1.0 / 3.0, 3),
        ("s", 0.093, -1.0 / 3.0, 3),
        ("c", 1.27, 2.0 / 3.0, 3),
        ("b", 4.18, -1.0 / 3.0, 3),
        ("t", 172.76, 2.0 / 3.0, 3),
    ]
    sigma = 0.0
    active = 0
    for _name, mass, charge, nc in fermions:
        if q_gev > mass:
            sigma += nc * charge * charge * math.log(q_gev / mass)
            active += 1
    alpha_inv_q = ALPHA0_INV - (2.0 / (3.0 * math.pi)) * sigma
    return alpha_inv_q, sigma, active


def finite_density_gate(e_gev: float = TEV_GEV) -> dict[str, float]:
    """Representative medium corrections for a TeV external photon."""

    e_ev = e_gev * 1.0e9
    # Plasma frequency: hbar omega_p = 3.713e-11 eV * sqrt(n_e/cm^-3).
    n_igm = 1.0e-7
    omega_p_ev = 3.713e-11 * math.sqrt(n_igm)
    plasma_dv = -0.5 * (omega_p_ev / e_ev) ** 2

    # Head-on gamma-gamma pair production threshold: epsilon >= m_e^2/E.
    epsilon_thr_ev = ME_EV * ME_EV / e_ev
    cmb_epsilon_ev = 2.701 * 8.617333262e-5 * 2.7255
    e_thr_cmb_ev = ME_EV * ME_EV / cmb_epsilon_ev

    return {
        "omega_p_ev": omega_p_ev,
        "plasma_dv": plasma_dv,
        "epsilon_thr_ev": epsilon_thr_ev,
        "cmb_mean_ev": cmb_epsilon_ev,
        "e_thr_cmb_tev": e_thr_cmb_ev / 1.0e12,
    }


def main() -> None:
    print("[1] Physical-polarization sum for the normalized null-chain leg")
    physical, axial, target = polarization_gate()
    print(f"    physical two-pol sum J.Pi.J = {physical:.12f}")
    print(f"    axial-gauge completeness    = {axial:.12f}")
    print(f"    conserved-current target -J^2 = {target:.12f}")
    assert abs(physical - target) < 1.0e-12
    assert abs(axial - target) < 1.0e-12
    print("    -> detector currents see one photon leg with two transverse polarizations.")

    print("\n[2] Spinor sums and Ward contraction")
    ward_trace, shell_err = spinor_ward_gate()
    print(f"    max on-shell error = {shell_err:.3e}")
    print(f"    spin-summed |k_mu J^mu|^2 trace = {ward_trace:.3e}")
    assert shell_err < 1.0e-14
    assert ward_trace < 1.0e-12
    print("    -> with the E_{1/2} spin frame, standard Dirac spin sums are compatible")
    print("       with the same normalized endpoint photon leg.")

    print("\n[3] Loop running is ordinary charged-field vacuum polarization")
    for q in (91.1876, 1000.0):
        alpha_inv_q, sigma, active = running_alpha_one_loop(q)
        print(
            f"    Q={q:8.3f} GeV: active thresholds={active}, "
            f"sum Nc Q_f^2 ln(Q/m_f)={sigma:.3f}, alpha^-1(Q)~{alpha_inv_q:.3f}"
        )
        assert alpha_inv_q < ALPHA0_INV
    print("    -> this is standard EFT running around a normalized external leg;")
    print("       it does not introduce a finer vacuum oscillator tower.")

    print("\n[4] Finite-density effects are medium corrections, not support failures")
    fd = finite_density_gate()
    print(f"    IGM plasma energy for n_e=1e-7 cm^-3: hbar omega_p={fd['omega_p_ev']:.3e} eV")
    print(f"    TeV photon plasma velocity shift: dv/v={fd['plasma_dv']:.3e}")
    print(f"    gamma-gamma target threshold at 1 TeV: epsilon_thr={fd['epsilon_thr_ev']:.3f} eV")
    print(f"    CMB mean photon energy: epsilon_CMB={fd['cmb_mean_ev']:.3e} eV")
    print(f"    CMB pair-production threshold: E_gamma~{fd['e_thr_cmb_tev']:.1f} TeV")
    assert abs(fd["plasma_dv"]) < 1.0e-45
    assert fd["epsilon_thr_ev"] > 0.1
    assert fd["e_thr_cmb_tev"] > 100.0
    print("    -> TeV propagation is governed by standard EBL opacity / plasma physics;")
    print("       those are environment self-energies on the same normalized leg.")

    n_service = math.ceil(TEV_GEV / LAMBDA_GEV)
    print(
        rf"""
[5] VERDICT
  CLOSED AS AN INTERFACE, not as a full phenomenology monograph:
    With Z_endpoint=1, the trans-Lambda null-chain event supplies exactly one
    normalized external photon leg carrying total P^mu.  The standard QED
    interfaces then apply:
      - two physical transverse polarizations for conserved detector currents;
      - standard Dirac spinor sums and Ward contractions using the E_{{1/2}} frame;
      - ordinary charged-field loop running of alpha(Q);
      - ordinary finite-density/opacity corrections in media.

    For a 1 TeV event the service ledger contains N={n_service} Lambda_QCD
    service/action units, but those are internal history refinements of one
    endpoint leg, not N independent soft-photon legs.

  STILL NOT CLAIMED:
    a precision gamma-ray propagation code, all higher-loop SM corrections, or a
    detector-specific S-matrix library.  Those are standard normalized-leg EFT
    calculations around the now-fixed endpoint representation, not remaining
    microscopic support/LSZ normalisation gates.
ALL ASSERTIONS PASSED -- normalized-leg QED interface checks pass."""
    )


if __name__ == "__main__":
    main()
