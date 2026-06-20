#!/usr/bin/env python3
r"""Schwarzschild greybody transfer for the finite horizon source ladder.

The black-hole channel already supplies a local KMS source ladder

    P(F) proportional to g_Q(F) exp[-beta F],

and the freeze-surface audit maps a beta-one source shell to asymptotic
frequencies

    omega r_s = (E_inf/T_H) / (4 pi) = F / (4 pi),

where r_s=2M and beta=1.  This script supplies the next exterior map: the
standard Schwarzschild spin/partial-wave barrier.  It solves

    d^2 psi/dr_*^2 + [omega^2 - V_{s l}(r)] psi = 0,

with

    V_{s l}=f(r)[l(l+1)/r^2 + (1-s^2) r_s/r^3],
    f=1-r_s/r,

for massless integer spin s=0,1,2.  Units set r_s=1.  The near-horizon
solution is normalized to a unit ingoing wave; at large r_* the incoming
amplitude A_in is read off and the greybody factor is Gamma=1/|A_in|^2.

This is not a new QEC principle.  It is the standard exterior Schwarzschild
transfer map applied to the already-derived finite source ladder.  Exit 0 means
the spin/partial-wave greybody barrier is now computed to finite numerical
grade; the remaining flux gates are the horizon attempt-rate normalization and
the premise transferring the 10/27 service bandwidth.
"""

from __future__ import annotations

from dataclasses import dataclass
import cmath
import math


TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}
POSITIVE_F = [f for f in TARGET_GQ if f > 0]
X0 = -26.0
X1 = 180.0
DX = 0.025
LMAX = 9


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def source_ladder(beta: float = 1.0, energy_weighted: bool = False) -> dict[int, float]:
    raw: dict[int, float] = {}
    for f, g in TARGET_GQ.items():
        weight = g * math.exp(-beta * f)
        if energy_weighted:
            weight *= f
        raw[f] = weight
    z = sum(raw.values())
    return dict(sorted((f, v / z) for f, v in raw.items()))


def tortoise(r: float) -> float:
    return r + math.log(r - 1.0)


def f_schw(r: float) -> float:
    return 1.0 - 1.0 / r


def potential(r: float, spin: int, ell: int) -> float:
    f = f_schw(r)
    return f * (ell * (ell + 1.0) / (r * r) + (1.0 - spin * spin) / (r**3))


def omega_for_f(strain_f: int, beta_shell: float = 1.0) -> float:
    """Return omega*r_s for the beta-shell source.

    Since T_H=1/(4*pi*r_s) and E/T_H=beta_shell*F, omega*r_s=beta_shell*F/(4*pi).
    """

    return beta_shell * strain_f / (4.0 * math.pi)


@dataclass(frozen=True)
class ModeResult:
    spin: int
    ell: int
    omega_rs: float
    gamma: float
    flux_error: float
    scaled_flux_error: float


def rk4_step(u: complex, v: complex, r: float, h: float, spin: int, ell: int, omega: float) -> tuple[complex, complex, float]:
    def rhs(uu: complex, vv: complex, rr: float) -> tuple[complex, complex, float]:
        return vv, (potential(rr, spin, ell) - omega * omega) * uu, f_schw(rr)

    k1u, k1v, k1r = rhs(u, v, r)
    k2u, k2v, k2r = rhs(u + 0.5 * h * k1u, v + 0.5 * h * k1v, r + 0.5 * h * k1r)
    k3u, k3v, k3r = rhs(u + 0.5 * h * k2u, v + 0.5 * h * k2v, r + 0.5 * h * k2r)
    k4u, k4v, k4r = rhs(u + h * k3u, v + h * k3v, r + h * k3r)
    return (
        u + (h / 6.0) * (k1u + 2.0 * k2u + 2.0 * k3u + k4u),
        v + (h / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v),
        r + (h / 6.0) * (k1r + 2.0 * k2r + 2.0 * k3r + k4r),
    )


def greybody_factor(spin: int, ell: int, omega: float) -> ModeResult:
    """Compute Gamma_{s ell}(omega) by outward integration from the horizon."""

    r = 1.0 + math.exp(X0 - 1.0)
    # Unit ingoing horizon wave exp(-i omega r_*).
    u = cmath.exp(-1j * omega * X0)
    v = -1j * omega * u

    n = int(round((X1 - X0) / DX))
    h = (X1 - X0) / n
    x = X0
    for _ in range(n):
        u, v, r = rk4_step(u, v, r, h, spin, ell, omega)
        x += h

    # At infinity: u=A_out exp(i omega x)+A_in exp(-i omega x).
    a_in = 0.5 * (u - v / (1j * omega)) * cmath.exp(1j * omega * x)
    a_out = 0.5 * (u + v / (1j * omega)) * cmath.exp(-1j * omega * x)
    gamma_raw = 1.0 / (abs(a_in) ** 2)
    gamma = max(0.0, min(1.0, gamma_raw))
    flux_error = abs((abs(a_in) ** 2 - abs(a_out) ** 2) - 1.0)
    # For almost-total reflection the two asymptotic amplitudes are enormous
    # and the raw difference loses absolute precision.  The quantity relevant
    # to the transmitted flux is the Wronskian error scaled by Gamma.
    scaled_flux_error = flux_error * gamma
    return ModeResult(
        spin=spin,
        ell=ell,
        omega_rs=omega,
        gamma=gamma,
        flux_error=flux_error,
        scaled_flux_error=scaled_flux_error,
    )


def spin_lmin(spin: int) -> int:
    return {0: 0, 1: 1, 2: 2}[spin]


def barrier_spectrum(spin: int, strain_f: int) -> tuple[float, dict[int, ModeResult]]:
    omega = omega_for_f(strain_f)
    modes: dict[int, ModeResult] = {}
    weighted_sum = 0.0
    for ell in range(spin_lmin(spin), LMAX + 1):
        result = greybody_factor(spin, ell, omega)
        modes[ell] = result
        weighted_sum += (2 * ell + 1) * result.gamma
    return weighted_sum, modes


def transferred_ladder(spin: int, source: dict[int, float]) -> tuple[dict[int, float], dict[int, float], float, float]:
    """Apply the standard partial-wave greybody spectral weight to each F line."""

    greybody: dict[int, float] = {}
    raw: dict[int, float] = {}
    max_scaled_flux_error = 0.0
    for strain_f, p in source.items():
        if strain_f == 0:
            greybody[strain_f] = 0.0
            raw[strain_f] = 0.0
            continue
        gsum, modes = barrier_spectrum(spin, strain_f)
        greybody[strain_f] = gsum
        raw[strain_f] = p * gsum
        max_scaled_flux_error = max(max_scaled_flux_error, *(m.scaled_flux_error for m in modes.values()))
    throughput = sum(raw.values())
    norm = throughput if throughput else 1.0
    return dict(sorted((f, v / norm) for f, v in raw.items())), dict(sorted(greybody.items())), throughput, max_scaled_flux_error


def print_mode_table() -> None:
    print("[1] Standard Regge-Wheeler greybody factors at beta-one source lines")
    for spin in (0, 1, 2):
        print(f"    spin {spin}:")
        for strain_f in (3, 4, 6, 9, 12):
            omega = omega_for_f(strain_f)
            gsum, modes = barrier_spectrum(spin, strain_f)
            lead = modes[spin_lmin(spin)]
            tail = modes[LMAX]
            print(
                f"      F={strain_f:2d}, omega*r_s={omega:.6f}, "
                f"sum_l(2l+1)Gamma={gsum:.6f}, "
                f"Gamma_lmin={lead.gamma:.6f}, Gamma_lmax={tail.gamma:.3e}, "
                f"scaled_flux_err<={max(m.scaled_flux_error for m in modes.values()):.2e}"
            )


def main() -> None:
    print("BLACK-HOLE SPIN/PARTIAL-WAVE GREYBODY TRANSFER")
    print("=" * 92)
    print(f"    integration: r_* in [{X0}, {X1}], dx={DX}, lmax={LMAX}")

    print_mode_table()

    print("\n[2] Transfer the finite horizon ladder through the exterior barrier")
    number_source = source_ladder(beta=1.0, energy_weighted=False)
    energy_source = source_ladder(beta=1.0, energy_weighted=True)
    print(f"    number source P(F): {number_source}")
    print(f"    energy source P_E(F): {energy_source}")

    throughputs: dict[int, float] = {}
    for spin in (0, 1, 2):
        transferred, greybody, throughput, err = transferred_ladder(spin, energy_source)
        throughputs[spin] = throughput
        print(f"\n    spin {spin} greybody sums G_s(F): {greybody}")
        print(f"    spin {spin} transferred energy-line distribution: {transferred}")
        print(f"    spin {spin} total spectral throughput (relative units): {throughput:.6f}")
        check(err < 2.0e-3, f"spin {spin} scattering integration conserves flux")
        check(transferred[0] == 0.0, f"spin {spin} zero-energy line carries no transferred flux")

    print("\n[3] Structural checks")
    check(throughputs[0] > throughputs[1] > throughputs[2], "higher spin is more strongly filtered at these low lines")
    for spin in (0, 1, 2):
        g3, _ = barrier_spectrum(spin, 3)
        g12, _ = barrier_spectrum(spin, 12)
        check(g12 > g3, f"spin {spin} barrier opens with increasing line energy")

    print(
        """
[4] VERDICT
    The exterior Schwarzschild greybody transfer is now a computed map rather
    than a placeholder eta_F.

    What closes:
      * The finite source lines F map to standard dimensionless frequencies
        omega r_s = F/(4 pi) at the beta-one freeze shell.
      * For each line, the spin/partial-wave Regge-Wheeler barrier gives a
        numerical spectral weight sum_l (2l+1) Gamma_{s,l}(omega).
      * The observed line distribution is therefore the local KMS ladder
        multiplied by this standard exterior transfer.

    What remains:
      * The absolute Hawking flux coefficient still depends on the horizon
        attempt rate Gamma0, i.e. the conditional 10/27 Landauer-Moore transfer.
      * The script uses the standard integer-spin Schwarzschild barrier.  A
        complete phenomenology paper would add fermions, rotation, charge, and
        the precise species/polarization map of the QEC radiation ledger.

ALL ASSERTIONS PASSED"""
    )


if __name__ == "__main__":
    main()
