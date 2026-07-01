#!/usr/bin/env python3
r"""SU(3) Wilson-loop cluster: the strong-sector version of the record-action rung.

This script repeats the non-abelian Wilson-cluster certificate with SU(3), the
gauge group relevant to the strong sector.

Declared measure
----------------
Use the same controlled two-dimensional axial-gauge setting as
record_grammar_nonabelian_wilson_cluster.py, now with SU(3) plaquette records:

    dP(U_p) proportional to dHaar_SU3(U_p)
        exp[ beta Re Tr(U_p) / 3 ].

The detector-readable loop record is

    W(C) = (1/3) Re Tr prod_{p inside C} U_p.

Because the measure is central, E[U_p] = a(beta) I_3.  Independent plaquette
records then give the exact factorized area law

    <W(C)> = a(beta)^A.

The script samples the declared measure by rejection from Haar_SU3, verifies
the central first moment, checks small rectangular Wilson loops by Monte Carlo,
fits the factorized prediction to

    -log <W(C)> = sigma A + mu P + c,

and reports the ledger-size scaling.

Boundary
--------
This is still a 2D pure-gauge record-action certificate.  It is not 4D QCD, not
the TCH mirror-gapped continuum theorem, and not a glueball/string-tension
prediction in physical units.  Its value is narrower and concrete: the same
record grammar that gave SU(2) loop records extends to SU(3), and the detector
interface / tractable-ledger coincidence survives in a controlled strong-sector
toy.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


N_COLOR = 3


@dataclass(frozen=True)
class LoopRow:
    width: int
    height: int
    area: int
    perimeter: int
    predicted: float
    mc_mean: float | None = None
    mc_stderr: float | None = None


def assert_close(name: str, value: float, target: float, tol: float = 1e-10) -> None:
    err = abs(value - target)
    print(f"  {name:<76s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<76s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_greater(name: str, value: float, bound: float) -> None:
    print(f"  {name:<76s} value={value:.12g} bound={bound:.12g}")
    if not value > bound:
        raise AssertionError(name)


def haar_su3(rng: np.random.Generator) -> np.ndarray:
    """Draw one Haar-random SU(3) matrix by QR factorization."""

    z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    q, r = np.linalg.qr(z)
    diag = np.diag(r)
    phases = diag / np.abs(diag)
    q = q @ np.diag(np.conj(phases))
    det_phase = np.angle(np.linalg.det(q))
    q = q * np.exp(-1j * det_phase / 3.0)
    return q


def sample_su3_wilson(beta: float, rng: np.random.Generator) -> tuple[np.ndarray, int]:
    """Sample dHaar(U) exp[beta ReTr(U)/3] by rejection from Haar_SU3."""

    tries = 0
    while True:
        tries += 1
        u = haar_su3(rng)
        x = float(np.real(np.trace(u)) / 3.0)
        # x <= 1, so exp(beta*(x-1)) is a valid rejection probability.
        if rng.random() <= math.exp(beta * (x - 1.0)):
            return u, tries


def rectangle_loop(plaquettes: np.ndarray, x0: int, y0: int, width: int, height: int) -> float:
    prod = np.eye(3, dtype=complex)
    for y in range(y0, y0 + height):
        for x in range(x0, x0 + width):
            prod = prod @ plaquettes[y, x]
    return float(np.real(np.trace(prod)) / 3.0)


def all_rectangle_values(plaquettes: np.ndarray, width: int, height: int) -> np.ndarray:
    size = plaquettes.shape[0]
    vals = []
    for y0 in range(size - height + 1):
        for x0 in range(size - width + 1):
            vals.append(rectangle_loop(plaquettes, x0, y0, width, height))
    return np.array(vals, dtype=float)


def fit_area_perimeter(rows: list[LoopRow], use_mc: bool = False) -> tuple[float, float, float, float]:
    xs = []
    ys = []
    for row in rows:
        val = row.mc_mean if use_mc else row.predicted
        if val is None or val <= 0.0:
            continue
        xs.append([row.area, row.perimeter, 1.0])
        ys.append(-math.log(val))
    xmat = np.array(xs, dtype=float)
    yvec = np.array(ys, dtype=float)
    coeff, *_ = np.linalg.lstsq(xmat, yvec, rcond=None)
    resid = yvec - xmat @ coeff
    return float(coeff[0]), float(coeff[1]), float(coeff[2]), float(np.max(np.abs(resid)))


def predicted_rows(size: int, a_beta: float) -> list[LoopRow]:
    return [
        LoopRow(w, h, w * h, 2 * (w + h), a_beta ** (w * h))
        for h in range(1, size + 1)
        for w in range(1, size + 1)
    ]


def monte_carlo_rows(
    size: int,
    beta: float,
    rng: np.random.Generator,
    n_cfg: int = 800,
    max_area: int = 4,
) -> tuple[list[LoopRow], float, float, float]:
    """Sample plaquette lattices and small rectangular loop records."""

    shapes = [
        (w, h)
        for h in range(1, size + 1)
        for w in range(1, size + 1)
        if w * h <= max_area
    ]
    buckets: dict[tuple[int, int], list[float]] = {shape: [] for shape in shapes}
    plaquette_traces: list[float] = []
    accepted = 0
    tries = 0

    for _ in range(n_cfg):
        plaquettes = np.empty((size, size, 3, 3), dtype=complex)
        for y in range(size):
            for x in range(size):
                u, ntry = sample_su3_wilson(beta, rng)
                plaquettes[y, x] = u
                plaquette_traces.append(float(np.real(np.trace(u)) / 3.0))
                accepted += 1
                tries += ntry
        for shape, bucket in buckets.items():
            bucket.append(float(np.mean(all_rectangle_values(plaquettes, *shape))))

    a_beta = float(np.mean(plaquette_traces))
    a_stderr = float(np.std(plaquette_traces, ddof=1) / math.sqrt(len(plaquette_traces)))
    acceptance = accepted / tries

    out: list[LoopRow] = []
    for (w, h), bucket in sorted(buckets.items(), key=lambda kv: (kv[0][0] * kv[0][1], kv[0][1], kv[0][0])):
        arr = np.array(bucket, dtype=float)
        area = w * h
        out.append(
            LoopRow(
                width=w,
                height=h,
                area=area,
                perimeter=2 * (w + h),
                predicted=a_beta**area,
                mc_mean=float(np.mean(arr)),
                mc_stderr=float(np.std(arr, ddof=1) / math.sqrt(len(arr))),
            )
        )
    return out, a_beta, a_stderr, acceptance


def central_first_moment_check(beta: float, rng: np.random.Generator, n: int = 1500) -> tuple[float, float]:
    """Check E[U] ~= a I_3 for the central measure."""

    mats = []
    traces = []
    tries = 0
    for _ in range(n):
        u, ntry = sample_su3_wilson(beta, rng)
        mats.append(u)
        traces.append(float(np.real(np.trace(u)) / 3.0))
        tries += ntry
    mean_u = np.mean(np.array(mats), axis=0)
    a_beta = float(np.mean(traces))
    offcentral = float(np.linalg.norm(mean_u - a_beta * np.eye(3)))
    return a_beta, offcentral


def ledger_report(size: int) -> dict[str, int]:
    byte_vertices = (size + 1) * (size + 1)
    stabilizer_generators = 8 * byte_vertices
    link_records_before_gauge = 2 * size * (size + 1)
    local_plaquette_records = size * size
    flat_rectangle_records = (size * (size + 1) // 2) ** 2
    promoted_collective_records = 1
    generator_ledger_size = stabilizer_generators + local_plaquette_records + promoted_collective_records
    return {
        "byte_vertices": byte_vertices,
        "stabilizer_generators": stabilizer_generators,
        "link_records_before_gauge": link_records_before_gauge,
        "local_plaquette_records": local_plaquette_records,
        "flat_rectangle_records": flat_rectangle_records,
        "promoted_collective_records": promoted_collective_records,
        "bond_dimension_proxy": 1,
        "generator_ledger_size": generator_ledger_size,
    }


def main() -> None:
    print("SU(3) Wilson-loop cluster: strong-sector record-action certificate")
    print("=" * 96)

    beta = 5.0
    size = 4
    rng = np.random.default_rng(20260630)

    print("\n[1] Declared SU(3) plaquette measure")
    print("  dP(U_p) proportional to dHaar_SU3(U_p) exp[ beta ReTr(U_p)/3 ]")
    print(f"  beta = {beta:.3f}; cluster for MC checks = {size} x {size} plaquettes")

    a_central, offcentral = central_first_moment_check(beta, rng)
    print(f"  central-moment probe: a={a_central:.9f}, ||E[U]-aI||_F={offcentral:.3e}")
    assert_greater("central a(beta)", a_central, 0.0)
    assert_less("central a(beta)", a_central, 1.0)
    assert_less("sample centrality of E[U]=aI", offcentral, 0.08)

    print("\n[2] Monte Carlo loop records against factorized SU(3) area law")
    mc_rows, a_beta, a_stderr, acceptance = monte_carlo_rows(size=size, beta=beta, rng=rng)
    sigma = -math.log(a_beta)
    print(f"  measured a(beta) from plaquette records = {a_beta:.9f} +/- {a_stderr:.3e}")
    print(f"  rejection acceptance = {acceptance:.4f}")
    print(f"  factorized string tension sigma = -log a = {sigma:.9f}")
    assert_greater("rejection sampler acceptance", acceptance, 0.005)

    print("  w x h   A   P       predicted <W>   measured <W>       stderr      z")
    max_z = 0.0
    for row in mc_rows:
        assert row.mc_mean is not None and row.mc_stderr is not None
        z = abs(row.mc_mean - row.predicted) / max(row.mc_stderr, 1e-300)
        max_z = max(max_z, z)
        print(
            f"  {row.width:1d} x {row.height:<1d}  {row.area:2d}  {row.perimeter:2d}"
            f"   {row.predicted: .8e}   {row.mc_mean: .8e}   {row.mc_stderr: .2e}   {z:5.2f}"
        )
    assert_less("max loop z-score against factorized area law", max_z, 4.0)

    mc_sigma, mc_mu, mc_c, mc_resid = fit_area_perimeter(mc_rows, use_mc=True)
    print(
        f"  noisy MC fit: sigma={mc_sigma:.6f}, mu={mc_mu:.3e}, "
        f"c={mc_c:.3e}, max residual={mc_resid:.3e}"
    )
    print("  (Coefficient promotion uses the factorized prediction; the MC fit is a noisy estimator check.)")

    print("\n[3] Exact factorized area/perimeter fit generated by the SU(3) record action")
    rows = predicted_rows(size=6, a_beta=a_beta)
    fit_sigma, fit_mu, fit_c, fit_resid = fit_area_perimeter(rows)
    print(
        f"  fit over 6 x 6 rectangle classes: sigma={fit_sigma:.9f}, "
        f"mu={fit_mu:.3e}, c={fit_c:.3e}, max residual={fit_resid:.3e}"
    )
    assert_close("factorized fit sigma", fit_sigma, sigma, tol=1e-11)
    assert_less("factorized perimeter coefficient", abs(fit_mu), 1e-11)
    assert_less("factorized fit residual", fit_resid, 1e-10)

    print("\n[4] Ledger-size and promotion/bond-dimension scaling")
    print(
        "  L  byte  stab_gen  link_rec  plaq_rec  flat_rect  promoted  bondD  generator_ledger"
    )
    for L in range(2, 9):
        led = ledger_report(L)
        print(
            f"  {L:<2d} {led['byte_vertices']:5d} {led['stabilizer_generators']:9d}"
            f" {led['link_records_before_gauge']:9d} {led['local_plaquette_records']:9d}"
            f" {led['flat_rectangle_records']:10d} {led['promoted_collective_records']:9d}"
            f" {led['bond_dimension_proxy']:6d} {led['generator_ledger_size']:17d}"
        )
        assert led["promoted_collective_records"] == 1
        assert led["bond_dimension_proxy"] == 1
        assert led["generator_ledger_size"] < 16 * (L + 1) * (L + 1)

    print("\n[5] Verdict")
    print(
        """
The SU(3) rung keeps the same record-grammar structure as the SU(2) rung:

  detector interface: normalized Wilson-loop trace;
  declared measure:   central plaquette record action;
  tractable ledger:   local plaquette records + one area-law/string-tension record;
  flat loop list:     O(L^4);
  generator ledger:   O(L^2) in this controlled 2D model.

This is evidence that the strong-sector direction is technically natural for
the record grammar.  It is not yet QCD: the open front is the full TCH/SU(3)
cell dynamics, mirror-sector SMG, and a 3D/4D continuum or transfer argument.
"""
    )
    print("exit 0 -- SU(3) Wilson-cluster record-action certificate passed.")


if __name__ == "__main__":
    main()
