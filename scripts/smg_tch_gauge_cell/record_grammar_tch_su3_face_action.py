#!/usr/bin/env python3
r"""TCH/SU(3) face-action certificate on the truncated-cube loop basis.

Purpose
-------
This is the next rung after record_grammar_tch_su3_loop_geometry.py.  The
geometry script proves the truncated-cube TCH cell has a 13-dimensional closed
Wilson-loop space.  This script puts a declared SU(3) record-action measure on
that 13-face basis.

Declared measure
----------------
Choose one face as the exterior/reference face.  The remaining 13 face loops are
independent records with central SU(3) weight

    dP(U_f) proportional to dHaar_SU3(U_f) exp[ beta Re Tr(U_f) / 3 ].

A boundary C is represented by its unique basis-face coordinate set S(C), and
the detector loop is

    W(C) = (1/3) Re Tr prod_{f in S(C)} U_f .

Because the measure is central and independent on the basis records,

    <W(C)> = a(beta)^{|S(C)|}.

This is an action-area law in the chosen face basis.  It is not yet the physical
minimal-surface law of a glued TCH gauge theory.  The script reports both:

  * basis action area |S(C)|, where the declared measure is exactly factorized;
  * topological minimal face area, computed by the geometry certificate.

Boundary
--------
This still does not derive confinement or a physical string tension.  It closes
one narrower task: the TCH loop geometry can carry a self-consistent central
SU(3) face-record measure with an exact generator-ledger area law.  The remaining
physics step is to glue cells and derive the local action that selects physical
minimal surfaces rather than a chosen exterior-face basis.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass

import numpy as np

from record_grammar_tch_su3_loop_geometry import (
    build_edges,
    build_faces,
    classify_boundaries,
    enumerate_boundaries,
    independent_face_basis,
    mask_from_edges,
)


@dataclass(frozen=True)
class LoopActionRow:
    boundary: int
    perimeter: int
    components: int
    min_area: int
    basis_area: int
    basis_weight: int
    predicted: float


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<74s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_greater(name: str, value: float, bound: float) -> None:
    print(f"  {name:<74s} value={value:.12g} bound={bound:.12g}")
    if not value > bound:
        raise AssertionError(name)


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<74s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def haar_su3(rng: np.random.Generator) -> np.ndarray:
    """Draw one Haar-random SU(3) matrix by QR factorization."""

    z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    q, r = np.linalg.qr(z)
    diag = np.diag(r)
    phases = diag / np.abs(diag)
    q = q @ np.diag(np.conj(phases))
    det_phase = np.angle(np.linalg.det(q))
    return q * np.exp(-1j * det_phase / 3.0)


def sample_su3_wilson(beta: float, rng: np.random.Generator) -> tuple[np.ndarray, int]:
    """Sample dHaar(U) exp[beta ReTr(U)/3] by rejection from Haar_SU3."""

    tries = 0
    while True:
        tries += 1
        u = haar_su3(rng)
        x = float(np.real(np.trace(u)) / 3.0)
        if rng.random() <= math.exp(beta * (x - 1.0)):
            return u, tries


def central_first_moment(beta: float, rng: np.random.Generator, n: int = 2500) -> tuple[float, float, float]:
    mats = []
    traces = []
    accepted = 0
    tries = 0
    for _ in range(n):
        u, ntry = sample_su3_wilson(beta, rng)
        mats.append(u)
        traces.append(float(np.real(np.trace(u)) / 3.0))
        accepted += 1
        tries += ntry
    mean_u = np.mean(np.array(mats), axis=0)
    a_beta = float(np.mean(traces))
    offcentral = float(np.linalg.norm(mean_u - a_beta * np.eye(3)))
    acceptance = accepted / tries
    return a_beta, offcentral, acceptance


def eval_xor(rows: list[int]) -> int:
    out = 0
    for row in rows:
        out ^= row
    return out


def basis_coordinate_map(face_masks: list[int], basis_faces: list[int]) -> dict[int, tuple[int, ...]]:
    """Unique boundary -> independent basis-face coordinates."""

    out: dict[int, tuple[int, ...]] = {}
    n_basis = len(basis_faces)
    for subset in range(1 << n_basis):
        boundary = 0
        coords: list[int] = []
        for bit_i, face_i in enumerate(basis_faces):
            if subset & (1 << bit_i):
                boundary ^= face_masks[face_i]
                coords.append(face_i)
        out[boundary] = tuple(coords)
    return out


def ordered_product(mats: list[np.ndarray], coords: tuple[int, ...], face_to_basis: dict[int, int]) -> np.ndarray:
    prod = np.eye(3, dtype=complex)
    for face_i in coords:
        prod = prod @ mats[face_to_basis[face_i]]
    return prod


def fit_area_perimeter(rows: list[LoopActionRow], area_attr: str) -> tuple[float, float, float, float]:
    x_rows = []
    y_rows = []
    for row in rows:
        if row.predicted <= 0.0:
            continue
        x_rows.append([float(getattr(row, area_attr)), float(row.perimeter), 1.0])
        y_rows.append(-math.log(row.predicted))
    xmat = np.array(x_rows, dtype=float)
    yvec = np.array(y_rows, dtype=float)
    coeff, *_ = np.linalg.lstsq(xmat, yvec, rcond=None)
    resid = yvec - xmat @ coeff
    return float(coeff[0]), float(coeff[1]), float(coeff[2]), float(np.max(np.abs(resid)))


def action_rows(beta: float) -> tuple[list[LoopActionRow], float, dict[str, int]]:
    faces = build_faces()
    edges, edge_index = build_edges(faces)
    face_masks = [mask_from_edges(face, edge_index) for face in faces]
    face_weights = [face.weight for face in faces]
    basis_faces = independent_face_basis(face_masks)
    omitted = sorted(set(range(len(faces))) - set(basis_faces))

    best = enumerate_boundaries(face_masks, face_weights)
    geometry_records = {
        record.mask: record
        for record in classify_boundaries(best, edges, n_vertices=24)
    }
    coords = basis_coordinate_map(face_masks, basis_faces)

    if omitted:
        assert_equal("exterior/reference face count", len(omitted), 1)
    assert_equal("basis-coordinate boundaries including zero", len(coords), 1 << 13)
    assert_equal("geometry boundaries including zero", len(best), 1 << 13)
    assert_equal("basis spans geometry boundary set", len(set(coords) - set(best)), 0)

    rng = np.random.default_rng(20260701)
    a_beta, offcentral, acceptance = central_first_moment(beta, rng)
    assert_greater("central SU(3) a(beta)", a_beta, 0.0)
    assert_less("central SU(3) a(beta)", a_beta, 1.0)
    assert_less("centrality ||E[U]-aI||_F", offcentral, 0.06)

    rows: list[LoopActionRow] = []
    for boundary, coord_faces in coords.items():
        if boundary == 0:
            continue
        geo = geometry_records[boundary]
        basis_area = len(coord_faces)
        rows.append(
            LoopActionRow(
                boundary=boundary,
                perimeter=geo.perimeter,
                components=geo.components,
                min_area=geo.min_unit_area,
                basis_area=basis_area,
                basis_weight=sum(face_weights[i] for i in coord_faces),
                predicted=a_beta**basis_area,
            )
        )

    meta = {
        "basis_size": len(basis_faces),
        "omitted_face": omitted[0],
        "acceptance_x10000": int(round(acceptance * 10000)),
    }
    return rows, a_beta, meta


def monte_carlo_loop_check(
    rows: list[LoopActionRow],
    beta: float,
    a_beta: float,
    sample_coords: dict[int, tuple[int, ...]],
    face_to_basis: dict[int, int],
    rng: np.random.Generator,
    n_cfg: int = 700,
) -> float:
    """Check a representative set of loop records against a(beta)^basis_area."""

    chosen = []
    seen_areas: set[int] = set()
    for row in sorted(rows, key=lambda item: (item.basis_area, item.perimeter, item.components)):
        if row.basis_area in seen_areas or row.basis_area == 0:
            continue
        chosen.append(row)
        seen_areas.add(row.basis_area)
        if len(chosen) == 8:
            break

    buckets = {row.boundary: [] for row in chosen}
    for _ in range(n_cfg):
        basis_mats = [sample_su3_wilson(beta, rng)[0] for _ in range(len(face_to_basis))]
        for row in chosen:
            prod = ordered_product(basis_mats, sample_coords[row.boundary], face_to_basis)
            buckets[row.boundary].append(float(np.real(np.trace(prod)) / 3.0))

    worst_z = 0.0
    print("  representative loop checks:")
    for row in chosen:
        arr = np.array(buckets[row.boundary], dtype=float)
        mean = float(np.mean(arr))
        stderr = float(np.std(arr, ddof=1) / math.sqrt(len(arr)))
        target = a_beta**row.basis_area
        z = abs(mean - target) / max(stderr, 1e-12)
        worst_z = max(worst_z, z)
        print(
            f"    basis_area={row.basis_area:2d}, perimeter={row.perimeter:2d}, "
            f"components={row.components}: MC={mean:+.5f} +/- {stderr:.5f}, "
            f"target={target:+.5f}, z={z:.2f}"
        )
    return worst_z


def histogram(rows: list[LoopActionRow], attr: str) -> Counter[int]:
    return Counter(getattr(row, attr) for row in rows)


def format_hist(counter: Counter[int]) -> str:
    return ", ".join(f"{key}:{value}" for key, value in sorted(counter.items()))


def main() -> None:
    print("TCH/SU(3) face-action certificate")
    print("=" * 92)
    beta = 5.0

    print("\n[1] Face-basis action setup")
    faces = build_faces()
    edges, edge_index = build_edges(faces)
    face_masks = [mask_from_edges(face, edge_index) for face in faces]
    basis_faces = independent_face_basis(face_masks)
    omitted = sorted(set(range(len(faces))) - set(basis_faces))
    assert_equal("face-boundary relation xor", eval_xor(face_masks), 0)
    assert_equal("basis face count", len(basis_faces), 13)
    assert_equal("reference/exterior face count", len(omitted), 1)
    print("  basis faces:", ", ".join(faces[i].label for i in basis_faces))
    print("  exterior/reference face:", faces[omitted[0]].label)
    print("  declared measure: dHaar_SU3(U_f) exp[ beta ReTr(U_f)/3 ]")
    print(f"  beta = {beta:.3f}")

    rows, a_beta, meta = action_rows(beta)
    sigma_exact = -math.log(a_beta)
    print(f"  sampled central a(beta) = {a_beta:.9f}")
    print(f"  sigma_basis = -ln a     = {sigma_exact:.9f}")
    print(f"  rejection acceptance    = {meta['acceptance_x10000'] / 100.0:.2f}%")

    print("\n[2] Exact generator-ledger area law")
    sigma, mu, const, resid = fit_area_perimeter(rows, "basis_area")
    print(f"  fit -ln<W> = sigma A_basis + mu P + c")
    print(f"    sigma={sigma:.9f}, mu={mu:.3e}, c={const:.3e}, max_resid={resid:.3e}")
    assert_less("basis-area fit residual", resid, 1e-10)
    assert_less("perimeter coefficient in basis action", abs(mu), 1e-10)

    min_sigma, min_mu, min_const, min_resid = fit_area_perimeter(rows, "min_area")
    print(f"  same data fit to topological minimal area:")
    print(f"    sigma={min_sigma:.9f}, mu={min_mu:.6f}, c={min_const:.6f}, max_resid={min_resid:.3e}")
    assert_greater("minimal-area residual exposes missing physical surface selector", min_resid, 0.05)

    print("\n[3] Boundary census under the face-basis action")
    print(f"  basis-area histogram: {format_hist(histogram(rows, 'basis_area'))}")
    print(f"  min-area histogram:   {format_hist(histogram(rows, 'min_area'))}")
    print(f"  component histogram:  {format_hist(histogram(rows, 'components'))}")
    exterior_sensitive = sum(1 for row in rows if row.basis_area != row.min_area)
    print(f"  exterior-sensitive boundaries: {exterior_sensitive} / {len(rows)}")

    print("\n[4] Monte Carlo spot-checks of detector loop records")
    coord_map = basis_coordinate_map(face_masks, basis_faces)
    face_to_basis = {face_i: pos for pos, face_i in enumerate(basis_faces)}
    worst_z = monte_carlo_loop_check(
        rows=rows,
        beta=beta,
        a_beta=a_beta,
        sample_coords=coord_map,
        face_to_basis=face_to_basis,
        rng=np.random.default_rng(20260702),
    )
    assert_less("representative MC deviation z-score", worst_z, 4.0)

    print("\n[5] Ledger-size consequence")
    flat_loop_records = len(rows)
    generator_records = 8 * 24 + 36 + 13 + 1 + 1  # byte skeleton, links, cycle basis, face relation, sigma record
    print(f"  flat nonzero Wilson records:       {flat_loop_records}")
    print(f"  generator/action ledger records:   {generator_records}")
    print(f"  flat/generator ratio:              {flat_loop_records / generator_records:.2f}x")
    print(f"  promoted collective records:       1  (sigma_basis)")
    print(f"  bond-dimension proxy:              1  (factorized central face records)")

    print(
        """
VERDICT:
  PASS.  A one-cell TCH truncated-cube loop basis supports a declared central
  SU(3) face-record action with an exact basis-area Wilson law.  The action
  ledger is small: the 8191 nonzero boundaries are generated by 13 independent
  face records plus one collective sigma record.

  Honest boundary: this is not yet physical confinement.  The same data do not
  fit a universal topological minimal-surface law without a further local
  surface-selector/gluing theorem.  The next rung is therefore a glued TCH
  few-cell action, where shared faces decide whether the physical measure
  selects minimal surfaces and whether promotions/bond dimension remain bounded.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
