"""Tong-Wang-You / Zeng-Zhu-Wang-You 3-4-5-0 validation model.

This module translates the published ITensor implementation linked from
arXiv:2202.12355 into the backend-independent term-list representation used by
this package. The MPS site order is the paper code's order:

    lattice orbital 0: flavors 3,4,5,0
    lattice orbital 1: flavors 3,4,5,0
    ...

One unit cell contains four lattice orbitals, so the full chain has 16*N
fermion modes for N unit cells.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi
from typing import Iterable, Literal

import numpy as np

from .charges import ChargeVector, Mode
from .terms import FermionOp, FermionTerm, add_hopping

FlavorName = Literal["3", "4", "5", "0"]
EdgeName = Literal["A", "B"]
MassKind = Literal["dirac", "majorana"]

FLAVOR_NAMES: tuple[FlavorName, ...] = ("3", "4", "5", "0")
FLAVOR_TO_INDEX: dict[FlavorName, int] = {name: i for i, name in enumerate(FLAVOR_NAMES)}

# The authors' ITensor code conserves two U(1) charges Q1,Q2, not total
# particle number. The physical 3-4-5-0 charge is Q1, but Q2 is an additional
# anomaly-free bookkeeping symmetry used by the reference implementation.
FLAVOR_CHARGES: dict[FlavorName, ChargeVector] = {
    "3": (3, 0),
    "4": (4, 5),
    "5": (5, 4),
    "0": (0, 3),
}


@dataclass(frozen=True)
class Model3450Spec:
    """Parameters for the 3-4-5-0 Chern-ladder validation model."""

    unit_cells: int = 20
    t1: float = 1.0
    t2: float = 0.5
    g1: float = 0.0
    g2: float = 0.0
    include_free: bool = True
    include_interaction: bool = True

    @property
    def lattice_orbitals(self) -> int:
        return 4 * self.unit_cells

    @property
    def edge_sites(self) -> int:
        return 2 * self.unit_cells

    @property
    def modes_per_orbital(self) -> int:
        return 4

    @property
    def total_modes(self) -> int:
        return self.lattice_orbitals * self.modes_per_orbital


def site_index_3450(orbital: int, flavor: int | FlavorName) -> int:
    """Return the zero-based fermion-mode index for a lattice orbital/flavor."""

    flavor_index = FLAVOR_TO_INDEX[flavor] if isinstance(flavor, str) else flavor
    if not 0 <= flavor_index < len(FLAVOR_NAMES):
        raise ValueError(f"flavor index out of range: {flavor_index}")
    return 4 * orbital + flavor_index


def orbital_from_edge(edge: EdgeName, edge_position: int) -> int:
    """Map an edge coordinate to the underlying ladder orbital."""

    if edge_position < 0:
        raise ValueError("edge_position must be non-negative")
    if edge == "A":
        return 2 * edge_position
    if edge == "B":
        return 2 * edge_position + 1
    raise ValueError(f"unknown edge {edge!r}")


def edge_site_index(edge: EdgeName, edge_position: int, flavor: int | FlavorName) -> int:
    return site_index_3450(orbital_from_edge(edge, edge_position), flavor)


def orbital_metadata(orbital: int) -> tuple[int, str, int]:
    """Return (unit_cell, edge, edge_position) for a ladder orbital."""

    unit_cell = orbital // 4
    sub = orbital % 4
    edge = "A" if sub in (0, 2) else "B"
    edge_position = orbital // 2
    return unit_cell, edge, edge_position


def build_3450_modes(spec: Model3450Spec) -> list[Mode]:
    if spec.unit_cells <= 0:
        raise ValueError("unit_cells must be positive")
    modes: list[Mode] = []
    for orbital in range(spec.lattice_orbitals):
        cell, edge, edge_position = orbital_metadata(orbital)
        for flavor_index, flavor_name in enumerate(FLAVOR_NAMES):
            modes.append(
                Mode(
                    index=site_index_3450(orbital, flavor_index),
                    cell=cell,
                    sector=f"edge_{edge}",
                    flavor=flavor_index,
                    charge=FLAVOR_CHARGES[flavor_name],
                )
            )
            # Store edge_position only implicitly to avoid changing the Mode API.
            _ = edge_position
    return modes


def build_3450_terms(spec: Model3450Spec) -> tuple[list[Mode], list[FermionTerm]]:
    modes = build_3450_modes(spec)
    terms: list[FermionTerm] = []
    if spec.include_free:
        build_3450_free_terms(spec, terms)
    if spec.include_interaction:
        build_3450_interaction_terms(spec, terms)
    return modes, terms


def build_3450_free_terms(spec: Model3450Spec, terms: list[FermionTerm]) -> None:
    """Add the Chern-ladder hopping terms from the published ITensor code."""

    phase = np.exp(1j * pi / 4)
    for unit in range(spec.unit_cells):
        a = 4 * unit
        b = 4 * unit + 1
        c = 4 * unit + 2
        d = 4 * unit + 3
        has_next = unit < spec.unit_cells - 1
        next_a = 4 * (unit + 1)
        next_b = 4 * (unit + 1) + 1

        for flavor_index, flavor_name in enumerate(FLAVOR_NAMES):
            # Flavors 3,4 have the Fig. 1 hopping pattern; flavors 5,0 use the
            # complex-conjugated pattern, giving the opposite Chern number.
            chern_sign = 1 if flavor_name in ("3", "4") else -1

            def p(power: int) -> complex:
                return spec.t1 * phase ** (chern_sign * power)

            def s(orbital: int) -> int:
                return site_index_3450(orbital, flavor_index)

            # Edge A, then edge B.
            add_hopping(terms, s(a), s(c), p(-1), f"3450.free.{flavor_name}.A.intra")
            if has_next:
                add_hopping(terms, s(c), s(next_a), p(+1), f"3450.free.{flavor_name}.A.next")

            add_hopping(terms, s(b), s(d), p(+1), f"3450.free.{flavor_name}.B.intra")
            if has_next:
                add_hopping(terms, s(d), s(next_b), p(-1), f"3450.free.{flavor_name}.B.next")

            # Inter-edge nearest and next-nearest links inside the unit cell.
            add_hopping(terms, s(a), s(d), -spec.t2, f"3450.free.{flavor_name}.diag.ad")
            add_hopping(terms, s(b), s(c), -spec.t2, f"3450.free.{flavor_name}.diag.bc")
            add_hopping(terms, s(a), s(b), p(+1), f"3450.free.{flavor_name}.rung.ab")
            add_hopping(terms, s(c), s(d), p(-1), f"3450.free.{flavor_name}.rung.cd")

            # Inter-cell diagonals.
            if has_next:
                add_hopping(terms, s(c), s(next_b), spec.t2, f"3450.free.{flavor_name}.diag.c_next_b")
                add_hopping(terms, s(d), s(next_a), spec.t2, f"3450.free.{flavor_name}.diag.d_next_a")


def six_fermion_term(
    coeff: complex,
    ops: Iterable[tuple[Literal["C", "Cd"], int]],
    label: str,
) -> FermionTerm:
    kind = {"C": "annihilate", "Cd": "create"}
    return FermionTerm(
        coeff=coeff,
        ops=tuple(FermionOp(kind[op], site) for op, site in ops),
        label=label,
    )


def build_3450_interaction_terms(spec: Model3450Spec, terms: list[FermionTerm]) -> None:
    """Add the Eq. (3) six-fermion mirror-edge interactions on edge B."""

    for edge_position in range(spec.edge_sites - 1):
        left = orbital_from_edge("B", edge_position)
        right = orbital_from_edge("B", edge_position + 1)

        def s(orbital: int, flavor: FlavorName) -> int:
            return site_index_3450(orbital, flavor)

        if spec.g1 != 0:
            term = six_fermion_term(
                spec.g1,
                (
                    ("C", s(left, "3")),
                    ("C", s(left, "5")),
                    ("Cd", s(left, "4")),
                    ("Cd", s(right, "4")),
                    ("C", s(left, "0")),
                    ("C", s(right, "0")),
                ),
                "3450.int.g1",
            )
            terms.append(term)
            terms.append(term.hc())

        if spec.g2 != 0:
            term = six_fermion_term(
                spec.g2,
                (
                    ("C", s(left, "3")),
                    ("C", s(right, "3")),
                    ("Cd", s(left, "5")),
                    ("Cd", s(right, "5")),
                    ("C", s(left, "4")),
                    ("C", s(left, "0")),
                ),
                "3450.int.g2",
            )
            terms.append(term)
            terms.append(term.hc())


def free_single_particle_matrix_3450(spec: Model3450Spec) -> np.ndarray:
    """Return the one-body matrix for H_free in the package MPS ordering."""

    free_spec = Model3450Spec(
        unit_cells=spec.unit_cells,
        t1=spec.t1,
        t2=spec.t2,
        g1=0.0,
        g2=0.0,
        include_free=True,
        include_interaction=False,
    )
    _, terms = build_3450_terms(free_spec)
    h = np.zeros((free_spec.total_modes, free_spec.total_modes), dtype=complex)
    for term in terms:
        if len(term.ops) != 2:
            continue
        left, right = term.ops
        if left.kind == "create" and right.kind == "annihilate":
            h[left.site, right.site] += term.coeff
    return h


def free_spectrum_3450(spec: Model3450Spec) -> tuple[np.ndarray, float, float]:
    """Return sorted one-body levels, Hermiticity residual, and half-fill gap."""

    h = free_single_particle_matrix_3450(spec)
    hermiticity = float(np.linalg.norm(h - h.conj().T))
    evals = np.linalg.eigvalsh(h)
    half = len(evals) // 2
    gap = float(evals[half] - evals[half - 1])
    return evals, hermiticity, gap


def initial_product_state_3450(spec: Model3450Spec) -> list[str]:
    """Half-filled product seed matching the published ITensor scripts."""

    state = ["empty" for _ in range(spec.total_modes)]
    for site in range(spec.total_modes):
        if site % 8 in (0, 1, 2, 3):
            state[site] = "full"
    return state


def single_fermion_correlation_term(
    edge: EdgeName,
    flavor: FlavorName,
    left_position: int,
    distance: int,
) -> list[tuple[str, int]]:
    left = edge_site_index(edge, left_position, flavor)
    right = edge_site_index(edge, left_position + distance, flavor)
    return [("Cd", left), ("C", right)]


def paper_edge_site_index(edge: EdgeName, unit_cell: int, flavor: int | FlavorName) -> int:
    """Site index used by the paper's correlation scripts.

    The published ITensor correlator code samples one boundary orbital per unit
    cell: orbital 4*x for edge A and orbital 4*x+1 for edge B.
    """

    if unit_cell < 0:
        raise ValueError("unit_cell must be non-negative")
    orbital = 4 * unit_cell + (0 if edge == "A" else 1)
    return site_index_3450(orbital, flavor)


def paper_single_fermion_correlation_term(
    edge: EdgeName,
    flavor: FlavorName,
    left_unit_cell: int,
    distance: int,
) -> list[tuple[str, int]]:
    left = paper_edge_site_index(edge, left_unit_cell, flavor)
    right = paper_edge_site_index(edge, left_unit_cell + distance, flavor)
    return [("Cd", left), ("C", right)]


MASS_PAIRS: tuple[tuple[FlavorName, FlavorName], ...] = (
    ("3", "5"),
    ("3", "0"),
    ("4", "5"),
    ("4", "0"),
)


def mass_correlation_term(
    edge: EdgeName,
    kind: MassKind,
    pair: tuple[FlavorName, FlavorName],
    left_position: int,
    distance: int,
) -> list[tuple[str, int]]:
    """Return a paper-style four-operator mass correlator term."""

    f_left, f_right = pair
    x = left_position
    y = left_position + distance
    xl = edge_site_index(edge, x, f_left)
    xr = edge_site_index(edge, x, f_right)
    yl = edge_site_index(edge, y, f_left)
    yr = edge_site_index(edge, y, f_right)
    if kind == "dirac":
        return [("Cd", xl), ("C", xr), ("Cd", yr), ("C", yl)]
    if kind == "majorana":
        return [("C", xl), ("C", xr), ("Cd", yr), ("Cd", yl)]
    raise ValueError(f"unknown mass kind {kind!r}")


def paper_mass_correlation_term(
    edge: EdgeName,
    kind: MassKind,
    pair: tuple[FlavorName, FlavorName],
    left_unit_cell: int,
    distance: int,
) -> list[tuple[str, int]]:
    f_left, f_right = pair
    x = left_unit_cell
    y = left_unit_cell + distance
    xl = paper_edge_site_index(edge, x, f_left)
    xr = paper_edge_site_index(edge, x, f_right)
    yl = paper_edge_site_index(edge, y, f_left)
    yr = paper_edge_site_index(edge, y, f_right)
    if kind == "dirac":
        return [("Cd", xl), ("C", xr), ("Cd", yr), ("C", yl)]
    if kind == "majorana":
        return [("C", xl), ("C", xr), ("Cd", yr), ("Cd", yl)]
    raise ValueError(f"unknown mass kind {kind!r}")


def centered_left_position(edge_sites: int, distance: int) -> int:
    if distance <= 0:
        raise ValueError("distance must be positive")
    if distance >= edge_sites:
        raise ValueError("distance must be smaller than the number of edge sites")
    return edge_sites // 2 - distance // 2 - 1


def centered_left_unit_cell(unit_cells: int, distance: int) -> int:
    if distance <= 0:
        raise ValueError("distance must be positive")
    if distance >= unit_cells:
        raise ValueError("distance must be smaller than the number of unit cells")
    return unit_cells // 2 - distance // 2 - 1
