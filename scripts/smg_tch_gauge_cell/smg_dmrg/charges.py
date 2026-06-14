"""Charge bookkeeping for the 1+1D SMG prototype."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

ChargeVector = tuple[int, ...]


@dataclass(frozen=True)
class Mode:
    """One spinless fermion mode in the MPS/Fock ordering."""

    index: int
    cell: int
    sector: str
    flavor: int
    charge: ChargeVector


def add_charge(a: ChargeVector, b: ChargeVector) -> ChargeVector:
    if len(a) != len(b):
        raise ValueError(f"charge rank mismatch: {len(a)} != {len(b)}")
    return tuple(x + y for x, y in zip(a, b))


def neg_charge(a: ChargeVector) -> ChargeVector:
    return tuple(-x for x in a)


def zero_charge(rank: int) -> ChargeVector:
    return tuple(0 for _ in range(rank))


def build_flavor_cartan_modes(
    cells: int,
    flavors: int,
    *,
    mirror_charge_sign: int = -1,
    include_total_number: bool = True,
) -> list[Mode]:
    """Build physical+mirror modes with simple abelian Cartan charges.

    This is a conservative DMRG blocking choice rather than a full non-abelian
    implementation. Each flavor carries one Cartan basis vector. Mirror charges
    default to the conjugate sign, matching the vector-like physical+mirror
    bookkeeping used in the feasibility scripts.
    """

    if cells <= 0:
        raise ValueError("cells must be positive")
    if flavors <= 0:
        raise ValueError("flavors must be positive")
    if mirror_charge_sign not in (-1, 1):
        raise ValueError("mirror_charge_sign must be -1 or +1")

    modes: list[Mode] = []
    per_cell = 2 * flavors
    for cell in range(cells):
        for sector_offset, sector in enumerate(("physical", "mirror")):
            sign = 1 if sector == "physical" else mirror_charge_sign
            for flavor in range(flavors):
                index = cell * per_cell + sector_offset * flavors + flavor
                cartan = [0 for _ in range(flavors)]
                cartan[flavor] = sign
                charge = tuple(([1] if include_total_number else []) + cartan)
                modes.append(Mode(index, cell, sector, flavor, charge))
    return modes


def charge_rank(modes: Iterable[Mode]) -> int:
    ranks = {len(mode.charge) for mode in modes}
    if len(ranks) != 1:
        raise ValueError(f"inconsistent charge ranks: {sorted(ranks)}")
    return ranks.pop()

