"""Model construction for the SMG DMRG prototype."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from .charges import ChargeVector, Mode, add_charge, build_flavor_cartan_modes
from .terms import FermionOp, FermionTerm, add_hopping, add_onsite_number


@dataclass(frozen=True)
class ModelSpec:
    """Parameters for a 1+1D physical+mirror chain."""

    cells: int = 2
    flavors: int = 4
    t: float = 1.0
    r: float = 1.0
    mass: float = 0.0
    u: float = 1.0
    mirror_charge_sign: int = -1
    opposite_wilson_for_mirror: bool = True
    quartic_sign: float = 1.0
    normalize_pair_groups: bool = True

    @property
    def modes_per_cell(self) -> int:
        return 2 * self.flavors

    @property
    def total_modes(self) -> int:
        return self.cells * self.modes_per_cell


def mode_index(spec: ModelSpec, cell: int, sector: str, flavor: int) -> int:
    if sector not in ("physical", "mirror"):
        raise ValueError(f"unknown sector {sector!r}")
    sector_offset = 0 if sector == "physical" else 1
    return cell * spec.modes_per_cell + sector_offset * spec.flavors + flavor


def build_model_terms(spec: ModelSpec) -> tuple[list[Mode], list[FermionTerm]]:
    modes = build_flavor_cartan_modes(
        spec.cells,
        spec.flavors,
        mirror_charge_sign=spec.mirror_charge_sign,
    )
    terms: list[FermionTerm] = []
    build_wilson_mirror_terms(spec, terms)
    build_charge_grouped_pair_smg_terms(spec, modes, terms)
    return modes, terms


def build_wilson_mirror_terms(spec: ModelSpec, terms: list[FermionTerm]) -> None:
    """Add the isolated Wilson/mirror kinetic convention.

    This is the convention-localised part of the build. It follows the scope
    file's schematic nearest-neighbour Wilson expression. The exact convention
    should be replaced once the Phase-1 validation reference is selected.
    """

    for cell in range(spec.cells):
        for sector in ("physical", "mirror"):
            wilson_sign = -1.0 if (
                sector == "mirror" and spec.opposite_wilson_for_mirror
            ) else 1.0
            for flavor in range(spec.flavors):
                site = mode_index(spec, cell, sector, flavor)
                onsite = spec.mass + wilson_sign * spec.r
                add_onsite_number(terms, site, onsite, f"{sector}.onsite")

    for cell in range(spec.cells - 1):
        for sector in ("physical", "mirror"):
            wilson_sign = -1.0 if (
                sector == "mirror" and spec.opposite_wilson_for_mirror
            ) else 1.0
            hop = -0.5 * spec.t - 0.5 * wilson_sign * spec.r
            for flavor in range(spec.flavors):
                i = mode_index(spec, cell, sector, flavor)
                j = mode_index(spec, cell + 1, sector, flavor)
                add_hopping(terms, i, j, hop, f"{sector}.wilson_hop")


def build_charge_grouped_pair_smg_terms(
    spec: ModelSpec,
    modes: list[Mode],
    terms: list[FermionTerm],
) -> None:
    """Add a Cartan-conserving antisymmetric mirror-sector quartic.

    Pair annihilators are grouped by their total abelian charge. Within each
    charge sector we add B_Q^dagger B_Q, with

        B_Q = sum_{i<j, q_i+q_j=Q} c_i c_j.

    This is a generic Cartan-neutral antisymmetric quartic scaffold. It is not
    the final SO(10) 120-channel Clebsch-Gordan tensor.
    """

    if spec.u == 0:
        return

    mode_by_site = {mode.index: mode for mode in modes}
    for cell in range(spec.cells):
        mirror_sites = [
            mode_index(spec, cell, "mirror", flavor) for flavor in range(spec.flavors)
        ]
        groups: dict[ChargeVector, list[tuple[int, int]]] = defaultdict(list)
        for pos, i in enumerate(mirror_sites):
            for j in mirror_sites[pos + 1 :]:
                pair_charge = add_charge(mode_by_site[i].charge, mode_by_site[j].charge)
                groups[pair_charge].append((i, j))

        for group_index, pairs in enumerate(groups.values()):
            norm = len(pairs) if spec.normalize_pair_groups else 1
            coeff = spec.quartic_sign * spec.u / norm
            for i, j in pairs:
                for k, l in pairs:
                    label = f"mirror.smg_pair_group_{group_index}"
                    if (i, j) == (k, l):
                        # c_i^dag c_j^dag c_j c_i = n_i n_j for i < j.
                        # Emitting the density form avoids repeated single-site
                        # operators in the MPO backend.
                        terms.append(
                            FermionTerm(
                                coeff,
                                (FermionOp("number", i), FermionOp("number", j)),
                                label,
                            )
                        )
                    else:
                        terms.append(
                            FermionTerm(
                                coeff,
                                (
                                    FermionOp("create", i),
                                    FermionOp("create", j),
                                    FermionOp("annihilate", l),
                                    FermionOp("annihilate", k),
                                ),
                                label,
                            )
                        )
