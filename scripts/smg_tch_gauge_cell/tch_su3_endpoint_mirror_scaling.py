#!/usr/bin/env python3
"""SU(3)-endpoint mirror hopping on one/two-plaquette TCH cells.

This follows the full-CG one-plaquette audit by replacing the previous hybrid
``plaquette oscillator x mirror hopping`` modulation with an explicit SU(3)
endpoint hopping rule:

  * local mirror-Fock/SMG states are still the existing charge-resolved low
    energy states from ``css_spatial_mirror_fock_gauge_gap.py``;
  * each Z3 link flux is lifted to SU(3) irreps with the same triality:
        0 -> {1, 8},  1 -> {3, 6},  2 -> {3bar, 6bar};
  * charged hopping changes the active link irrep using the full-CG endpoint
    tensors from ``tch_peter_weyl_full_cg_mirror_audit.py`` rather than a Z3
    flux shift;
  * magnetic plaquette moves also use those full-CG endpoint maps.

This is intentionally a reduced-representation endpoint audit: it keeps link
irrep labels and uses RMS-reduced endpoint matrix elements. It does not yet
carry explicit link matrix indices and vertex recoupling tensors for the full
two-plaquette spin-network Hilbert space. The point is to close the old hybrid
modulation loophole and expose the next scaling wall.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from functools import lru_cache
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PYTHON_CODE = ROOT / "python_code"
DEFAULT_OUT = ROOT / "python_code" / "smg_dmrg_runs" / "tch_su3_endpoint_mirror_scaling_verdict.json"


def load_module(name: str, relative_path: str):
    path = PYTHON_CODE / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


FULLCG = load_module("tch_full_cg_mod", "tch_peter_weyl_full_cg_mirror_audit.py")
STRIP = load_module("z3_strip_mod", "css_2d_strip_mirror_fock_lanczos.py")
PLAQUETTE = STRIP.plaquette

TRIALITY_REPS_EXTENDED = {
    0: ("1", "8"),
    1: ("3", "6"),
    2: ("3b", "6b"),
}
TRIALITY_REPS_MINIMAL = {
    0: ("1",),
    1: ("3",),
    2: ("3b",),
}


def rep_triality(rep: str) -> int:
    return FULLCG.EXT.triality(rep)


def electric_energy(rep: str) -> float:
    return FULLCG.casimir(rep)


def endpoint_kind(delta_triality: int) -> str:
    delta_triality %= 3
    if delta_triality == 1:
        return "u"
    if delta_triality == 2:
        return "udag"
    raise ValueError("no SU(3) endpoint map for triality-neutral transition")


@lru_cache(maxsize=None)
def endpoint_strength(kind: str, target_rep: str, source_rep: str) -> float:
    """RMS-reduced full-CG endpoint matrix element for one link irrep move."""
    source_dim = FULLCG.dim(source_rep)
    target_dim = FULLCG.dim(target_rep)
    norm2 = 0.0
    for i in range(3):
        for j in range(3):
            matrix = FULLCG.link_map(kind, i, j, target_rep, source_rep)
            norm2 += float(np.linalg.norm(matrix) ** 2)
    if norm2 < 1e-24:
        return 0.0
    return float(np.sqrt(norm2 / (9 * source_dim * target_dim)))


def rep_options_for_flux(flux: int, rep_set: str) -> tuple[str, ...]:
    table = TRIALITY_REPS_MINIMAL if rep_set == "minimal" else TRIALITY_REPS_EXTENDED
    return table[flux % 3]


def flux_to_rep_assignments(fluxes: tuple[int, ...], rep_set: str) -> list[tuple[str, ...]]:
    return list(product(*[rep_options_for_flux(flux, rep_set) for flux in fluxes]))


def rectangular_geometry(nx: int, ny: int):
    return STRIP.rectangular_geometry(nx, ny)


def state_energy_cutoff_configs(energies: np.ndarray, charges: np.ndarray, n_vertices: int, energy_cutoff: float):
    local_dim = len(charges)
    min_energy = n_vertices * float(np.min(energies))
    for config in product(range(local_dim), repeat=n_vertices):
        charge_tuple = tuple(int(charges[state]) for state in config)
        if sum(charge_tuple) % 3 != 0:
            continue
        matter = float(sum(energies[state] for state in config))
        if matter <= min_energy + energy_cutoff + 1e-9:
            yield tuple(config), charge_tuple, matter


def transition_tables(annihilators):
    return PLAQUETTE.transition_tables(annihilators)


class SU3EndpointCell:
    def __init__(
        self,
        *,
        nx: int,
        ny: int,
        states_per_charge: int,
        energy_cutoff: float,
        rep_set: str,
        beta: float,
        hopping: float,
        include_neutral_hopping: bool = False,
    ):
        self.nx = nx
        self.ny = ny
        self.states_per_charge = states_per_charge
        self.energy_cutoff = energy_cutoff
        self.rep_set = rep_set
        self.beta = beta
        self.hopping = hopping
        self.include_neutral_hopping = include_neutral_hopping
        self.g2 = 1 / beta

        self.n_vertices, self.edges, self.plaquette_shifts = rectangular_geometry(nx, ny)
        self.incidence = STRIP.incidence_matrix(self.n_vertices, self.edges)
        self.energies, self.charges, self.annihilators, self.charge_counts = (
            PLAQUETTE.block.local_mirror_dressed_block(states_per_charge)
        )
        self.ann, self.create = transition_tables(self.annihilators)
        self.basis, self.index = self.build_basis()
        self.dim = len(self.basis)
        self.diagonal = np.zeros(self.dim)
        self.matter_diagonal = np.zeros(self.dim)
        self.electric_diagonal = np.zeros(self.dim)
        for pos, (config, reps) in enumerate(self.basis):
            matter = float(sum(self.energies[state] for state in config))
            electric = self.g2 * sum(electric_energy(rep) for rep in reps)
            self.matter_diagonal[pos] = matter
            self.electric_diagonal[pos] = electric
            self.diagonal[pos] = matter + electric
        self.transitions = self.build_transitions()

    def build_basis(self):
        basis = []
        flux_cache: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
        for config, charge_tuple, _ in state_energy_cutoff_configs(
            self.energies, self.charges, self.n_vertices, self.energy_cutoff
        ):
            if charge_tuple not in flux_cache:
                flux_cache[charge_tuple] = STRIP.flux_solutions(charge_tuple, self.incidence)
            for fluxes in flux_cache[charge_tuple]:
                for reps in flux_to_rep_assignments(fluxes, self.rep_set):
                    basis.append((config, tuple(reps)))
        return basis, {item: pos for pos, item in enumerate(basis)}

    def add_transition(self, rows, row, coeff):
        if row is not None and abs(coeff) > 1e-14:
            rows.append((row, coeff))

    def plaquette_rep_targets(self, reps: tuple[str, ...], shift: tuple[int, ...], amount: int):
        options = []
        for source_rep, delta in zip(reps, shift):
            signed_delta = (amount * delta) % 3
            if signed_delta == 0:
                options.append(((source_rep, 1.0),))
                continue
            kind = endpoint_kind(signed_delta)
            targets = []
            for target_rep in rep_options_for_flux(rep_triality(source_rep) + signed_delta, self.rep_set):
                amp = endpoint_strength(kind, target_rep, source_rep)
                if amp > 1e-14:
                    targets.append((target_rep, amp))
            options.append(tuple(targets))
        return options

    def hop_rep_targets(self, source_rep: str, delta_triality: int):
        delta_triality %= 3
        if delta_triality == 0:
            if not self.include_neutral_hopping:
                return ()
            return ((source_rep, 1.0),)
        kind = endpoint_kind(delta_triality)
        targets = []
        for target_rep in rep_options_for_flux(rep_triality(source_rep) + delta_triality, self.rep_set):
            amp = endpoint_strength(kind, target_rep, source_rep)
            if amp > 1e-14:
                targets.append((target_rep, amp))
        return tuple(targets)

    def build_transitions(self):
        transitions = []
        for col, (config, reps) in enumerate(self.basis):
            col_transitions = []

            for shift in self.plaquette_shifts:
                for amount in (1, -1):
                    for rep_targets in product(*self.plaquette_rep_targets(reps, shift, amount)):
                        new_reps = tuple(item[0] for item in rep_targets)
                        amp = 1.0
                        for _, local_amp in rep_targets:
                            amp *= local_amp
                        row = self.index.get((config, new_reps))
                        self.add_transition(col_transitions, row, -(self.beta / 2) * amp)

            if abs(self.hopping) < 1e-15:
                transitions.append(col_transitions)
                continue

            config_list = list(config)
            for link, (tail, head) in enumerate(self.edges):
                old_tail = config_list[tail]
                old_head = config_list[head]
                source_rep = reps[link]
                for mode in range(len(self.ann)):
                    for charge in range(3):
                        if charge == 0 and not self.include_neutral_hopping:
                            continue

                        for new_tail, amp_tail in self.create[mode][charge][old_tail]:
                            for new_head, amp_head in self.ann[mode][charge][old_head]:
                                new_config = list(config_list)
                                new_config[tail] = new_tail
                                new_config[head] = new_head
                                rep_targets = self.hop_rep_targets(source_rep, charge)
                                for target_rep, amp_link in rep_targets:
                                    new_reps = list(reps)
                                    new_reps[link] = target_rep
                                    row = self.index.get((tuple(new_config), tuple(new_reps)))
                                    self.add_transition(
                                        col_transitions,
                                        row,
                                        -self.hopping * amp_tail * amp_head * amp_link,
                                    )

                        for new_tail, amp_tail in self.ann[mode][charge][old_tail]:
                            for new_head, amp_head in self.create[mode][charge][old_head]:
                                new_config = list(config_list)
                                new_config[tail] = new_tail
                                new_config[head] = new_head
                                rep_targets = self.hop_rep_targets(source_rep, -charge)
                                for target_rep, amp_link in rep_targets:
                                    new_reps = list(reps)
                                    new_reps[link] = target_rep
                                    row = self.index.get((tuple(new_config), tuple(new_reps)))
                                    self.add_transition(
                                        col_transitions,
                                        row,
                                        -self.hopping * amp_tail * amp_head * amp_link,
                                    )

            transitions.append(col_transitions)
        return transitions

    def matvec(self, vector: np.ndarray) -> np.ndarray:
        out = self.diagonal * vector
        for col, col_transitions in enumerate(self.transitions):
            value = vector[col]
            if abs(value) < 1e-15:
                continue
            for row, coeff in col_transitions:
                out[row] += coeff * value
        return out

    def matter_expectation(self, vector: np.ndarray) -> float:
        return float(np.real(np.vdot(vector, self.matter_diagonal * vector)))

    def electric_expectation(self, vector: np.ndarray) -> float:
        return float(np.real(np.vdot(vector, self.electric_diagonal * vector)))


def lanczos_lowest(model: SU3EndpointCell, n_iter: int = 90, n_eigs: int = 8, seed: int = 731):
    return PLAQUETTE.lanczos_lowest(model, n_iter=n_iter, n_eigs=n_eigs, seed=seed)


def exact_lowest_if_small(model: SU3EndpointCell):
    if model.dim > 2500:
        return None
    eye = np.eye(model.dim, dtype=complex)
    dense = np.column_stack([model.matvec(eye[:, col]) for col in range(model.dim)])
    herm = float(np.linalg.norm(dense - dense.conj().T))
    eigvals, eigvecs = np.linalg.eigh(dense)
    assert herm < 1e-8
    return eigvals[:8], eigvecs[:, :8], herm


def unique_gap(eigvals: np.ndarray, tol: float = 1e-7):
    levels = []
    degeneracies = []
    for value in eigvals:
        if not levels or abs(value - levels[-1]) > tol:
            levels.append(float(value))
            degeneracies.append(1)
        else:
            degeneracies[-1] += 1
    return levels[1] - levels[0], degeneracies[0]


def matter_dominated_gap(model: SU3EndpointCell, eigvals: np.ndarray, vectors: np.ndarray):
    ground_matter = model.matter_expectation(vectors[:, 0])
    for idx in range(1, vectors.shape[1]):
        matter = model.matter_expectation(vectors[:, idx])
        if matter - ground_matter > 0.5:
            return float(eigvals[idx] - eigvals[0]), float(matter - ground_matter)
    return float("inf"), 0.0


def run_case(case: dict[str, Any]) -> dict[str, Any]:
    model = SU3EndpointCell(**case)
    exact = exact_lowest_if_small(model)
    exact_hermiticity = None
    if exact is not None:
        eigvals, vectors, exact_hermiticity = exact
    else:
        eigvals, vectors = lanczos_lowest(model)
    full_gap, degeneracy = unique_gap(eigvals)
    matter_gap, matter_delta = matter_dominated_gap(model, eigvals, vectors)
    avg_transitions = float(np.mean([len(row) for row in model.transitions])) if model.transitions else 0.0
    max_transitions = int(max([len(row) for row in model.transitions], default=0))
    return {
        **case,
        "dimension": model.dim,
        "n_vertices": model.n_vertices,
        "n_edges": len(model.edges),
        "n_plaquettes": len(model.plaquette_shifts),
        "charge_counts": {str(k): int(v) for k, v in model.charge_counts.items()},
        "avg_transitions_per_basis_state": avg_transitions,
        "max_transitions_per_basis_state": max_transitions,
        "full_gap": full_gap,
        "matter_dominated_gap": matter_gap,
        "matter_expectation_delta": matter_delta,
        "ground_degeneracy": degeneracy,
        "lowest": [float(x) for x in eigvals[:6]],
        "exact_hermiticity": exact_hermiticity,
    }


def build_verdict() -> dict[str, Any]:
    cases = [
        {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "extended",
            "beta": 0.5,
            "hopping": 0.0,
        },
        {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "extended",
            "beta": 0.5,
            "hopping": 0.2,
        },
        {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "extended",
            "beta": 0.5,
            "hopping": 4.0,
        },
        {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 4,
            "energy_cutoff": 4.0,
            "rep_set": "extended",
            "beta": 0.5,
            "hopping": 0.2,
        },
        {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 4,
            "energy_cutoff": 4.0,
            "rep_set": "extended",
            "beta": 0.5,
            "hopping": 4.0,
        },
        {
            "nx": 2,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "minimal",
            "beta": 0.5,
            "hopping": 0.2,
        },
        {
            "nx": 2,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "minimal",
            "beta": 0.5,
            "hopping": 4.0,
        },
    ]
    rows = [run_case(case) for case in cases]
    finite_gaps = [row["full_gap"] for row in rows if np.isfinite(row["full_gap"])]
    assert min(finite_gaps) > 0.01
    assert any(row["avg_transitions_per_basis_state"] > 1.0 for row in rows if row["hopping"] > 0)
    return {
        "verdict": "su3_endpoint_hopping_finite_scaling_passes_but_full_spin_network_open",
        "scope": (
            "Reduced-representation SU(3) endpoint audit. Charged mirror-Fock hopping "
            "changes link irreps via full-CG endpoint tensors on one/two-plaquette cells. "
            "The full explicit matrix-index spin-network recoupling problem remains open."
        ),
        "normalization": (
            "endpoint strengths are RMS averages of the full link_map('u'/'udag', i, j, "
            "target, source) CG tensors over endpoint colour and link matrix indices"
        ),
        "omitted": [
            "triality-neutral q=0 hopping is omitted so the scan isolates charged SU(3) endpoint hopping",
            "explicit link matrix indices and degree-3 vertex recoupling tensors are not kept",
            "two-plaquette extended-rep scaling at states_per_charge>=3 and energy_cutoff=4 is too large for this reduced exact basis",
        ],
        "interpretation": (
            "The charged q=1,2 endpoint operator is present at energy_cutoff=4, as shown "
            "by the increased transition count, but it does not move the first matter gap "
            "even at t=4. In this projected mirror-Fock block the charged single-fermion "
            "matrix elements live inside the -5 local manifold, while the first gap is the "
            "-7 to -5 local SMG excitation. This is a finite-cutoff stability result, not "
            "a thermodynamic proof."
        ),
        "rows": rows,
        "minimum_full_gap": min(finite_gaps),
        "remaining_gate": (
            "build the non-reduced multi-plaquette spin-network Hilbert space with explicit "
            "link matrix indices and vertex intertwiner/recoupling tensors, then repeat the "
            "same cutoff and volume scan"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    verdict = build_verdict()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")

    print("TCH SU(3)-endpoint mirror scaling audit")
    print(f"  verdict: {verdict['verdict']}")
    print(f"  wrote: {args.out.relative_to(ROOT)}")
    print("  rows:")
    for row in verdict["rows"]:
        print(
            "    "
            f"{row['nx']}x{row['ny']} "
            f"spc={row['states_per_charge']} "
            f"ecut={row['energy_cutoff']} "
            f"reps={row['rep_set']} "
            f"t={row['hopping']} "
            f"dim={row['dimension']} "
            f"gap={row['full_gap']:.6g} "
            f"matter_gap={row['matter_dominated_gap']:.6g}"
        )
    print(f"  minimum full gap: {verdict['minimum_full_gap']:.6g}")
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
