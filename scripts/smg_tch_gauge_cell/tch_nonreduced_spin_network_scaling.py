#!/usr/bin/env python3
"""Non-reduced SU(3) spin-network mirror-hopping audit.

This is the next finite artifact after ``tch_su3_endpoint_mirror_scaling.py``.
The reduced endpoint audit kept only link irrep labels and RMS-reduced CG
endpoint strengths. Here the basis carries explicit vertex-intertwiner labels
and matrix elements are computed by tensor contractions over:

  * local mirror-Fock matter CG maps,
  * full-CG SU(3) link endpoint maps, and
  * vertex singlet/recoupling tensors for the actual graph degree.

The link matrix indices are not averaged into a scalar endpoint strength. They
are contracted explicitly in the local edge operator. This is still finite and
truncated: local mirror states are represented by their minimal SU(3) charge
lift q=0,1,2 -> 1,3,3bar, and the scan is small enough for exact/Krylov ED.
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
DEFAULT_OUT = ROOT / "python_code" / "smg_dmrg_runs" / "tch_nonreduced_spin_network_scaling_verdict.json"


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

REPS_MINIMAL = ("1", "3", "3b")
REPS_EXTENDED = ("1", "3", "3b", "6", "6b", "8")
MATTER_BY_CHARGE = {0: "1", 1: "3", 2: "3b"}


def rep_triality(rep: str) -> int:
    return FULLCG.EXT.triality(rep)


def dual_rep(rep: str) -> str:
    return FULLCG.dual(rep)


def rep_dim(rep: str) -> int:
    return FULLCG.dim(rep)


def matter_rep(charge: int) -> str:
    return MATTER_BY_CHARGE[int(charge) % 3]


def rep_options_for_flux(flux: int, rep_set: str) -> tuple[str, ...]:
    if rep_set == "minimal":
        table = {0: ("1",), 1: ("3",), 2: ("3b",)}
    else:
        table = {0: ("1", "8"), 1: ("3", "6"), 2: ("3b", "6b")}
    return table[flux % 3]


def rectangular_geometry(nx: int, ny: int):
    return STRIP.rectangular_geometry(nx, ny)


@lru_cache(maxsize=None)
def invariant_basis(rep_tuple: tuple[str, ...]) -> tuple[np.ndarray, ...]:
    dims = [rep_dim(rep) for rep in rep_tuple]
    total_dim = int(np.prod(dims))
    blocks = FULLCG.generator_sum(rep_tuple)
    stacked = np.vstack(blocks)
    _, singular, vh = np.linalg.svd(stacked)
    nullity = int(np.sum(singular < 1e-9))
    if nullity == 0:
        return ()
    vecs = vh.conj().T[:, -nullity:]
    # SVD already gives an orthonormal nullspace. Fix phases for deterministic JSON.
    tensors = []
    for col in range(nullity):
        vec = vecs[:, col]
        pivot = int(np.argmax(np.abs(vec)))
        if abs(vec[pivot]) > 0:
            vec = vec * np.conj(vec[pivot]) / abs(vec[pivot])
        tensors.append(vec.reshape(tuple(dims)))
    return tuple(tensors)


@lru_cache(maxsize=None)
def matter_component(factor_rep: str, source_rep: str, target_rep: str, color: int) -> np.ndarray:
    if target_rep not in FULLCG.fuse_targets(factor_rep, source_rep):
        return np.zeros((rep_dim(target_rep), rep_dim(source_rep)), complex)
    embedding = FULLCG.cg_embedding(factor_rep, source_rep, target_rep)
    source_dim = rep_dim(source_rep)
    target_dim = rep_dim(target_rep)
    tensor = embedding.reshape((rep_dim(factor_rep), source_dim, target_dim), order="C")
    return tensor[color, :, :].T.copy()


def charge_factor(delta_charge: int) -> str:
    delta_charge %= 3
    if delta_charge == 1:
        return "3"
    if delta_charge == 2:
        return "3b"
    raise ValueError("neutral matter transitions are omitted in this audit")


def link_kind(delta_triality: int) -> str:
    delta_triality %= 3
    if delta_triality == 1:
        return "u"
    if delta_triality == 2:
        return "udag"
    raise ValueError("neutral link transitions are omitted in this audit")


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


class SpinNetworkCell:
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
        n_lanczos: int = 80,
        progress_every: int = 0,
    ):
        self.nx = nx
        self.ny = ny
        self.states_per_charge = states_per_charge
        self.energy_cutoff = energy_cutoff
        self.rep_set = rep_set
        self.beta = beta
        self.hopping = hopping
        self.n_lanczos = n_lanczos
        self.progress_every = progress_every
        self.g2 = 1 / beta

        self.n_vertices, self.edges, self.plaquette_shifts = rectangular_geometry(nx, ny)
        self.incidence = STRIP.incidence_matrix(self.n_vertices, self.edges)
        self.vertex_edges = self.build_vertex_edges()
        self.energies, self.charges, self.annihilators, self.charge_counts = (
            PLAQUETTE.block.local_mirror_dressed_block(states_per_charge)
        )
        self.ann, self.create = PLAQUETTE.transition_tables(self.annihilators)
        self.basis, self.index = self.build_basis()
        self.dim = len(self.basis)
        if self.progress_every:
            print(
                f"  built basis: nx={nx} ny={ny} spc={states_per_charge} "
                f"ecut={energy_cutoff} reps={rep_set} dim={self.dim}",
                flush=True,
            )
        self.diagonal = np.zeros(self.dim)
        self.matter_diagonal = np.zeros(self.dim)
        self.electric_diagonal = np.zeros(self.dim)
        for pos, state in enumerate(self.basis):
            config, reps, _ = state
            matter = float(sum(self.energies[site] for site in config))
            electric = self.g2 * sum(FULLCG.casimir(rep) for rep in reps)
            self.matter_diagonal[pos] = matter
            self.electric_diagonal[pos] = electric
            self.diagonal[pos] = matter + electric
        self.transitions = self.build_transitions()

    def build_vertex_edges(self) -> list[list[tuple[int, str]]]:
        rows: list[list[tuple[int, str]]] = [[] for _ in range(self.n_vertices)]
        for edge, (tail, head) in enumerate(self.edges):
            rows[tail].append((edge, "tail"))
            rows[head].append((edge, "head"))
        return rows

    def vertex_rep_tuple(
        self,
        config: tuple[int, ...],
        reps: tuple[str, ...],
        vertex: int,
    ) -> tuple[str, ...]:
        out = [matter_rep(int(self.charges[config[vertex]]))]
        for edge, side in self.vertex_edges[vertex]:
            rep = reps[edge]
            out.append(rep if side == "tail" else dual_rep(rep))
        return tuple(out)

    def build_basis(self):
        basis = []
        flux_cache: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
        for config, charge_tuple, _ in state_energy_cutoff_configs(
            self.energies, self.charges, self.n_vertices, self.energy_cutoff
        ):
            if charge_tuple not in flux_cache:
                flux_cache[charge_tuple] = STRIP.flux_solutions(charge_tuple, self.incidence)
            for fluxes in flux_cache[charge_tuple]:
                rep_choices = [rep_options_for_flux(flux, self.rep_set) for flux in fluxes]
                for reps in product(*rep_choices):
                    vertex_bases = [
                        invariant_basis(self.vertex_rep_tuple(config, tuple(reps), vertex))
                        for vertex in range(self.n_vertices)
                    ]
                    if any(len(rows) == 0 for rows in vertex_bases):
                        continue
                    for intertwiners in product(*[range(len(rows)) for rows in vertex_bases]):
                        basis.append((tuple(config), tuple(reps), tuple(int(i) for i in intertwiners)))
        return basis, {item: pos for pos, item in enumerate(basis)}

    def vertex_tensor(self, state, vertex: int) -> np.ndarray:
        config, reps, intertwiners = state
        rep_tuple = self.vertex_rep_tuple(config, reps, vertex)
        return invariant_basis(rep_tuple)[intertwiners[vertex]]

    def vertex_open_matrix(
        self,
        vertex: int,
        active_edge: int,
        target_state,
        source_state,
        matter_map: np.ndarray,
    ) -> np.ndarray:
        target_tensor = self.vertex_tensor(target_state, vertex)
        source_tensor = self.vertex_tensor(source_state, vertex)
        incident = self.vertex_edges[vertex]
        active_pos = [idx for idx, (edge, _) in enumerate(incident) if edge == active_edge][0]
        active_axis = active_pos + 1

        target_order = [0, active_axis] + [axis for axis in range(1, target_tensor.ndim) if axis != active_axis]
        source_order = [0, active_axis] + [axis for axis in range(1, source_tensor.ndim) if axis != active_axis]
        tt = np.transpose(target_tensor, target_order)
        ss = np.transpose(source_tensor, source_order)
        other_shape_t = tt.shape[2:]
        other_shape_s = ss.shape[2:]
        if other_shape_t != other_shape_s:
            return np.zeros((tt.shape[1], ss.shape[1]), complex)
        tt = tt.reshape(tt.shape[0], tt.shape[1], int(np.prod(other_shape_t, dtype=int)))
        ss = ss.reshape(ss.shape[0], ss.shape[1], int(np.prod(other_shape_s, dtype=int)))
        return np.einsum("mao,mn,nbo->ab", tt.conj(), matter_map, ss, optimize=True)

    @lru_cache(maxsize=None)
    def edge_operator_element(
        self,
        edge: int,
        charge: int,
        target_state,
        source_state,
    ) -> complex:
        source_config, source_reps, _ = source_state
        target_config, target_reps, _ = target_state
        tail, head = self.edges[edge]
        source_tail_rep = matter_rep(int(self.charges[source_config[tail]]))
        target_tail_rep = matter_rep(int(self.charges[target_config[tail]]))
        source_head_rep = matter_rep(int(self.charges[source_config[head]]))
        target_head_rep = matter_rep(int(self.charges[target_config[head]]))

        delta_tail = (self.charges[target_config[tail]] - self.charges[source_config[tail]]) % 3
        delta_head = (self.charges[target_config[head]] - self.charges[source_config[head]]) % 3
        if (delta_tail + delta_head) % 3 != 0:
            return 0.0j
        if delta_tail == 0 or delta_head == 0:
            return 0.0j

        link_delta = (rep_triality(target_reps[edge]) - rep_triality(source_reps[edge])) % 3
        if link_delta == 0:
            return 0.0j
        kind = link_kind(link_delta)

        total = 0.0j
        factor_tail = charge_factor(int(delta_tail))
        factor_head = charge_factor(int(delta_head))
        for color_tail in range(3):
            tail_map = matter_component(factor_tail, source_tail_rep, target_tail_rep, color_tail)
            if not np.any(np.abs(tail_map) > 1e-14):
                continue
            tail_open = self.vertex_open_matrix(tail, edge, target_state, source_state, tail_map)
            if not np.any(np.abs(tail_open) > 1e-14):
                continue
            for color_head in range(3):
                head_map = matter_component(factor_head, source_head_rep, target_head_rep, color_head)
                if not np.any(np.abs(head_map) > 1e-14):
                    continue
                head_open = self.vertex_open_matrix(head, edge, target_state, source_state, head_map)
                if not np.any(np.abs(head_open) > 1e-14):
                    continue
                link = FULLCG.link_map(
                    kind,
                    color_tail,
                    color_head,
                    target_reps[edge],
                    source_reps[edge],
                )
                if not np.any(np.abs(link) > 1e-14):
                    continue
                total += np.einsum("ac,bd,abcd->", tail_open, head_open, link, optimize=True)
        return total

    def add_transition(self, rows, row, coeff):
        if row is not None and abs(coeff) > 1e-14:
            rows.append((row, coeff))

    def build_transitions(self):
        transitions = []
        for col, source_state in enumerate(self.basis):
            if self.progress_every and col and col % self.progress_every == 0:
                print(f"  transitions: {col}/{len(self.basis)}", flush=True)
            config, reps, intertwiners = source_state
            col_transitions = []
            config_list = list(config)
            for edge, (tail, head) in enumerate(self.edges):
                old_tail = config_list[tail]
                old_head = config_list[head]
                for mode in range(len(self.ann)):
                    for charge in (1, 2):
                        # create at tail, annihilate at head.
                        for new_tail, amp_tail in self.create[mode][charge][old_tail]:
                            for new_head, amp_head in self.ann[mode][charge][old_head]:
                                new_config = list(config_list)
                                new_config[tail] = new_tail
                                new_config[head] = new_head
                                self.emit_hop_targets(
                                    col_transitions,
                                    edge,
                                    charge,
                                    charge,
                                    tuple(new_config),
                                    reps,
                                    intertwiners,
                                    source_state,
                                    -self.hopping * amp_tail * amp_head,
                                )

                        # annihilate at tail, create at head.
                        for new_tail, amp_tail in self.ann[mode][charge][old_tail]:
                            for new_head, amp_head in self.create[mode][charge][old_head]:
                                new_config = list(config_list)
                                new_config[tail] = new_tail
                                new_config[head] = new_head
                                self.emit_hop_targets(
                                    col_transitions,
                                    edge,
                                    charge,
                                    -charge,
                                    tuple(new_config),
                                    reps,
                                    intertwiners,
                                    source_state,
                                    -self.hopping * amp_tail * amp_head,
                                )
            transitions.append(col_transitions)
        return transitions

    def emit_hop_targets(
        self,
        col_transitions,
        edge: int,
        charge: int,
        link_delta: int,
        new_config: tuple[int, ...],
        source_reps: tuple[str, ...],
        source_intertwiners: tuple[int, ...],
        source_state,
        matter_amp: complex,
    ) -> None:
        if abs(matter_amp) < 1e-14:
            return
        link_source_rep = source_reps[edge]
        for target_rep in (REPS_MINIMAL if self.rep_set == "minimal" else REPS_EXTENDED):
            if rep_triality(target_rep) != (rep_triality(link_source_rep) + link_delta) % 3:
                continue
            target_reps = list(source_reps)
            target_reps[edge] = target_rep
            target_reps = tuple(target_reps)
            changed_vertices = set(self.edges[edge])
            vertex_ranges = []
            for vertex in range(self.n_vertices):
                if vertex in changed_vertices:
                    rows = invariant_basis(self.vertex_rep_tuple(new_config, target_reps, vertex))
                    if not rows:
                        return
                    vertex_ranges.append(range(len(rows)))
                else:
                    vertex_ranges.append((source_intertwiners[vertex],))
            for target_intertwiners in product(*vertex_ranges):
                target_state = (new_config, target_reps, tuple(int(i) for i in target_intertwiners))
                row = self.index.get(target_state)
                if row is None:
                    continue
                coeff = matter_amp * self.edge_operator_element(edge, charge, target_state, source_state)
                self.add_transition(col_transitions, row, coeff)

    def matvec(self, vector: np.ndarray) -> np.ndarray:
        out = self.diagonal * vector
        for col, rows in enumerate(self.transitions):
            value = vector[col]
            if abs(value) < 1e-15:
                continue
            for row, coeff in rows:
                out[row] += coeff * value
        return out

    def matter_expectation(self, vector: np.ndarray) -> float:
        return float(np.real(np.vdot(vector, self.matter_diagonal * vector)))


def dense_lowest_if_small(model: SpinNetworkCell):
    if model.dim > 1000:
        return None
    eye = np.eye(model.dim, dtype=complex)
    dense = np.column_stack([model.matvec(eye[:, col]) for col in range(model.dim)])
    dense = (dense + dense.conj().T) / 2
    herm = float(np.linalg.norm(dense - dense.conj().T))
    eigvals, eigvecs = np.linalg.eigh(dense)
    return eigvals[:8], eigvecs[:, :8], herm


def lanczos_lowest(model: SpinNetworkCell, n_iter: int | None = None, n_eigs: int = 8, seed: int = 911):
    return PLAQUETTE.lanczos_lowest(model, n_iter=n_iter or model.n_lanczos, n_eigs=n_eigs, seed=seed)


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


def matter_dominated_gap(model: SpinNetworkCell, eigvals: np.ndarray, vectors: np.ndarray):
    ground = model.matter_expectation(vectors[:, 0])
    for idx in range(1, vectors.shape[1]):
        matter = model.matter_expectation(vectors[:, idx])
        if matter - ground > 0.5:
            return float(eigvals[idx] - eigvals[0]), float(matter - ground)
    return float("inf"), 0.0


def run_case(case: dict[str, Any]) -> dict[str, Any]:
    model = SpinNetworkCell(**case)
    exact = dense_lowest_if_small(model)
    if exact is None:
        eigvals, eigvecs = lanczos_lowest(model)
        herm = None
    else:
        eigvals, eigvecs, herm = exact
    full_gap, degeneracy = unique_gap(eigvals)
    matter_gap, matter_delta = matter_dominated_gap(model, eigvals, eigvecs)
    avg_transitions = float(np.mean([len(row) for row in model.transitions])) if model.transitions else 0.0
    max_transitions = int(max([len(row) for row in model.transitions], default=0))
    max_vertex_mult = 0
    for state in model.basis[: min(model.dim, 1000)]:
        config, reps, _ = state
        max_vertex_mult = max(
            max_vertex_mult,
            max(len(invariant_basis(model.vertex_rep_tuple(config, reps, v))) for v in range(model.n_vertices)),
        )
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
        "dense_hermiticity": herm,
        "sampled_max_vertex_multiplicity": max_vertex_mult,
    }


def build_verdict() -> dict[str, Any]:
    cases = [
        {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 0.0,
            "rep_set": "minimal",
            "beta": 0.5,
            "hopping": 0.2,
        },
        {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 2.0,
            "rep_set": "minimal",
            "beta": 0.5,
            "hopping": 0.2,
        },
        {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "minimal",
            "beta": 0.5,
            "hopping": 0.2,
            "n_lanczos": 50,
        },
        {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "minimal",
            "beta": 0.5,
            "hopping": 4.0,
            "n_lanczos": 50,
        },
        {
            "nx": 2,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 2.0,
            "rep_set": "minimal",
            "beta": 0.5,
            "hopping": 0.2,
            "n_lanczos": 50,
        },
    ]
    rows = [run_case(case) for case in cases]
    assert min(row["full_gap"] for row in rows) > 0.01
    assert any(row["avg_transitions_per_basis_state"] > 0 for row in rows)
    return {
        "verdict": "nonreduced_spin_network_scaling_constructed_gap_survives_finite_scan",
        "scope": (
            "Finite non-reduced SU(3) spin-network audit with explicit link endpoint "
            "indices in operator contractions and explicit vertex intertwiner labels. "
            "Local mirror states use the minimal charge lift q=0,1,2 -> 1,3,3bar."
        ),
        "omitted": [
            "neutral q=0 hopping remains omitted so the audit isolates charged SU(3) endpoint hopping",
            "full extended-rep and two-plaquette cutoff-4 scaling are still deep-box jobs for this non-reduced basis",
            "the local mirror-Fock block is charge-resolved, not a full microscopic SU(3) multiplet derivation",
            "magnetic plaquette transitions are not included in this hopping-focused scan",
        ],
        "rows": rows,
        "minimum_full_gap": min(row["full_gap"] for row in rows),
        "remaining_gate": (
            "add full plaquette Wilson moves to the same non-reduced multi-plaquette basis "
            "and push the extended-rep/two-plaquette cutoff scan on the deep box"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    verdict = build_verdict()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")

    print("TCH non-reduced spin-network mirror scaling audit")
    print(f"  verdict: {verdict['verdict']}")
    print(f"  wrote: {args.out.relative_to(ROOT)}")
    for row in verdict["rows"]:
        print(
            "    "
            f"{row['nx']}x{row['ny']} "
            f"spc={row['states_per_charge']} "
            f"ecut={row['energy_cutoff']} "
            f"reps={row['rep_set']} "
            f"t={row['hopping']} "
            f"dim={row['dimension']} "
            f"avgT={row['avg_transitions_per_basis_state']:.3g} "
            f"gap={row['full_gap']:.6g} "
            f"matter_gap={row['matter_dominated_gap']:.6g}"
        )
    print(f"  minimum full gap: {verdict['minimum_full_gap']:.6g}")
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
