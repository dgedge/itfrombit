#!/usr/bin/env python3
"""
Charge-block/Krylov audit for the two-plaquette mirror-Fock/SMG strip.

This ports the two-site charge-block representation to the 2D strip:

  * 2x1 plaquette strip: 6 sites, 7 links, 2 plaquette loops,
  * exact Z3 Gauss law solved by two independent loop fluxes,
  * no Python dictionary over individual local-state configurations in matvec,
  * local tensors are grouped by the six site charges,
  * hopping acts by small charge-sector matrices along tensor axes,
  * Lanczos is run directly against the block matvec.

The acceptance criterion is deliberately stricter than "nonzero gap":
the gap must visibly respond when t is made O(1). If it stays flat, the run is
declared truncation-suppressed rather than evidence for mirror-gap stability.
"""

import importlib.util
from itertools import product
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPATIAL_PATH = HERE / "css_spatial_mirror_fock_gauge_gap.py"
STRIP_PATH = HERE / "css_2d_strip_mirror_fock_lanczos.py"

N_SITES = 6
EDGES = [
    (0, 1),  # lower left horizontal
    (1, 2),  # lower right horizontal
    (3, 4),  # upper left horizontal
    (4, 5),  # upper right horizontal
    (0, 3),  # left vertical
    (1, 4),  # middle vertical
    (2, 5),  # right vertical
]
LEFT_LOOP = 0
RIGHT_LOOP = 1


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


spatial = load_module("spatial", SPATIAL_PATH)
old_strip = load_module("old_strip", STRIP_PATH)


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def fluxes_from_charges(charges, left_loop, right_loop):
    q0, q1, q2, q3, _, q5 = charges
    f4 = left_loop
    f6 = right_loop
    f0 = (-q0 - f4) % 3
    f1 = (q2 + f6) % 3
    f5 = (-q1 + f0 - f1) % 3
    f2 = (f4 - q3) % 3
    f3 = (q5 - f6) % 3
    return (f0, f1, f2, f3, f4, f5, f6)


def shift_link(fluxes, link, amount):
    out = list(fluxes)
    out[link] = (out[link] + amount) % 3
    return tuple(out)


def apply_axis(tensor, operator, axis):
    moved = np.moveaxis(tensor, axis, 0)
    transformed = np.tensordot(operator, moved, axes=(1, 0))
    return np.moveaxis(transformed, 0, axis)


def apply_two_axes(tensor, first_axis, first_op, second_axis, second_op):
    out = apply_axis(tensor, first_op, first_axis)
    return apply_axis(out, second_op, second_axis)


class ChargeBlockStripHamiltonian:
    def __init__(self, states_per_charge, beta, hopping, build_matter=True):
        self.states_per_charge = states_per_charge
        self.beta = beta
        self.hopping = hopping
        self.g2 = 1 / beta
        self.energies, self.charges, self.annihilators, _ = (
            spatial.local_mirror_dressed_block(states_per_charge)
        )
        self.sectors = [np.where(self.charges == charge)[0] for charge in range(3)]
        self.sector_pos = {}
        for charge, indices in enumerate(self.sectors):
            for pos, index in enumerate(indices):
                self.sector_pos[int(index)] = (charge, pos)
        self.sector_energies = [self.energies[indices] for indices in self.sectors]
        self.ann = self.restrict_annihilators()
        self.shape = (states_per_charge,) * N_SITES
        self.block_size = states_per_charge**N_SITES
        self.block_keys = []
        for charges in product(range(3), repeat=N_SITES):
            if sum(charges) % 3 != 0:
                continue
            for left_loop in range(3):
                for right_loop in range(3):
                    self.block_keys.append((tuple(charges), left_loop, right_loop))
        self.block_index = {key: pos for pos, key in enumerate(self.block_keys)}
        self.dim = len(self.block_keys) * self.block_size
        self.diagonal = self.build_diagonal(include_electric=True)
        self.matter_diagonal = (
            self.build_diagonal(include_electric=False) if build_matter else None
        )
        self.transitions = self.build_charge_transitions()

    def restrict_annihilators(self):
        restricted = {}
        for mode, components in enumerate(self.annihilators):
            for charge, matrix in enumerate(components):
                for target_charge in range(3):
                    rows = self.sectors[target_charge]
                    for source_charge in range(3):
                        cols = self.sectors[source_charge]
                        block = matrix[np.ix_(rows, cols)]
                        if np.linalg.norm(block) > 1e-12:
                            restricted[(mode, charge, target_charge, source_charge)] = block
        return restricted

    def block_slice(self, block_id):
        start = block_id * self.block_size
        return slice(start, start + self.block_size)

    def build_diagonal(self, include_electric):
        diagonal = np.empty(self.dim, dtype=float)
        site_mesh = np.indices(self.shape, sparse=True)
        charge_energy_cache = {}
        for charges in {key[0] for key in self.block_keys}:
            values = np.zeros(self.shape, dtype=float)
            for axis, charge in enumerate(charges):
                values = values + self.sector_energies[charge][site_mesh[axis]]
            charge_energy_cache[charges] = values.ravel()

        for block_id, (charges, left_loop, right_loop) in enumerate(self.block_keys):
            values = charge_energy_cache[charges].copy()
            if include_electric:
                fluxes = fluxes_from_charges(charges, left_loop, right_loop)
                electric = self.g2 * sum(spatial.electric_energy(flux) for flux in fluxes)
                values += electric
            diagonal[self.block_slice(block_id)] = values
        return diagonal

    def target_key_after_link_shift(self, charges, left_loop, right_loop, link, shift):
        fluxes = fluxes_from_charges(charges, left_loop, right_loop)
        shifted = shift_link(fluxes, link, shift)
        return shifted[4], shifted[6]

    def build_charge_transitions(self):
        transitions = [[] for _ in self.block_keys]
        for block_id, (charges, left_loop, right_loop) in enumerate(self.block_keys):
            # Magnetic plaquette terms. The independent loops are f4 and f6.
            for amount in [1, -1]:
                target = self.block_index[(charges, (left_loop - amount) % 3, right_loop)]
                transitions[block_id].append((target, None, None, None, None, -(self.beta / 2)))
                target = self.block_index[(charges, left_loop, (right_loop + amount) % 3)]
                transitions[block_id].append((target, None, None, None, None, -(self.beta / 2)))

            if abs(self.hopping) < 1e-15:
                continue

            for link, (tail, head) in enumerate(EDGES):
                tail_charge = charges[tail]
                head_charge = charges[head]
                for mode in range(6):
                    for charge in range(3):
                        # c_tail^dag c_head, accompanied by U_link^charge.
                        new_tail_charge = (tail_charge - charge) % 3
                        new_head_charge = (head_charge + charge) % 3
                        tail_ann = self.ann.get((mode, charge, tail_charge, new_tail_charge))
                        head_ann = self.ann.get((mode, charge, new_head_charge, head_charge))
                        if tail_ann is not None and head_ann is not None:
                            new_charges = list(charges)
                            new_charges[tail] = new_tail_charge
                            new_charges[head] = new_head_charge
                            new_left, new_right = self.target_key_after_link_shift(
                                tuple(new_charges), left_loop, right_loop, link, charge
                            )
                            target = self.block_index[(tuple(new_charges), new_left, new_right)]
                            transitions[block_id].append(
                                (
                                    target,
                                    tail,
                                    tail_ann.conj().T,
                                    head,
                                    head_ann,
                                    -self.hopping,
                                )
                            )

                        # c_tail c_head^dag, accompanied by U_link^-charge.
                        new_tail_charge = (tail_charge + charge) % 3
                        new_head_charge = (head_charge - charge) % 3
                        tail_ann = self.ann.get((mode, charge, new_tail_charge, tail_charge))
                        head_ann = self.ann.get((mode, charge, head_charge, new_head_charge))
                        if tail_ann is not None and head_ann is not None:
                            new_charges = list(charges)
                            new_charges[tail] = new_tail_charge
                            new_charges[head] = new_head_charge
                            new_left, new_right = self.target_key_after_link_shift(
                                tuple(new_charges), left_loop, right_loop, link, -charge
                            )
                            target = self.block_index[(tuple(new_charges), new_left, new_right)]
                            transitions[block_id].append(
                                (
                                    target,
                                    tail,
                                    tail_ann,
                                    head,
                                    head_ann.conj().T,
                                    -self.hopping,
                                )
                            )
        return transitions

    def matvec(self, vector):
        out = self.diagonal * vector
        for block_id, block_transitions in enumerate(self.transitions):
            source = vector[self.block_slice(block_id)].reshape(self.shape)
            if not np.any(np.abs(source) > 1e-15):
                continue
            for target, axis_a, op_a, axis_b, op_b, coeff in block_transitions:
                target_slice = self.block_slice(target)
                if op_a is None:
                    out[target_slice] += coeff * source.ravel()
                else:
                    moved = apply_two_axes(source, axis_a, op_a, axis_b, op_b)
                    out[target_slice] += coeff * moved.ravel()
        return out

    def matter_expectation(self, vector):
        if self.matter_diagonal is None:
            raise RuntimeError("matter_diagonal was not built for this model")
        return float(np.real(np.vdot(vector, self.matter_diagonal * vector)))


def lanczos_lowest(model, n_iter=42, n_eigs=12, seed=31415):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=model.dim) + 1j * rng.normal(size=model.dim)
    q = q / np.linalg.norm(q)
    q_prev = np.zeros_like(q)
    beta_prev = 0.0
    q_vectors = []
    alphas = []
    betas = []

    for _ in range(min(n_iter, model.dim)):
        q_vectors.append(q)
        v = model.matvec(q)
        alpha = float(np.real(np.vdot(q, v)))
        v = v - alpha * q - beta_prev * q_prev
        for basis_vector in q_vectors:
            v -= np.vdot(basis_vector, v) * basis_vector
        beta = float(np.linalg.norm(v))
        alphas.append(alpha)
        if beta < 1e-12:
            break
        betas.append(beta)
        q_prev = q
        q = v / beta
        beta_prev = beta

    tridiag = np.diag(alphas)
    for index in range(len(alphas) - 1):
        tridiag[index, index + 1] = betas[index]
        tridiag[index + 1, index] = betas[index]
    evals, evecs = np.linalg.eigh(tridiag)
    q_matrix = np.column_stack(q_vectors)
    vectors = q_matrix @ evecs[:, :n_eigs]
    return evals[:n_eigs], vectors


def format_gap(value):
    if np.isfinite(value):
        return f"{value:<.6g}"
    return "not in low window"


def low_window_matter_gap(model, evals, vectors):
    ground_matter = model.matter_expectation(vectors[:, 0])
    for index in range(1, vectors.shape[1]):
        delta = model.matter_expectation(vectors[:, index]) - ground_matter
        if delta > 0.5:
            return float(evals[index] - evals[0])
    return float("inf")


def old_to_block_vector(old_model, block_model, old_vector):
    block_vector = np.zeros(block_model.dim, dtype=complex)
    for old_pos, (config, fluxes) in enumerate(old_model.basis):
        charges = []
        local_positions = []
        for state in config:
            charge, pos = block_model.sector_pos[int(state)]
            charges.append(charge)
            local_positions.append(pos)
        key = (tuple(charges), fluxes[4], fluxes[6])
        block_id = block_model.block_index[key]
        local_flat = np.ravel_multi_index(tuple(local_positions), block_model.shape)
        block_vector[block_model.block_slice(block_id)][local_flat] = old_vector[old_pos]
    return block_vector


def agreement_check():
    hd("A. Matvec Agreement With Enumerated Strip")
    old_model = old_strip.StripHamiltonian(beta=0.5, hopping=0.2)
    block_model = ChargeBlockStripHamiltonian(states_per_charge=2, beta=0.5, hopping=0.2)
    assert old_model.dim == block_model.dim
    rng = np.random.default_rng(11)
    old_vector = rng.normal(size=old_model.dim) + 1j * rng.normal(size=old_model.dim)
    old_vector = old_vector / np.linalg.norm(old_vector)
    block_vector = old_to_block_vector(old_model, block_model, old_vector)
    old_result = old_to_block_vector(old_model, block_model, old_model.matvec(old_vector))
    block_result = block_model.matvec(block_vector)
    residual = float(np.linalg.norm(old_result - block_result) / np.linalg.norm(old_result))
    print(f"  dim={block_model.dim}, relative matvec residual={residual:.3e}")
    assert residual < 1e-10


def hermiticity_probe(model, seed=23):
    rng = np.random.default_rng(seed)
    left = rng.normal(size=model.dim) + 1j * rng.normal(size=model.dim)
    right = rng.normal(size=model.dim) + 1j * rng.normal(size=model.dim)
    mismatch = np.vdot(left, model.matvec(right)) - np.vdot(model.matvec(left), right)
    return float(abs(mismatch) / (np.linalg.norm(left) * np.linalg.norm(right)))


def gap_row(states_per_charge, beta, hopping, n_iter):
    model = ChargeBlockStripHamiltonian(
        states_per_charge=states_per_charge,
        beta=beta,
        hopping=hopping,
    )
    herm = hermiticity_probe(model)
    assert herm < 1e-10
    evals, vectors = lanczos_lowest(model, n_iter=n_iter)
    full_gap = float(evals[1] - evals[0])
    matter_gap = low_window_matter_gap(model, evals, vectors)
    return model.dim, full_gap, matter_gap, herm


def responsiveness_scan():
    hd("B. Two-Plaquette t-Responsiveness Scan")
    print("  beta=0.5, two plaquettes, exact Z3 Gauss law")
    print("  n/q   t      dim       full gap   matter gap        herm")
    rows = []
    cases = [
        (2, 40, [0.2, 1.0, 4.0]),
        (3, 24, [0.2, 4.0]),
    ]
    for states_per_charge, n_iter, hoppings in cases:
        for hopping in hoppings:
            dim, gap, m_gap, herm = gap_row(
                states_per_charge=states_per_charge,
                beta=0.5,
                hopping=hopping,
                n_iter=n_iter,
            )
            rows.append((states_per_charge, hopping, dim, gap, m_gap))
            print(
                f"  {states_per_charge:<5d} "
                f"{hopping:<6g} "
                f"{dim:<9d} "
                f"{gap:<10.6g} "
                f"{format_gap(m_gap):<17s} "
                f"{herm:.1e}"
            )

    by_cutoff = {}
    for states_per_charge, hopping, _, gap, _ in rows:
        by_cutoff.setdefault(states_per_charge, {})[hopping] = gap
    for states_per_charge, gaps in by_cutoff.items():
        signed_drop = gaps[0.2] - gaps[4.0]
        response = abs(signed_drop)
        relative = response / max(abs(gaps[0.2]), 1e-12)
        print(
            f"  response n/q={states_per_charge}: "
            f"signed_drop={signed_drop:.6g}, |response|={response:.6g}, "
            f"relative={relative:.3f}"
        )
    low_drop = abs(by_cutoff[2][0.2] - by_cutoff[2][4.0])
    high_drop = abs(by_cutoff[3][0.2] - by_cutoff[3][4.0])
    assert high_drop > 5 * low_drop
    assert high_drop > 0.05
    return rows


def main():
    hd("Scope")
    print(
        "This is a larger-cutoff two-plaquette strip audit. It tests whether "
        "the Krylov gap responds to hopping after porting the charge-block "
        "representation into the 2D strip."
    )
    agreement_check()
    rows = responsiveness_scan()
    hd("Verdict")
    print(
        "The charge-block strip reproduces the old enumerated n/q=2 matvec and "
        "then shows the required qualitative behavior: at larger local cutoff "
        "the two-plaquette gap visibly responds to O(1) hopping. A flat gap is "
        "therefore no longer accepted as evidence of stability."
    )
    print(f"\nminimum scanned full gap: {min(row[3] for row in rows):.6g}")
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
