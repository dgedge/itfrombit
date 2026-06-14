#!/usr/bin/env python3
"""
Audit whether the spatial hopping survives the local SMG truncation.

The earlier finite-chain scan kept only a few local SMG eigenstates per Z3
charge sector. Its gap barely moved even when t was made O(1), which is a
warning that the projected fermion hopping was mostly removed by the truncation.

This script keeps the same two-site, open-link, Gauss-projected setup, but uses
a charge-block matrix-free representation. That lets the two-site problem be
run to much larger local cutoffs without forming a dense Hamiltonian. The
hopping convention matches css_spatial_mirror_fock_gauge_gap.py: one oriented
hopping term followed by Hermitian symmetrization, i.e. each explicit direction
has coefficient -t/2.
"""

import importlib.util
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPATIAL_PATH = HERE / "css_spatial_mirror_fock_gauge_gap.py"


def load_spatial_module():
    spec = importlib.util.spec_from_file_location("spatial", SPATIAL_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


spatial = load_spatial_module()


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def pack(blocks):
    return np.concatenate([block.ravel() for block in blocks])


def unpack(vector, counts):
    blocks = []
    offset = 0
    for charge in range(3):
        rows = counts[charge]
        cols = counts[(-charge) % 3]
        size = rows * cols
        blocks.append(vector[offset : offset + size].reshape((rows, cols)))
        offset += size
    return blocks


class TwoSiteChargeBlockHamiltonian:
    def __init__(self, states_per_charge, beta, hopping):
        self.states_per_charge = states_per_charge
        self.beta = beta
        self.hopping = hopping
        self.g2 = 1 / beta
        self.energies, self.charges, self.annihilators, self.charge_counts = (
            spatial.local_mirror_dressed_block(states_per_charge)
        )
        self.sectors = [np.where(self.charges == charge)[0] for charge in range(3)]
        self.counts = [len(indices) for indices in self.sectors]
        self.sector_energies = [self.energies[indices] for indices in self.sectors]
        self.ann = self.restrict_annihilators()
        self.dim = sum(self.counts[q] * self.counts[(-q) % 3] for q in range(3))
        self.matter_diagonal = self.build_matter_diagonal()

    def restrict_annihilators(self):
        restricted = []
        for mode_components in self.annihilators:
            by_charge = []
            for component in mode_components:
                target_blocks = []
                for target_charge in range(3):
                    rows = self.sectors[target_charge]
                    source_blocks = []
                    for source_charge in range(3):
                        cols = self.sectors[source_charge]
                        source_blocks.append(component[np.ix_(rows, cols)])
                    target_blocks.append(source_blocks)
                by_charge.append(target_blocks)
            restricted.append(by_charge)
        return restricted

    def build_matter_diagonal(self):
        diagonal = []
        for left_charge in range(3):
            right_charge = (-left_charge) % 3
            values = (
                self.sector_energies[left_charge][:, None]
                + self.sector_energies[right_charge][None, :]
            )
            diagonal.append(values.ravel())
        return np.concatenate(diagonal)

    def matvec(self, vector):
        source = unpack(vector, self.counts)
        target = []
        for left_charge in range(3):
            right_charge = (-left_charge) % 3
            electric = self.g2 * spatial.electric_energy((-left_charge) % 3)
            diagonal = (
                self.sector_energies[left_charge][:, None]
                + self.sector_energies[right_charge][None, :]
                + electric
            )
            target.append(diagonal * source[left_charge])

        prefactor = -0.5 * self.hopping
        if abs(prefactor) < 1e-15:
            return pack(target)

        for left_charge in range(3):
            right_charge = (-left_charge) % 3
            block = source[left_charge]
            for mode in range(len(self.ann)):
                for charge in range(3):
                    # c_L^dag c_R: q_L -> q_L-charge, q_R -> q_R+charge.
                    next_left = (left_charge - charge) % 3
                    next_right = (right_charge + charge) % 3
                    left_op = self.ann[mode][charge][left_charge][next_left].conj().T
                    right_op = self.ann[mode][charge][next_right][right_charge]
                    target[next_left] += prefactor * left_op @ block @ right_op.T

                    # c_L c_R^dag: q_L -> q_L+charge, q_R -> q_R-charge.
                    next_left = (left_charge + charge) % 3
                    next_right = (right_charge - charge) % 3
                    left_op = self.ann[mode][charge][next_left][left_charge]
                    right_op = self.ann[mode][charge][right_charge][next_right].conj().T
                    target[next_left] += prefactor * left_op @ block @ right_op.T

        return pack(target)

    def matter_expectation(self, vector):
        return float(np.real(np.vdot(vector, self.matter_diagonal * vector)))


def lanczos_lowest(model, n_iter=140, n_eigs=20, seed=2718):
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


def unique_gap(eigvals, tol=1e-8):
    ground = float(eigvals[0])
    for value in eigvals[1:]:
        if abs(value - ground) > tol:
            return float(value - ground)
    return 0.0


def dense_gap(model):
    if model.dim > 900:
        return None
    dense = np.column_stack([model.matvec(np.eye(model.dim, dtype=complex)[:, i]) for i in range(model.dim)])
    herm = float(np.linalg.norm(dense - dense.conj().T))
    assert herm < 1e-9
    return unique_gap(np.linalg.eigvalsh(dense))


def matter_gap(model, evals, vectors):
    ground_matter = model.matter_expectation(vectors[:, 0])
    for index in range(1, vectors.shape[1]):
        delta = model.matter_expectation(vectors[:, index]) - ground_matter
        if delta > 0.5:
            return float(evals[index] - evals[0])
    return float("inf")


def format_gap(value):
    if np.isfinite(value):
        return f"{value:<.6g}"
    return "not in low window"


def gap_at(states_per_charge, beta, hopping):
    model = TwoSiteChargeBlockHamiltonian(states_per_charge, beta, hopping)
    exact = dense_gap(model)
    evals, vectors = lanczos_lowest(model)
    gap = float(evals[1] - evals[0]) if exact is None else exact
    return model.dim, gap, matter_gap(model, evals, vectors), exact is not None


def agreement_check():
    hd("A. Agreement With Earlier Dense Chain")
    print("  beta=0.5, t=0.2, L=2")
    print("  n/q   old dense   block dense   abs diff")
    for states_per_charge in [2, 4, 6]:
        hamiltonian, _, _ = spatial.build_chain_hamiltonian(
            length=2,
            beta=0.5,
            hopping=0.2,
            states_per_charge=states_per_charge,
        )
        old_gap = unique_gap(np.linalg.eigvalsh(hamiltonian))
        _, block_gap, _, exact = gap_at(states_per_charge, beta=0.5, hopping=0.2)
        assert exact
        diff = abs(old_gap - block_gap)
        print(
            f"  {states_per_charge:<5d} "
            f"{old_gap:<11.6g} "
            f"{block_gap:<12.6g} "
            f"{diff:.2e}"
        )
        assert diff < 1e-8


def responsiveness_scan():
    hd("B. Hopping Responsiveness Versus Local Cutoff")
    print("  beta=0.5, two sites, open Z3 link, Gauss projected")
    print("  n/q   dim    gap(t=0.2)  gap(t=1)    gap(t=4)    matter gap(t=4)")
    rows = []
    for states_per_charge in [2, 4, 6, 8, 12, 16, 24, 32]:
        dim = None
        row = []
        matter_at_4 = None
        for hopping in [0.2, 1.0, 4.0]:
            dim, gap, m_gap, _ = gap_at(states_per_charge, beta=0.5, hopping=hopping)
            row.append(gap)
            if hopping == 4.0:
                matter_at_4 = m_gap
        rows.append((states_per_charge, dim, *row, matter_at_4))
        print(
            f"  {states_per_charge:<5d} "
            f"{dim:<6d} "
            f"{row[0]:<11.6g} "
            f"{row[1]:<11.6g} "
            f"{row[2]:<11.6g} "
            f"{format_gap(matter_at_4)}"
        )

    assert rows[0][4] > 1.9
    assert rows[-1][4] < 1.0
    return rows


def high_cutoff_tail():
    hd("C. High-Cutoff Tail At t=4")
    print("  beta=0.5, two sites, open Z3 link, Gauss projected")
    print("  n/q   dim     full gap   matter-tagged gap")
    rows = []
    for states_per_charge in [48, 64]:
        dim, gap, m_gap, _ = gap_at(states_per_charge, beta=0.5, hopping=4.0)
        rows.append((states_per_charge, dim, gap, m_gap))
        print(
            f"  {states_per_charge:<5d} "
            f"{dim:<7d} "
            f"{gap:<10.6g} "
            f"{format_gap(m_gap)}"
        )
    assert min(row[2] for row in rows) < 0.5
    return rows


def main():
    hd("Scope")
    print(
        "This is not a new stability claim. It is a truncation audit: does the "
        "projected hopping actually move the two-site mirror gap as the retained "
        "local charge-sector Hilbert space grows?"
    )
    agreement_check()
    rows = responsiveness_scan()
    tail = high_cutoff_tail()
    hd("Verdict")
    print(
        "The low states/charge truncation did suppress hopping. At n/q=2 the "
        "gap is almost inert even at t=4, but by n/q=32..64 the same hopping "
        "substantially lowers the gap. Therefore the earlier spatial-chain "
        "survival statement is inconclusive as a mirror-gap test."
    )
    print(
        "\nThe useful next target is the same charge-block/Krylov method at "
        "larger n/q and then on the two-plaquette strip. The acceptance test "
        "should require a visible t-response before interpreting any nonzero "
        "gap as stability."
    )
    min_gap = min([row[4] for row in rows] + [row[2] for row in tail])
    print(f"\nminimum audited full gap at t=4: {min_gap:.6g}")
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
