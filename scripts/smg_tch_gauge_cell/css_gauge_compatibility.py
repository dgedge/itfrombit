#!/usr/bin/env python3
"""
Gauge-compatibility audit for the cell-level CSS SMG construction.

Context:
  smg_construction.py proves an exact local result: the self-dual [8,4,4]
  CSS Hamiltonian H = -sum Z_stab - sum X_stab has a unique, gapped,
  entangled ground state with no low-weight bilinear order.

This script checks the next question that matters for the chiral lattice
gauge problem:

  Do the bare cell X-stabilizers commute with the encoded SM gauge actions?

Answer:
  No. As operators acting only on the matter-cell register, the X-stabilizers
  do not commute with the encoded U(1)_Y, SU(2)_L, or SU(3)_c actions used in
  r2_smg_operator.py. This is not by itself a contradiction of item 13: if the
  X-half is really the gauge-sector Wilson/plaquette half, it must be dressed
  with gauge-link degrees of freedom. But it means the local CSS gap is not yet
  the dynamically gauged non-abelian mirror interaction.

The useful sharpened target is therefore:

  Construct gauge-dressed CSS X operators whose matter-cell reduction contains
  the X-half used in smg_construction.py and whose full commutators with the
  gauged constraints vanish.

The dressing used here is the minimal algebraic compensator construction:

  1. Close each bare stabilizer under the adjoint action [G, .] of the encoded
     gauge checks.
  2. Add a link/plaquette compensator Hilbert space carrying the dual orbit
     representation.
  3. Form D = sum_a O_a x |a><vac|, where O_0 is the original stabilizer
     component and {O_a} is its gauge orbit.

Then [G_matter + G_link, D] = 0 by the dual-representation construction.
That exact zero is generic: it is not the physics test. The useful content is
that the adjoint orbit closes finitely with small numerical residual, and the
explicit block-sparse commutator below makes the cancellation auditable rather
than hiding it in a scalar rep - rep helper.

This remains a representation-level dressing, not a dynamical Wilson-action
construction and not a proof that the dressed interaction preserves the SMG gap.

numpy; self-asserting.
"""

import itertools
import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1, -1]).astype(complex)
P0 = np.array([[1, 0], [0, 0]], complex)
P1 = np.array([[0, 0], [0, 1]], complex)
SR = np.array([[0, 1], [0, 0]], complex)
SL = np.array([[0, 0], [1, 0]], complex)

G0, G1, LQ, C0, C1, I3, CHI, W = range(8)


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def op(singles):
    out = np.array([[1]], complex)
    for bit in range(8):
        out = np.kron(out, singles.get(bit, I2))
    return out


def x_codeword(codeword):
    return op({bit: X for bit, value in enumerate(codeword) if value})


def z_codeword(codeword):
    return op({bit: Z for bit, value in enumerate(codeword) if value})


def ctrl(control_bit, control_value, target_bit, target_op):
    active = P0 if control_value == 0 else P1
    inactive = P1 if control_value == 0 else P0
    return op({control_bit: active, target_bit: target_op}) + op(
        {control_bit: inactive, target_bit: I2}
    )


def commutator(a, b):
    return a @ b - b @ a


def norm(a):
    return float(np.linalg.norm(a))


def idx(bits):
    return int("".join(map(str, bits)), 2)


def charge_zf(bits):
    return 1 if bits[I3] == 0 else -1


def sum_colour_projectors(bits):
    return -3 if (bits[C0], bits[C1]) == (0, 0) else -1


def electric_charge(bits):
    return 0.5 * charge_zf(bits) + (1 / 3) * sum_colour_projectors(bits) + 0.5


def weak_t3(bits):
    if bits[CHI] != 0:
        return 0.0
    return 0.5 if bits[I3] == 0 else -0.5


def build_gauge_actions():
    all_bits = [tuple(bits) for bits in itertools.product([0, 1], repeat=8)]

    hypercharge = np.diag(
        [electric_charge(bits) - weak_t3(bits) for bits in all_bits]
    ).astype(complex)

    weak_t3_action = np.diag([weak_t3(bits) for bits in all_bits]).astype(complex)
    weak_raise = ctrl(CHI, 0, I3, SR)
    weak_lower = ctrl(CHI, 0, I3, SL)

    colour_cycle = np.zeros((256, 256), complex)
    cycle = {
        (0, 1): (1, 0),
        (1, 0): (1, 1),
        (1, 1): (0, 1),
        (0, 0): (0, 0),
    }
    for bits in all_bits:
        target = list(bits)
        target[C0], target[C1] = cycle[(bits[C0], bits[C1])]
        colour_cycle[idx(tuple(target)), idx(bits)] = 1

    return {
        "U(1)_Y": hypercharge,
        "SU(2)_L T3": weak_t3_action,
        "SU(2)_L T+": weak_raise,
        "SU(2)_L T-": weak_lower,
        "SU(3)_c cycle": colour_cycle,
        "SU(3)_c cycle^-1": colour_cycle.conj().T,
    }


def report_commutators(label, operators, gauge_actions):
    hd(label)
    for op_name, operator in operators.items():
        print(op_name)
        for gauge_name, gauge_op in gauge_actions.items():
            value = norm(commutator(operator, gauge_op))
            print(f"  ||[{op_name}, {gauge_name}]|| = {value:.6g}")


def orthonormal_orbit(seed, gauge_actions, tol=1e-9):
    """Return an orthonormal Hilbert-Schmidt basis for the adjoint gauge orbit."""
    basis_vecs = []
    basis_mats = []
    queue = []

    def add(matrix):
        vector = matrix.reshape(-1).astype(complex)
        # Modified Gram-Schmidt; repeat once because the orbit bases are moderately large.
        for _ in range(2):
            for basis in basis_vecs:
                vector -= np.vdot(basis, vector) * basis
        size = np.linalg.norm(vector)
        if size <= tol:
            return False
        vector = vector / size
        basis_vecs.append(vector)
        basis_mats.append(vector.reshape(seed.shape))
        queue.append(basis_mats[-1])
        return True

    add(seed)
    cursor = 0
    while cursor < len(queue):
        matrix = queue[cursor]
        cursor += 1
        for gauge_op in gauge_actions.values():
            add(commutator(gauge_op, matrix))

    return basis_vecs, basis_mats


def representation_matrix(basis_vecs, basis_mats, gauge_op):
    """Matrix R with [G, O_j] = sum_i R[i,j] O_i."""
    basis = np.vstack(basis_vecs)  # rows are orthonormal basis vectors
    commutator_vectors = np.column_stack(
        [commutator(gauge_op, matrix).reshape(-1) for matrix in basis_mats]
    )
    rep = basis.conj() @ commutator_vectors
    residuals = commutator_vectors - basis.T @ rep
    max_residual = float(np.max(np.linalg.norm(residuals, axis=0)))
    return rep, max_residual


def explicit_dressed_commutator_residual(basis_mats, rep, gauge_op):
    """
    Build the occupied blocks of [G_matter + G_link, D].

    Link basis convention:
      |0> is the vacuum/anchor state.
      |j+1> is the compensator state paired with matter-orbit operator O_j.

    D = sum_j O_j x |j+1><0|.
    G_link|0> = 0.
    G_link|j+1> = sum_i L[i,j]|i+1>, with L = -R^T.

    A dense matrix would have dimension 256*(orbit_dim+1), so this routine keeps
    the occupied link-column blocks, each a 256x256 matter matrix. This is still
    the generic dual-representation cancellation, but it is an explicit block
    commutator of constructed operators rather than norm(rep - rep).
    """
    link_rep = -rep.T
    dim = len(basis_mats)
    blocks = {
        row: np.zeros_like(basis_mats[0])
        for row in range(dim)
    }

    # Matter part: sum_j [G_matter, O_j] x |j+1><0|.
    for row, matrix in enumerate(basis_mats):
        blocks[row] += commutator(gauge_op, matrix)

    # Link part: sum_j O_j x [G_link, |j+1><0|].
    for source, matrix in enumerate(basis_mats):
        for target in range(dim):
            coeff = link_rep[target, source]
            if abs(coeff) > 1e-14:
                blocks[target] += coeff * matrix

    block_norms = [norm(block) for block in blocks.values()]
    full_norm = float(np.sqrt(sum(value * value for value in block_norms)))
    return max(block_norms), full_norm


def report_dressed_orbits(label, operators, gauge_actions):
    hd(label)
    max_full_residual = 0.0
    max_closure_residual = 0.0
    orbit_dims = []

    for op_name, operator in operators.items():
        basis_vecs, basis_mats = orthonormal_orbit(operator, gauge_actions)
        orbit_dims.append(len(basis_mats))
        print(f"{op_name}")
        print(
            f"  matter orbit dimension = {len(basis_mats)}; "
            f"link compensator dimension = {len(basis_mats) + 1}"
        )
        for gauge_name, gauge_op in gauge_actions.items():
            rep, closure_residual = representation_matrix(basis_vecs, basis_mats, gauge_op)
            max_block_residual, full_residual = explicit_dressed_commutator_residual(
                basis_mats, rep, gauge_op
            )
            max_closure_residual = max(max_closure_residual, closure_residual)
            max_full_residual = max(max_full_residual, full_residual)
            print(
                f"  dressed ||[D, {gauge_name}_matter + {gauge_name}_link]|| "
                f"(block-sparse) = {full_residual:.3e}; "
                f"max block = {max_block_residual:.3e}; "
                f"orbit closure residual = {closure_residual:.3e}"
            )

    print(
        f"\n  max orbit dimension = {max(orbit_dims)}; "
        f"max closure residual = {max_closure_residual:.3e}; "
        f"max explicit dressed commutator residual = {max_full_residual:.3e}"
    )
    assert max_closure_residual < 1e-7
    assert max_full_residual < 1e-10
    return max_full_residual


def main():
    # RM(1,3) generator rows: the self-dual [8,4,4] extended Hamming code.
    generators = [
        [1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 1, 1, 1, 1],
        [0, 0, 1, 1, 0, 0, 1, 1],
        [0, 1, 0, 1, 0, 1, 0, 1],
    ]

    x_stabs = [x_codeword(row) for row in generators]
    z_stabs = [z_codeword(row) for row in generators]
    h_x = -sum(x_stabs)
    h_z = -sum(z_stabs)
    h_css = h_x + h_z
    z_chi_z_w = op({CHI: Z, W: Z})

    gauge_actions = build_gauge_actions()

    report_commutators(
        "Bare matter-cell X stabilizers versus encoded gauge actions",
        {
            f"Xstab{i} support={[bit for bit, value in enumerate(row) if value]}": stab
            for i, (row, stab) in enumerate(zip(generators, x_stabs))
        },
        gauge_actions,
    )

    report_commutators(
        "Hamiltonian-level commutators",
        {
            "H_X = -sum Xstab": h_x,
            "H_Z = -sum Zstab": h_z,
            "H_CSS = H_X + H_Z": h_css,
            "Z_chi Z_W anti-mirror": z_chi_z_w,
        },
        gauge_actions,
    )

    dressed_x_residual = report_dressed_orbits(
        "Gauge-dressed X stabilizers: full matter-plus-link commutators",
        {
            f"Dressed Xstab{i} support={[bit for bit, value in enumerate(row) if value]}": stab
            for i, (row, stab) in enumerate(zip(generators, x_stabs))
        },
        gauge_actions,
    )

    hx_norms = [norm(commutator(h_x, gauge_op)) for gauge_op in gauge_actions.values()]
    hcss_norms = [norm(commutator(h_css, gauge_op)) for gauge_op in gauge_actions.values()]
    anti_mirror_norms = [
        norm(commutator(z_chi_z_w, gauge_op)) for gauge_op in gauge_actions.values()
    ]

    assert any(value > 1e-9 for value in hx_norms)
    assert any(value > 1e-9 for value in hcss_norms)
    assert all(value < 1e-9 for value in anti_mirror_norms)
    assert dressed_x_residual < 1e-10

    hd("Verdict")
    print(
        "The cell-level CSS X-half is an exact local SMG interaction, but as a bare "
        "matter-register operator it is not gauge invariant under the encoded SM "
        "actions. The anti-mirror diagnostic Z_chi Z_W is gauge invariant, but it is "
        "only a diagonal stabilizer, not the entangling X-half that produces the "
        "unique CSS gap."
    )
    print(
        "\nThe dual-orbit compensator construction above produces gauge-dressed "
        "X-stabilizer operators whose explicit block-sparse matter-plus-link "
        "commutators vanish to numerical precision. This exact zero is generic "
        "once the dual link representation is chosen; the nontrivial data are "
        "the finite orbit dimensions and closure residuals."
    )
    print(
        "\nRemaining gates: identify these compensator orbit/link spaces with actual "
        "TCH gauge-cell Wilson/plaquette degrees of freedom, decide how the Z-half "
        "is represented as Gauss-law constraints or is similarly dressed, make the "
        "dressed terms Hermitian in the full gauge algebra, and verify that the "
        "dressed interaction still has the CSS symmetric gap under background Wilson "
        "twists or dynamical gauge fluctuations."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
