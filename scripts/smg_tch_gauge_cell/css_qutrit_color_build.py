#!/usr/bin/env python3
"""
Qutrit-colour progression of the CSS/SU(3) compensator.

Purpose:
  css_matter_coupled_gauss_proxy.py showed that the Pauli/Z2 colour shadow
  recovers the CSS matter gap toward weak coupling but loses it at small beta
  / large g^2, the confining colour regime. ANCHOR §2.8 already flags the
  structural issue: the explicit dual mapping between the F2-closed colour
  algebra and the one-hot Z_{c_i} projector basis is open.

This script makes the next minimal qutrit-colour build:

  A. Local colour algebra:
     The two colour bits contain one lepton state |00> plus three colour states
     |01>, |10>, |11>. On the coloured subspace, define a qutrit clock

         Q = |r><r| + omega |g><g| + omega^2 |b><b|

     and reconstruct the one-hot projectors P_r, P_g, P_b by Fourier inversion.
     This supplies the missing local Z3 projector algebra, conditional on the
     coloured-subspace projector P_col.

  B. Why Pauli colour flips are not the colour operator:
     X_C0, X_C1, and X_C0 X_C1 are F2 translations. They exchange the lepton
     state with a colour state, so they do not preserve the qutrit colour
     subspace. The colour-cycle operator r -> g -> b -> r does preserve it and
     is the local Z3 orbit generator.

  C. Qutrit colour-orbit compensator for CSS X:
     Decompose each CSS X stabilizer under the adjoint Z3 colour-cycle action,

         X_g = sum_q X_g^(q),        C X_g^(q) C^dag = omega^q X_g^(q),

     and dress it with a qutrit clock compensator:

         Xhat_g = sum_q X_g^(q) x Z_comp^q.

     Xhat_g is invariant under C_matter x S_comp, Hermitian, and squares to one.

  D. Full local dressed CSS check:
     The Z-half must be dressed by the same qutrit orbit. Bare Z stabilizers do
     not commute with the qutrit-dressed X-half. Same-qutrit dressed Z and X
     do commute, leaving a threefold qutrit gauge/orbit degeneracy. The combined
     colour-cycle Gauss projector selects a unique physical ground state.

  E. Strong-coupling bookkeeping:
     Reprint the SU(3) Pauli-shadow matter gaps next to the qutrit-colour
     centre/orbit bookkeeping. The qutrit column is explicitly not a computed
     full colour matter gap; it is the local CSS stabilizer value after removing
     binary sign sectors. The full strong-coupling colour question is deferred
     to an open-link gauge Hamiltonian.

Scope:
  This is still not the full TCH SU(3) Wilson-link construction. It supplies a
  concrete qutrit colour algebra and compensator target for item 13 / §2.8, and
  it tells the next implementation what the colour-cell operator must reduce to.

numpy; self-asserting.
"""

import itertools

import numpy as np


G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
GENERATORS = [
    [1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 1, 1, 0, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 0, 1],
]
ALL_BITS = [tuple(bits) for bits in itertools.product([0, 1], repeat=8)]
OMEGA = np.exp(2j * np.pi / 3)

COLOR_LABELS = {
    (0, 0): "ell",
    (0, 1): "r",
    (1, 0): "g",
    (1, 1): "b",
}
COLOR_ORDER = ["ell", "r", "g", "b"]
COLOR_INDEX = {label: index for index, label in enumerate(COLOR_ORDER)}


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def idx(bits):
    return int("".join(map(str, bits)), 2)


def color_label(bits):
    return COLOR_LABELS[(bits[C0], bits[C1])]


def local_color_projectors():
    projectors = {}
    for label in COLOR_ORDER:
        matrix = np.zeros((4, 4), complex)
        matrix[COLOR_INDEX[label], COLOR_INDEX[label]] = 1
        projectors[label] = matrix
    return projectors


def lift_color_operator(local):
    matrix = np.zeros((256, 256), complex)
    for bits in ALL_BITS:
        source_label = color_label(bits)
        source = COLOR_INDEX[source_label]
        for target, coeff in enumerate(local[:, source]):
            if abs(coeff) < 1e-14:
                continue
            target_bits = list(bits)
            target_label = COLOR_ORDER[target]
            if target_label == "ell":
                target_bits[C0], target_bits[C1] = 0, 0
            elif target_label == "r":
                target_bits[C0], target_bits[C1] = 0, 1
            elif target_label == "g":
                target_bits[C0], target_bits[C1] = 1, 0
            elif target_label == "b":
                target_bits[C0], target_bits[C1] = 1, 1
            matrix[idx(tuple(target_bits)), idx(bits)] += coeff
    return matrix


def local_pauli_flip(bit):
    matrix = np.zeros((4, 4), complex)
    for c0, c1 in itertools.product([0, 1], repeat=2):
        source_label = COLOR_LABELS[(c0, c1)]
        target = [c0, c1]
        target[bit] ^= 1
        target_label = COLOR_LABELS[tuple(target)]
        matrix[COLOR_INDEX[target_label], COLOR_INDEX[source_label]] = 1
    return matrix


def local_color_cycle():
    matrix = np.zeros((4, 4), complex)
    cycle = {"ell": "ell", "r": "g", "g": "b", "b": "r"}
    for source, target in cycle.items():
        matrix[COLOR_INDEX[target], COLOR_INDEX[source]] = 1
    return matrix


def local_color_shift_qutrit():
    matrix = np.zeros((3, 3), complex)
    for state in range(3):
        matrix[(state + 1) % 3, state] = 1
    return matrix


def local_color_clock_partial(projectors):
    return (
        projectors["r"]
        + OMEGA * projectors["g"]
        + (OMEGA**2) * projectors["b"]
    )


def qutrit_ops():
    clock = np.diag([1, OMEGA, OMEGA**2]).astype(complex)
    shift = local_color_shift_qutrit()
    return clock, shift


def local_color_algebra_report():
    hd("A. F2 Colour Bits -> Qutrit One-Hot Projector Algebra")
    projectors = local_color_projectors()
    p_col = projectors["r"] + projectors["g"] + projectors["b"]
    q_clock = local_color_clock_partial(projectors)

    reconstructed = {
        "r": (p_col + q_clock + q_clock.conj().T) / 3,
        "g": (p_col + (OMEGA**2) * q_clock + OMEGA * q_clock.conj().T) / 3,
        "b": (p_col + OMEGA * q_clock + (OMEGA**2) * q_clock.conj().T) / 3,
    }

    max_reconstruction = 0.0
    for label in ["r", "g", "b"]:
        residual = float(np.linalg.norm(reconstructed[label] - projectors[label]))
        max_reconstruction = max(max_reconstruction, residual)
        print(f"  ||P_{label} - Fourier_qutrit(P_{label})|| = {residual:.3e}")

    x_c0 = local_pauli_flip(0)
    x_c1 = local_pauli_flip(1)
    x_both = x_c0 @ x_c1
    lepton = projectors["ell"]
    leakage = {
        "X_C0": float(np.linalg.norm(p_col @ x_c0 @ lepton) + np.linalg.norm(lepton @ x_c0 @ p_col)),
        "X_C1": float(np.linalg.norm(p_col @ x_c1 @ lepton) + np.linalg.norm(lepton @ x_c1 @ p_col)),
        "X_C0X_C1": float(np.linalg.norm(p_col @ x_both @ lepton) + np.linalg.norm(lepton @ x_both @ p_col)),
    }
    for name, value in leakage.items():
        print(f"  coloured/lepton leakage of {name:8s} = {value:.6g}")

    cycle = local_color_cycle()
    qutrit_shift = local_color_shift_qutrit()
    colored_block = cycle[1:, 1:]
    cycle_preserves_col = float(np.linalg.norm((np.eye(4) - p_col) @ cycle @ p_col))
    shift_residual = float(np.linalg.norm(colored_block - qutrit_shift))
    print(f"  colour cycle preserves coloured qutrit subspace = {cycle_preserves_col:.3e}")
    print(f"  coloured block equals qutrit shift              = {shift_residual:.3e}")

    assert max_reconstruction < 1e-12
    assert min(leakage.values()) > 1.0
    assert cycle_preserves_col < 1e-12
    assert shift_residual < 1e-12

    print(
        "  -> The one-hot colour projectors are a qutrit Fourier algebra on the "
        "coloured subspace. Pauli F2 flips are not colour rotations because they "
        "leak between lepton and coloured sectors; the Z3 colour cycle is the "
        "local qutrit orbit generator."
    )


def x_codeword(codeword):
    matrix = np.zeros((256, 256), complex)
    row = tuple(codeword)
    for bits in ALL_BITS:
        target = tuple(bit ^ flip for bit, flip in zip(bits, row))
        matrix[idx(target), idx(bits)] = 1
    return matrix


def z_codeword(codeword):
    z = np.diag([1, -1]).astype(complex)
    identity = np.eye(2, dtype=complex)
    matrix = np.array([[1]], complex)
    for value in codeword:
        matrix = np.kron(matrix, z if value else identity)
    return matrix


def color_cycle_unitary():
    return lift_color_operator(local_color_cycle())


def adjoint_components(operator, unitary):
    components = []
    powers = [np.linalg.matrix_power(unitary, n) for n in range(3)]
    inv_powers = [power.conj().T for power in powers]
    for charge in range(3):
        component = np.zeros_like(operator)
        for n in range(3):
            component += (OMEGA ** (-charge * n)) * powers[n] @ operator @ inv_powers[n]
        components.append(component / 3)
    return components


def qutrit_dress_components(components):
    clock, _ = qutrit_ops()
    out = np.zeros((256 * 3, 256 * 3), complex)
    for charge, component in enumerate(components):
        out += np.kron(component, np.linalg.matrix_power(clock, charge))
    return out


def color_orbit_compensator_report():
    hd("B. Qutrit Colour-Orbit Compensator for CSS X Stabilizers")
    cycle = color_cycle_unitary()
    _, compensator_shift = qutrit_ops()
    total_cycle = np.kron(cycle, compensator_shift)
    identity = np.eye(256 * 3, dtype=complex)

    max_reconstruction = 0.0
    max_covariance = 0.0
    max_dressed = 0.0
    max_hermiticity = 0.0
    max_involution = 0.0

    for index, codeword in enumerate(GENERATORS):
        xstab = x_codeword(codeword)
        components = adjoint_components(xstab, cycle)
        dressed = qutrit_dress_components(components)

        reconstruction = float(np.linalg.norm(sum(components) - xstab))
        bare_residual = float(np.linalg.norm(cycle @ xstab @ cycle.conj().T - xstab))
        dressed_residual = float(
            np.linalg.norm(total_cycle @ dressed @ total_cycle.conj().T - dressed)
        )
        hermiticity = float(np.linalg.norm(dressed - dressed.conj().T))
        involution = float(np.linalg.norm(dressed @ dressed - identity))
        charge_norms = [float(np.linalg.norm(component)) for component in components]

        for charge, component in enumerate(components):
            covariance = float(
                np.linalg.norm(
                    cycle @ component @ cycle.conj().T
                    - (OMEGA**charge) * component
                )
            )
            max_covariance = max(max_covariance, covariance)

        max_reconstruction = max(max_reconstruction, reconstruction)
        max_dressed = max(max_dressed, dressed_residual)
        max_hermiticity = max(max_hermiticity, hermiticity)
        max_involution = max(max_involution, involution)

        print(
            f"  Xstab{index}: bare orbit residual={bare_residual:.6g}; "
            f"component norms={np.round(charge_norms, 6)}; "
            f"reconstruct={reconstruction:.3e}; dressed={dressed_residual:.3e}; "
            f"Herm={hermiticity:.3e}; square={involution:.3e}"
        )

    assert max_reconstruction < 1e-10
    assert max_covariance < 1e-10
    assert max_dressed < 1e-10
    assert max_hermiticity < 1e-10
    assert max_involution < 1e-10

    print(
        "  -> The colour-orbit qutrit compensator is an explicit alternative to "
        "a scalar center-only repair: it treats the three colour states as a "
        "local qutrit orbit and makes the CSS X-half equivariant under that "
        "Z3 action."
    )


def full_dressed_css_report():
    hd("C. Full Same-Qutrit Dressed CSS Stabilizer Check")
    cycle = color_cycle_unitary()
    _, compensator_shift = qutrit_ops()
    total_cycle = np.kron(cycle, compensator_shift)
    identity = np.eye(256 * 3, dtype=complex)

    dressed_z = [
        qutrit_dress_components(adjoint_components(z_codeword(row), cycle))
        for row in GENERATORS
    ]
    dressed_x = [
        qutrit_dress_components(adjoint_components(x_codeword(row), cycle))
        for row in GENERATORS
    ]
    bare_z = [np.kron(z_codeword(row), np.eye(3, dtype=complex)) for row in GENERATORS]

    max_bare_z_vs_dressed_x = 0.0
    max_dressed_comm = 0.0
    max_dressed_invariance = 0.0
    all_dressed = dressed_z + dressed_x
    for z_stab in bare_z:
        for x_stab in dressed_x:
            max_bare_z_vs_dressed_x = max(
                max_bare_z_vs_dressed_x,
                float(np.linalg.norm(z_stab @ x_stab - x_stab @ z_stab)),
            )
    for left in all_dressed:
        max_dressed_invariance = max(
            max_dressed_invariance,
            float(np.linalg.norm(total_cycle @ left @ total_cycle.conj().T - left)),
        )
        for right in all_dressed:
            max_dressed_comm = max(
                max_dressed_comm,
                float(np.linalg.norm(left @ right - right @ left)),
            )

    projector = identity.copy()
    for stabilizer in all_dressed:
        projector = projector @ ((identity + stabilizer) / 2)
    projector_error = float(np.linalg.norm(projector @ projector - projector))
    ground_rank = round(float(np.real(np.trace(projector))))

    gauss_projector = (identity + total_cycle + total_cycle @ total_cycle) / 3
    gauss_error = float(np.linalg.norm(gauss_projector @ gauss_projector - gauss_projector))
    gauss_comm = float(np.linalg.norm(projector @ gauss_projector - gauss_projector @ projector))
    physical_ground_rank = round(float(np.real(np.trace(projector @ gauss_projector))))

    orbit_lock_strength = 2.0
    hamiltonian = -sum(all_dressed) - orbit_lock_strength * (
        total_cycle + total_cycle.conj().T
    ) / 2
    eigvals = np.linalg.eigvalsh(hamiltonian)
    unique = []
    for value in eigvals[:20]:
        if not unique or abs(value - unique[-1]) > 1e-8:
            unique.append(float(value))
    locked_gap = unique[1] - unique[0]
    locked_degeneracy = int(np.sum(np.isclose(eigvals, eigvals[0], atol=1e-8)))

    print(f"  max ||[bare Z, dressed X]||              = {max_bare_z_vs_dressed_x:.6g}")
    print(f"  max same-qutrit dressed stabilizer comm  = {max_dressed_comm:.3e}")
    print(f"  max same-qutrit dressed orbit residual   = {max_dressed_invariance:.3e}")
    print(f"  dressed-CSS projector rank before Gauss  = {ground_rank}")
    print(f"  ||P_CSS^2-P_CSS||                        = {projector_error:.3e}")
    print(f"  ||P_Gauss^2-P_Gauss||                    = {gauss_error:.3e}")
    print(f"  ||[P_CSS,P_Gauss]||                      = {gauss_comm:.3e}")
    print(f"  physical ground rank after Gauss         = {physical_ground_rank}")
    print(
        f"  with orbit-lock lambda={orbit_lock_strength:g}: "
        f"ground degeneracy={locked_degeneracy}, gap={locked_gap:.6g}"
    )

    assert max_bare_z_vs_dressed_x > 1.0
    assert max_dressed_comm < 1e-10
    assert max_dressed_invariance < 1e-10
    assert projector_error < 1e-10
    assert ground_rank == 3
    assert gauss_error < 1e-10
    assert gauss_comm < 1e-10
    assert physical_ground_rank == 1
    assert locked_degeneracy == 1
    assert locked_gap > 1.9

    print(
        "  -> Dressing only X is insufficient. Dressing Z and X by the same "
        "qutrit colour orbit preserves the CSS commuting algebra. The residual "
        "threefold qutrit degeneracy is gauge/orbit redundancy, removed by the "
        "combined colour-cycle Gauss projector or by an equivalent commuting "
        "orbit-lock term."
    )


def su3_reps(level):
    return [(p, q) for p in range(level + 1) for q in range(level + 1 - p)]


def su3_casimir(p, q):
    return (p * p + q * q + p * q + 3 * p + 3 * q) / 3


def su3_character_gap(level, beta):
    reps = su3_reps(level)
    index = {rep: i for i, rep in enumerate(reps)}
    dim = len(reps)
    mult_fund = np.zeros((dim, dim), complex)
    for col, (p, q) in enumerate(reps):
        for target in [(p + 1, q), (p - 1, q + 1), (p, q - 1)]:
            if target in index:
                mult_fund[index[target], col] += 1

    g2 = 1 / beta
    electric = np.diag([g2 * su3_casimir(*rep) for rep in reps]).astype(complex)
    magnetic = -(beta / 2) * (mult_fund + mult_fund.conj().T)
    eigvals = np.linalg.eigvalsh(electric + magnetic)[:2]
    return float(eigvals[1] - eigvals[0])


def pauli_shadow_no_flip_matter_gap(beta, level):
    # Same normalized Wilson-character source model as css_matter_coupled_gauss_proxy.py.
    def eig0(source):
        reps = su3_reps(level)
        index = {rep: i for i, rep in enumerate(reps)}
        dim = len(reps)
        mult_fund = np.zeros((dim, dim), complex)
        for col, (p, q) in enumerate(reps):
            for target in [(p + 1, q), (p - 1, q + 1), (p, q - 1)]:
                if target in index:
                    mult_fund[index[target], col] += 1
        wilson = (mult_fund + mult_fund.conj().T) / 6
        g2 = 1 / beta
        hamiltonian = np.diag([g2 * su3_casimir(*rep) for rep in reps]).astype(complex)
        hamiltonian -= (beta + source) * wilson
        return float(np.linalg.eigvalsh(hamiltonian)[0])

    ground = eig0(4)
    # closed flux P_X=+ and remove the global sign partner leaves sectors with sum x = 0.
    no_flip_x_gap = eig0(0) - ground
    return min(2.0, no_flip_x_gap)


def strong_coupling_bookkeeping_report():
    hd("D. Strong-Coupling Matter-Gap Bookkeeping")
    print(
        "  The Pauli-shadow column is the old Z2 no-flip SU(3) diagnostic. "
        "The qutrit-local column is not a computed full-colour gap: it is the "
        "local CSS stabilizer value after the Z3 colour/orbit bookkeeping has "
        "removed binary sign sectors."
    )
    print("  beta       trunc  Pauli no-flip matter   qutrit-local CSS value")

    rows = [(0.5, 18), (1, 18), (2, 20), (5, 24), (10, 28), (20, 32), (30, 36)]
    pauli_mins = []
    qutrit_mins = []
    for beta, level in rows:
        pauli_gap = pauli_shadow_no_flip_matter_gap(beta, level)
        qutrit_local_css_value = 2.0
        pauli_mins.append(pauli_gap)
        qutrit_mins.append(qutrit_local_css_value)
        print(
            f"  {beta:<10g} "
            f"{level:<6d} "
            f"{pauli_gap:<23.6g} "
            f"{qutrit_local_css_value:<15.6g}"
        )

    assert min(pauli_mins) < 1.0
    assert min(qutrit_mins) == 2.0
    print(
        "  -> This is bookkeeping, not a full strong-coupling proof. It shows "
        "what value the local dressed CSS algebra would have if the only problem "
        "were the binary sign sectors. The actual colour-dynamical gap must be "
        "computed in an open-link gauge Hamiltonian."
    )


def main():
    hd("Scope")
    print(
        "This progresses the qutrit-colour build from a center-only repair to a "
        "local colour-orbit algebra: one-hot projectors from a qutrit clock, "
        "Z3 colour-cycle compensation for CSS X, and a strong-coupling comparison "
        "against the old Pauli shadow."
    )
    local_color_algebra_report()
    color_orbit_compensator_report()
    full_dressed_css_report()
    strong_coupling_bookkeeping_report()
    hd("Verdict")
    print(
        "The natural progression is viable at the algebraic/proxy level: keep the "
        "F2 register as the lepton-plus-three-colour label, but add a qutrit "
        "operator algebra on the coloured subspace. The colour-cell operator must "
        "use the Z3 colour cycle/clock algebra, not Pauli colour-bit flips, and "
        "both CSS halves must be dressed by the same qutrit orbit charge."
    )
    print(
        "\nNext hard test: replace the bookkeeping in section D by an open-link "
        "SU(3) qutrit-colour Wilson/TCH Hamiltonian with vertex Gauss projection, "
        "then scan the confining beta<=1 regime for the actual mirror matter gap."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
