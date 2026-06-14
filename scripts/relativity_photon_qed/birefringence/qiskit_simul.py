import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector


def build_circlette_lattice():
	"""Defines the 4.8.8 Archimedean coupling map for 12 qubits."""
	octagon = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0)]
	square1 = [(0, 8), (8, 9), (9, 1)]
	square2 = [(4, 10), (10, 11), (11, 5)]
	return octagon + square1 + square2


def build_evolution_step(edges, theta, num_qubits=12):
	"""
	Builds the unitary operator U for a single time step.
	By encapsulating this in a single circuit, we can use Qiskit's
	.inverse() safely, bypassing any manual gate-ordering bugs.
	"""
	qc = QuantumCircuit(num_qubits)
	for i in range(num_qubits):
		qc.rx(theta, i)
		qc.rz(theta, i)
	for edge in edges:
		qc.cx(edge[0], edge[1])
	return qc


def run_sanity_check(edges):
	"""Change 1: Sanity-check run (t=1, no perturbation)."""
	print("Running Sanity Check (t=1, no perturbation)...")
	theta = np.pi / 4
	step_circ = build_evolution_step(edges, theta)

	qc = QuantumCircuit(12)
	qc.compose(step_circ, inplace=True)
	qc.compose(step_circ.inverse(), inplace=True)

	initial_state = Statevector.from_label('0' * 12)
	final_state = initial_state.evolve(qc)
	fidelity = np.abs(initial_state.inner(final_state)) ** 2

	assert np.isclose(fidelity, 1.0), f"Sanity check failed: {fidelity}"
	print(f"✓ Sanity Check Passed! U_dagger perfectly inverts U (Fidelity = {fidelity:.1f})\n")


def simulate_echo_fast(edges, max_t, theta, perturb_qubit):
	"""
	Computes the Loschmidt Echo efficiently using the mathematical identity:
	F(t) = |<psi_t | X | psi_t>|^2
	This bypasses Issue 5 (computational waste) completely!
	"""
	echoes = []
	current_state = Statevector.from_label('0' * 12)
	step_circ = build_evolution_step(edges, theta)

	# Circuit to apply the Z perturbation
	qc_z = QuantumCircuit(12)
	qc_z.z(perturb_qubit)

	for t in range(1, max_t + 1):
		# Evolve forward by one step: |psi_t> = U |psi_{t-1}>
		current_state = current_state.evolve(step_circ)

		# Calculate expectation value <psi_t | X | psi_t>
		state_z = current_state.evolve(qc_z)
		expectation_val = np.real(current_state.inner(state_z))

		# Fidelity is the squared expectation value
		fidelity = expectation_val ** 2
		echoes.append(fidelity)

	return echoes


def main():
	edges = build_circlette_lattice()
	run_sanity_check(edges)

	max_t = 20

	# Change 2: Parameterise and sweep theta
	thetas = {
		r'$\pi/16$ (Weak)': np.pi / 16,
		r'$\pi/8$ (Theory)': np.pi / 8,
		r'$\pi/4$ (T-gate Limit)': np.pi / 4,
		r'$\pi/2$ (Clifford Limit)': np.pi / 2
	}

	# Change 3: Vary perturbation location
	targets = {
		0: 'Qubit 0 (Octagon-Square Vertex)',
		2: 'Qubit 2 (Octagon Edge)',
		8: 'Qubit 8 (Square Corner)'
	}

	fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
	fig.suptitle('Circlette Lattice Loschmidt Echo Dynamics', fontsize=16)

	print("Simulating echoes across geometry and rotation angles...")
	for ax, (q_idx, title) in zip(axes, targets.items()):
		ax.set_title(title)
		ax.set_xlabel('Time Depth (t)')
		ax.set_ylim(-0.05, 1.05)
		ax.grid(True, alpha=0.3)
		if ax == axes[0]:
			ax.set_ylabel('Echo Fidelity')

		for theta_label, theta_val in thetas.items():
			fidelities = simulate_echo_fast(edges, max_t, theta_val, q_idx)
			ax.plot(range(1, max_t + 1), fidelities, marker='o', markersize=4, label=f'θ = {theta_label}')

		ax.legend()

	plt.tight_layout()
	plt.show()


if __name__ == "__main__":
	main()