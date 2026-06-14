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


def build_unconstrained_step(edges, theta, num_qubits=12):
	"""The original, violently chaotic circuit."""
	qc = QuantumCircuit(num_qubits)
	for i in range(num_qubits):
		qc.rx(theta, i)
		qc.rz(theta, i)
	for edge in edges:
		qc.cx(edge[0], edge[1])
	return qc


def build_constrained_step(edges, theta, num_qubits=12):
	"""
	The Z-Magnetization Conserving Circuit.
	Simulates a particle-hopping tight-binding model on the lattice.
	"""
	qc = QuantumCircuit(num_qubits)
	# Local chemical potential (conserves Z)
	for i in range(num_qubits):
		qc.rz(theta, i)
	# XX+YY hopping strictly conserves the number of excitations
	for edge in edges:
		qc.rxx(theta, edge[0], edge[1])
		qc.ryy(theta, edge[0], edge[1])
	return qc


def simulate_head_to_head(edges, max_t, theta, perturb_qubit, is_constrained):
	"""Simulates the Loschmidt Echo for either architecture."""
	# STATE PREP: We need particles for the constrained model to hop!
	# Initialize a half-filled lattice (6 particles on alternating nodes)
	qc_init = QuantumCircuit(12)
	qc_init.x([0, 2, 4, 6, 8, 10])
	current_state = Statevector.from_instruction(qc_init)

	if is_constrained:
		step_circ = build_constrained_step(edges, theta)
	else:
		step_circ = build_unconstrained_step(edges, theta)

	# The Butterfly Signal (Z-gate phase kick)
	qc_z = QuantumCircuit(12)
	qc_z.z(perturb_qubit)

	echoes = []

	for t in range(1, max_t + 1):
		# Evolve forward
		current_state = current_state.evolve(step_circ)

		# Calculate fast expectation value of Z
		state_z = current_state.evolve(qc_z)
		expectation_val = np.real(current_state.inner(state_z))

		# Fidelity = squared expectation value
		fidelity = expectation_val ** 2
		echoes.append(fidelity)

	return echoes


def main():
	edges = build_circlette_lattice()
	max_t = 25  # Extended time depth to watch structural dynamics

	# We will test your theoretical target angle
	theta_val = np.pi / 8

	targets = {
		0: 'Qubit 0 (Octagon-Square Vertex)',
		2: 'Qubit 2 (Octagon Edge)',
		8: 'Qubit 8 (Square Corner)'
	}

	fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
	fig.suptitle(r'Circlette Lattice Echo Dynamics ($\theta = \pi/8$): Unconstrained vs Conserved Z', fontsize=16)

	print("Running Head-to-Head dynamics (This will take a few seconds)...")
	for ax, (q_idx, title) in zip(axes, targets.items()):
		ax.set_title(title)
		ax.set_xlabel('Time Depth (t)')
		ax.set_ylim(-0.05, 1.05)
		ax.grid(True, alpha=0.3)
		if ax == axes[0]:
			ax.set_ylabel('Echo Fidelity')

		# 1. Run Unconstrained (The baseline chaotic decay)
		fid_un = simulate_head_to_head(edges, max_t, theta_val, q_idx, is_constrained=False)
		ax.plot(range(1, max_t + 1), fid_un, marker='x', markersize=4,
		        color='gray', linestyle=':', label='Unconstrained (Chaotic)')

		# 2. Run Constrained (The structured codeword model)
		fid_con = simulate_head_to_head(edges, max_t, theta_val, q_idx, is_constrained=True)
		ax.plot(range(1, max_t + 1), fid_con, marker='o', markersize=4,
		        color='indigo', linestyle='-', label='Constrained (Particle Hopping)')

		ax.legend()

	plt.tight_layout()
	plt.show()


if __name__ == "__main__":
	main()