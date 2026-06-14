import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector


def build_488_lattice():
	"""Experiment: The 4.8.8 Circlette Patch"""
	octagon = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0)]
	square1 = [(0, 8), (8, 9), (9, 1)]
	square2 = [(4, 10), (10, 11), (11, 5)]
	return octagon + square1 + square2


def build_1d_ring_control():
	"""Control: A generic 12-qubit 1D ring (Perfect symmetry)."""
	return [(i, (i + 1) % 12) for i in range(12)]


def build_hopping_step(edges, theta):
	"""Pure XX+YY hopping to cleanly probe the adjacency matrix."""
	qc = QuantumCircuit(12)
	for edge in edges:
		qc.rxx(theta, edge[0], edge[1])
		qc.ryy(theta, edge[0], edge[1])
	return qc


def simulate_single_particle_walk(edges, max_t, theta, start_node):
	"""Injects EXACTLY ONE particle and measures return probability."""
	qc_init = QuantumCircuit(12)
	qc_init.x(start_node)
	current_state = Statevector.from_instruction(qc_init)

	step_circ = build_hopping_step(edges, theta)

	qc_z = QuantumCircuit(12)
	qc_z.z(start_node)

	probabilities = []
	for t in range(1, max_t + 1):
		current_state = current_state.evolve(step_circ)

		# Calculate <psi | Z | psi>
		state_z = current_state.evolve(qc_z)
		exp_z = np.real(current_state.inner(state_z))

		# Convert Z expectation (+1 for |0>, -1 for |1>) to probability of finding the particle
		prob_1 = (1.0 - exp_z) / 2.0
		probabilities.append(prob_1)

	return probabilities


def main():
	edges_488 = build_488_lattice()
	edges_1d = build_1d_ring_control()

	max_t = 30
	theta_val = np.pi / 8

	# We will inject at node 8 (Square Corner) and measure the geometric echoes
	target_node = 8

	fig, ax = plt.subplots(figsize=(12, 6))
	ax.set_title(f'Geometric Specificity: 4.8.8 vs 1D Ring Control\n(Single-Particle Injection & Measurement at Node {target_node})', fontsize=14)
	ax.set_xlabel('Time Depth (t)')
	ax.set_ylabel('Return Probability (Echo)')
	ax.set_ylim(-0.05, 1.05)
	ax.grid(True, alpha=0.3)

	# Run Generic 1D Ring Control
	# (Note: Node 8 on a 1D ring is geometrically identical to any other node)
	prob_1d = simulate_single_particle_walk(edges_1d, max_t, theta_val, start_node=target_node)
	ax.plot(range(1, max_t + 1), prob_1d, marker='x', markersize=5,
	        color='gray', linestyle=':', label='Control: Generic 1D Ring')

	# Run 4.8.8 Circlette Lattice
	prob_488 = simulate_single_particle_walk(edges_488, max_t, theta_val, start_node=target_node)
	ax.plot(range(1, max_t + 1), prob_488, marker='o', markersize=5,
	        color='indigo', linestyle='-', label='Experiment: 4.8.8 Circlette (Square Corner)')

	ax.legend()
	plt.tight_layout()
	plt.show()


if __name__ == "__main__":
	main()