import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from collections import deque

plt.rcParams.update({
    "pdf.fonttype": 42,          # TrueType, not Type 3 — required by most journals and arXiv
    "ps.fonttype": 42,
    "font.family": "serif",
    "font.size": 9,              # 9–10 pt is standard for journal figures
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.2,
    "figure.figsize": (3.4, 2.6), # single-column width for most physics journals
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

def build_square(n):
	edges = set()

	def idx(x, y):
		return y * n + x

	for y in range(n):
		for x in range(n):
			if x < n - 1: edges.add((idx(x, y), idx(x + 1, y)))
			if y < n - 1: edges.add((idx(x, y), idx(x, y + 1)))
	return list(edges), n * n, idx(n // 2, n // 2)


def build_hex(n):
	edges = set()

	def idx(x, y):
		return y * n + x

	for y in range(n):
		for x in range(n):
			if x < n - 1: edges.add((idx(x, y), idx(x + 1, y)))
			if y < n - 1 and (x + y) % 2 == 0: edges.add((idx(x, y), idx(x, y + 1)))
	return list(edges), n * n, idx(n // 2, n // 2)


def build_488(grid_size):
	edges = set()

	def idx(x, y, i):
		return (y * grid_size + x) * 4 + i

	for y in range(grid_size):
		for x in range(grid_size):
			u0, u1, u2, u3 = idx(x, y, 0), idx(x, y, 1), idx(x, y, 2), idx(x, y, 3)
			edges.add((u0, u1));
			edges.add((u1, u2));
			edges.add((u2, u3));
			edges.add((u3, u0))
			if x < grid_size - 1: edges.add((u1, idx(x + 1, y, 3)))
			if y < grid_size - 1: edges.add((u2, idx(x, y + 1, 0)))
	return list(edges), grid_size * grid_size * 4, idx(grid_size // 2, grid_size // 2, 0)


def get_echo(edges, N, start_node, max_t=25, theta=np.pi / 8, d_target=6):
	A = np.zeros((N, N))
	adj = {i: [] for i in range(N)}
	for u, v in edges:
		A[u, v] = A[v, u] = 1.0
		adj[u].append(v)
		adj[v].append(u)

	# BFS to find tripwires exactly d=6 hops away
	dist = {start_node: 0}
	q = deque([start_node])
	while q:
		curr = q.popleft()
		for nbr in adj[curr]:
			if nbr not in dist:
				dist[nbr] = dist[curr] + 1
				q.append(nbr)

	tripwires = [n for n, d in dist.items() if d == d_target]
	U_step = expm(-1j * theta * A)  # Continuous time exact operator

	echoes = []
	state_t = np.zeros(N, dtype=complex)
	state_t[start_node] = 1.0

	for t in range(1, max_t + 1):
		state_t = U_step @ state_t
		t_echoes = []
		for tw in tripwires:
			# The mathematical shortcut for perfectly reversing a 1-particle phase flip
			p_tw = np.abs(state_t[tw]) ** 2
			t_echoes.append((1 - 2 * p_tw) ** 2)
		echoes.append(np.mean(t_echoes))

	return echoes, len(tripwires)


def main():
	print("Building lattices and running Continuous-Time Echoes...")
	sq_edges, sq_N, c_sq = build_square(30)
	hx_edges, hx_N, c_hx = build_hex(30)
	G48_edges, G48_N, c_48 = build_488(15)

	echo_sq, count_sq = get_echo(sq_edges, sq_N, c_sq)
	echo_hx, count_hx = get_echo(hx_edges, hx_N, c_hx)
	echo_48, count_48 = get_echo(G48_edges, G48_N, c_48)

	plt.figure(figsize=(10, 6))
	plt.title("Quantum Butterfly Echo: Isolating Multi-Velocity Dispersion\n(Continuous Time Walk, Averaged over all tripwires exactly d=6 hops away)", fontsize=12)
	plt.plot(range(1, 26), echo_sq, label=f'Square (z=4, avg over {count_sq} nodes)', marker='s', color='#1f77b4', linestyle=':')
	plt.plot(range(1, 26), echo_hx, label=f'Hex (z=3, avg over {count_hx} nodes)', marker='^', color='#ff7f0e', linestyle='--')
	plt.plot(range(1, 26), echo_48, label=f'4.8.8 (z=3, avg over {count_48} nodes)', marker='o', color='#2ca02c', linestyle='-', linewidth=2)

	plt.xlabel("Time Depth (t)")
	plt.ylabel("Return Fidelity (Echo)")
	plt.axhline(1.0, color='k', linestyle='-', alpha=0.2)
	plt.legend(loc='lower right')
	plt.grid(True, alpha=0.3)
	plt.tight_layout()
	plt.savefig("Quantum_Butterfly_Echo.pdf", bbox_inches="tight")
	plt.show()

if __name__ == "__main__":
	main()