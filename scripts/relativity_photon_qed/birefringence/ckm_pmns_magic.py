import numpy as np
from scipy.linalg import hadamard
import matplotlib.pyplot as plt
# --- Configuration ---
DELTA = 2.0 / 9.0

# Pre-compute the 256x256 Walsh-Hadamard matrix for lightning-fast Pauli expectations
# H_MAT[s, z] = (-1)^(s dot z)
H_MAT = hadamard(256)


def build_walk_operator(delta=DELTA, random_phases=False, theta=None):
	"""Step 1: Build the coherent sum of CNOTs (Kraus-like object).

	Phase mode selection (mutually exclusive):
	  - random_phases=True: each A_k gets uniform random phase in [0, 2π)
	  - theta=<value>:     each A_k = √(δ/7) · exp(ikθ)  -- parametric sweep
	  - default:           physical case, θ = π/4
	"""
	U = np.zeros((256, 256), dtype=np.complex128)

	if random_phases:
		phases = np.random.uniform(0, 2 * np.pi, 7)
	elif theta is not None:
		phases = [k * theta for k in range(1, 8)]
	else:
		phases = [k * np.pi / 4.0 for k in range(1, 8)]

	for k in range(8):
		c = (2 - k) % 8
		t = (5 - k) % 8

		amp = np.sqrt(1.0 - delta) if k == 0 else np.sqrt(delta / 7.0) * np.exp(1j * phases[k - 1])

		for s in range(256):
			s_out = s ^ (1 << t) if (s & (1 << c)) else s
			U[s_out, s] += amp

	return U


def compute_stabiliser_renyi_entropy(v):
	"""
	Step 7: Computes Stabiliser 2-Renyi Entropy.
	Vectorised using the Fast Walsh-Hadamard Transform for ~1ms execution.
	"""
	purity_sum = 0.0
	for x in range(256):
		s_out = np.arange(256) ^ x
		f_x = np.conj(v[s_out]) * v
		# Calculates <v | X^x Z^z | v> for all z simultaneously
		expectations = H_MAT @ f_x
		purity_sum += np.sum(np.abs(expectations) ** 4)

	# E_g = -log2( Sum_{P} |<v|P|v>|^4 / 2^N )
	magic = -np.log2(purity_sum / 256.0)
	return max(0.0, magic)  # Floor at 0 for numerical noise

# tried radius 2 then 1
def get_neighborhood(anchor, radius=1):
	"""Step 4: Build state-space neighbourhoods"""
	return [s for s in range(256) if bin(s ^ anchor).count('1') <= radius]


def run_spectral_analysis(delta=DELTA, random_phases=False):
	# --- Step 1 & 2 ---
	U = build_walk_operator(delta, random_phases)
	U_dag = U.conj().T

	if not random_phases:
		unitarity_dev = np.linalg.norm(U_dag @ U - np.eye(256), ord='fro')
		print(f"[Sanity] ||U†U - I||_F = {unitarity_dev:.4f} (Must be > 0 for δ > 0)")

	M1 = U_dag @ U
	M2 = M1 @ M1

	if not random_phases:
		hermitian_err = np.linalg.norm(M2 - M2.conj().T)
		print(f"[Sanity] Hermiticity error = {hermitian_err:.2e} (Must be ~0)")

	# --- SANITY CHECK: H_up Projection ---
	# TODO: Developer must insert the 9 physical up-quark state integers here
	# up_quarks = [12, 20, 28, 14, 22, 30, 13, 21, 29]
	# Build 3x3 by summing over colour within each generation
	gens = [[12, 20, 28], [13, 21, 29], [14, 22, 30]]
	H_up_3 = np.zeros((3, 3), dtype=complex)
	for i, gen_i in enumerate(gens):
		for j, gen_j in enumerate(gens):
			H_up_3[i, j] = sum(M2[a, b] for a in gen_i for b in gen_j)
	labels = ["Gen 0 (τ, n=0)", "Gen 1 (e, n=1)", "Gen 2 (μ, n=2)"]
	if not random_phases:
		print("\n--- H_up Correctness Check ---")
		H_up_diag = np.diag(H_up_3).real
		ratios = H_up_diag / H_up_diag[0]
		print(f"\nH_up diagonal: {H_up_diag}")
		print(f"Ratios vs Part 04 (target 1.000, 1.014, 1.120): {ratios}")

	#for i, label in enumerate(labels):
	# 	print(f"  {label}: {H_up_3[i, i].real:.4f}")
	#print("  Part 04 target:        1.0710, 1.0860, 1.1990")
	#print(f"  Ratio 2/1: computed {H_up_3[1, 1].real / H_up_3[0, 0].real:.4f}, target 1.0140")
	#print(f"  Ratio 3/1: computed {H_up_3[2, 2].real / H_up_3[0, 0].real:.4f}, target 1.1195")
	#print("Part 04 target (approx):", [1.071, 1.086, 1.199])

	# --- Step 3: Diagonalise ---
	evals, evecs = np.linalg.eigh(M2)
	idx = np.argsort(evals)[::-1]  # Sort descending
	evals = evals[idx]
	evecs = evecs[:, idx]

	# --- Step 4: Define vR Pseudocodeword Anchors ---
	# TODO: Developer must insert actual 8-bit integers for the 3 generation anchors
	anchors = [192, 194, 193]  # <--- PLACEHOLDERS
	neighborhoods = [get_neighborhood(a, radius=1) for a in anchors]   # raius wa 2
	for g, N_g in enumerate(neighborhoods):
		weights = np.array([np.sum(np.abs(evecs[:, j][N_g]) ** 2) for j in range(256)])
		weights /= weights.sum()
		entropies = np.array([compute_stabiliser_renyi_entropy(evecs[:, j]) for j in range(256)])
		mean_E_g = np.sum(weights * entropies)
		# print(f"Gen {g}: ensemble E_g = {mean_E_g:.4f}")

	# --- Step 5, 6 & 8: Extraction ---
	results = []
	if not random_phases:
		print(f"Eigenvalue spectrum summary:")
		print(f"  Max: {evals[0]:.4f}")
		print(f"  Near-1 count: {np.sum(np.abs(evals - 1.0) < 0.1)}")
		print(f"  Min: {evals[-1]:.4f}")
		print(f"  Mean: {np.mean(evals):.4f}")
	for g, N_g in enumerate(neighborhoods):
		best_lam, best_W, best_E = None, 0.0, None

		for j in range(256):
			lam = evals[j].real
			if not (0.8 < lam < 1.2):
				continue

			v = evecs[:, j]
			W_g = np.sum(np.abs(v[N_g]) ** 2)
			others_W = [np.sum(np.abs(v[neighborhoods[other_g]]) ** 2)
			            for other_g in range(3) if other_g != g]

			if W_g > 0.5 and all(w < 0.25 for w in others_W):
				# strict generation-labelling
				if best_lam is None or lam > best_lam:
					best_lam = lam
					best_W = W_g
					best_E = compute_stabiliser_renyi_entropy(v)

		if best_lam is not None:
			if not random_phases:
				print(f"Gen {g}: λ = {best_lam:.5f} | W_g = {best_W:.5f} | E_g (Magic) = {best_E:.5f}")
			results.append((best_lam, best_W, best_E))
		else:
			if not random_phases:
				print(f"Gen {g}: NULL (No localised eigenmode found)")
			results.append(None)

	return results

def compute_ensemble_magic_at_theta(theta, delta=DELTA, anchors=(192, 194, 193), radius=2):
	"""Experiment A helper: compute ensemble E_g for the three generations
	using a parameterised phase exp(ikθ) instead of exp(ikπ/4).

	Returns: (E_gen0, E_gen1, E_gen2) — the three ensemble magic values.
	"""
	U = build_walk_operator(delta=delta, random_phases=False, theta=theta)
	M2 = (U.conj().T @ U) @ (U.conj().T @ U)

	evals, evecs = np.linalg.eigh(M2)
	idx = np.argsort(evals)[::-1]
	evals = evals[idx]
	evecs = evecs[:, idx]

	# Pre-compute entropy for every eigenvector once (reused across all 3 generations)
	entropies = np.array([compute_stabiliser_renyi_entropy(evecs[:, j]) for j in range(256)])

	neighborhoods = [get_neighborhood(a, radius=radius) for a in anchors]

	ensemble_E = []
	for g, N_g in enumerate(neighborhoods):
		weights = np.array([np.sum(np.abs(evecs[:, j][N_g]) ** 2) for j in range(256)])
		weights /= weights.sum()
		mean_E_g = np.sum(weights * entropies)
		ensemble_E.append(mean_E_g)

	return tuple(ensemble_E)

def run_experiment_a():
	"""Experiment A: phase-angle sweep across θ ∈ [0, π].

	For each θ, compute ensemble E_g for all three generations.
	Output: a table that can be tabulated or plotted to show the θ-dependence
	of the gen-0 magic signature.
	"""
	print("\n=== EXPERIMENT A: PHASE-ANGLE SWEEP ===")
	print("Sweeping θ in exp(ikθ) across [0, π] in 13 steps")
	print(f"{'theta':>8s}  {'theta/π':>8s}  {'E_gen0':>8s}  {'E_gen1':>8s}  {'E_gen2':>8s}  notes")
	print("-" * 70)

	theta_values = [
		0.0,
		np.pi / 16,
		np.pi / 8,
		3 * np.pi / 16,
		np.pi / 4,        # physical case
		5 * np.pi / 16,
		3 * np.pi / 8,
		7 * np.pi / 16,
		np.pi / 2,
		5 * np.pi / 8,
		3 * np.pi / 4,
		7 * np.pi / 8,
		np.pi,
	]

	results = []
	for theta in theta_values:
		E_g = compute_ensemble_magic_at_theta(theta)
		results.append((theta, E_g))

		note = " <-- physical (π/4)" if abs(theta - np.pi / 4) < 1e-9 else ""
		print(f"{theta:>8.4f}  {theta/np.pi:>8.4f}  "
		      f"{E_g[0]:>8.4f}  {E_g[1]:>8.4f}  {E_g[2]:>8.4f}{note}")

	# Identify which θ gives peak E_gen0 — informative for interpretation
	gen0_values = [r[1][0] for r in results]
	peak_idx = int(np.argmax(gen0_values))
	peak_theta = results[peak_idx][0]
	print(f"\nPeak gen-0 magic occurs at θ = {peak_theta:.4f} ({peak_theta/np.pi:.4f}·π)")
	print(f"Peak gen-0 E_g = {gen0_values[peak_idx]:.4f}")

	if abs(peak_theta - np.pi / 4) < np.pi / 16:
		print("→ Peak is at or near the physical π/4 value: phase-driven interpretation supported")
	else:
		print("→ Peak is well away from π/4: the physical case is not at a phase resonance")

	plt.figure(figsize=(8, 5))
	plt.plot([r[0] / np.pi for r in results], [r[1][0] for r in results], 'o-', label='Gen 0 (τ)')
	plt.plot([r[0] / np.pi for r in results], [r[1][1] for r in results], 's-', label='Gen 1 (e)')
	plt.plot([r[0] / np.pi for r in results], [r[1][2] for r in results], '^-', label='Gen 2 (μ)')
	plt.axvline(0.25, color='red', linestyle='--', alpha=0.5, label='Physical θ = π/4')
	plt.xlabel('θ / π')
	plt.ylabel('Ensemble E_g (magic)')
	plt.title('Phase-angle sweep: ensemble magic per generation')
	plt.legend()
	plt.grid(True, alpha=0.3)
	plt.savefig('experiment_a_phase_sweep.png', dpi=150, bbox_inches='tight')
	plt.show()
	return results

def diagnostic_1_basis_state_magic():
	"""Diagnostic 1: compute the magic of each basis state and tabulate
	by νR neighbourhood. Tests whether the random-phase gen-0 collapse
	is an encoding artefact (gen-0 basis states are individually low-magic)
	or something emergent (basis states are similar but eigenmode weighting
	differs)."""

	print("\n=== DIAGNOSTIC 1: PER-BASIS-STATE MAGIC ===")

	anchors = [192, 194, 193]
	labels = ["Gen 0 (τ, anchor 192)", "Gen 1 (e, anchor 194)", "Gen 2 (μ, anchor 193)"]
	neighborhoods = [get_neighborhood(a, radius=2) for a in anchors]

	# Compute magic of every single basis state |s⟩
	# A basis state has all amplitude on one position, so its magic is
	# computable from compute_stabiliser_renyi_entropy on a one-hot vector
	basis_state_magic = np.zeros(256)
	for s in range(256):
		v = np.zeros(256, dtype=np.complex128)
		v[s] = 1.0
		basis_state_magic[s] = compute_stabiliser_renyi_entropy(v)

	# Spoiler check: pure basis states are stabiliser states by definition,
	# so their magic should all be exactly 0. If this returns non-zero values,
	# the formula or normalisation has a problem (which would be Insight 4).
	print(f"Sanity check: max magic of a basis state = {basis_state_magic.max():.6f}")
	print(f"             (should be 0 — basis states are stabilisers)")
	print(f"             min = {basis_state_magic.min():.6f}, mean = {basis_state_magic.mean():.6f}")

	# That tells us whether Insight 4 (normalisation artefact) is in play.
	# But basis states being all stabiliser doesn't tell us about the
	# *neighbourhood structure* under random superposition. So let's also
	# do the more interesting test:

	print("\n--- Magic of random superpositions within each neighbourhood ---")
	print("(For each neighbourhood, draw 200 random unit vectors supported")
	print(" only on that neighbourhood, compute their average magic.)")
	print(f"{'Neighbourhood':<32s}  {'mean magic':>12s}  {'std':>10s}  {'size':>5s}")
	print("-" * 70)

	rng = np.random.default_rng(42)  # reproducible
	for g, (label, N_g) in enumerate(zip(labels, neighborhoods)):
		entropies_in_neighbourhood = []
		for _ in range(200):
			# Random complex vector supported only on this neighbourhood
			v = np.zeros(256, dtype=np.complex128)
			amps = rng.standard_normal(len(N_g)) + 1j * rng.standard_normal(len(N_g))
			amps /= np.linalg.norm(amps)
			v[N_g] = amps
			entropies_in_neighbourhood.append(compute_stabiliser_renyi_entropy(v))

		mean_e = np.mean(entropies_in_neighbourhood)
		std_e = np.std(entropies_in_neighbourhood)
		print(f"{label:<32s}  {mean_e:>12.4f}  {std_e:>10.4f}  {len(N_g):>5d}")

	print("\nReading guide:")
	print("  - If the three neighbourhoods give similar mean magic (within 10-20%):")
	print("    the gen-0 random-phase collapse is NOT an encoding-level property.")
	print("    It must be an emergent feature of the eigenmode weighting under")
	print("    the walk dynamics — Insight 2 or 3 territory.")
	print("  - If the gen-0 neighbourhood gives substantially lower mean magic")
	print("    (factor of 2 or more) than the other two:")
	print("    the encoding-level Pauli structure of gen-0 IS asymmetric.")
	print("    This would be Insight 1 — small but real new observation.")


def diagnostic_2b_eigenmode_structure():
	"""Diagnostic 2B: examine the structure of representative eigenmodes
	from each generation neighbourhood. Distinguishes 'concentrated
	amplitude on few basis states' (combinatorial) from 'spread amplitude
	with structured cancellation' (emergent)."""

	print("\n=== DIAGNOSTIC 2B: EIGENMODE STRUCTURE ===")
	print("For each generation, find the eigenmode most concentrated in its")
	print("νR neighbourhood and inspect where its amplitude actually lives.\n")

	# Build the physical walk and its M² operator
	U = build_walk_operator(delta=DELTA, random_phases=False)
	M2 = (U.conj().T @ U) @ (U.conj().T @ U)
	evals, evecs = np.linalg.eigh(M2)
	idx = np.argsort(evals)[::-1]
	evals = evals[idx]
	evecs = evecs[:, idx]

	anchors = [192, 194, 193]
	labels = ["Gen 0 (τ, anchor 192)", "Gen 1 (e, anchor 194)", "Gen 2 (μ, anchor 193)"]
	neighborhoods = [get_neighborhood(a, radius=2) for a in anchors]

	for g, (label, N_g) in enumerate(zip(labels, neighborhoods)):
		# Find the eigenmode with the largest weight in this neighbourhood
		weights = np.array([np.sum(np.abs(evecs[:, j][N_g]) ** 2) for j in range(256)])
		best_j = int(np.argmax(weights))
		v = evecs[:, best_j]

		# Concentration metric: participation ratio
		# PR = 1 / Σ|v_s|⁴, ranges from 1 (single basis state) to 256 (uniform)
		probs = np.abs(v) ** 2
		participation_ratio = 1.0 / np.sum(probs ** 2)

		# Magic of this eigenmode
		magic = compute_stabiliser_renyi_entropy(v)

		# Top-10 basis states by amplitude
		top_idx = np.argsort(probs)[::-1][:10]

		print(f"--- {label} ---")
		print(f"  Best eigenmode index:    {best_j}")
		print(f"  Eigenvalue λ:            {evals[best_j]:.4f}")
		print(f"  Weight in neighbourhood: {weights[best_j]:.4f}")
		print(f"  Magic (2-Rényi):         {magic:.4f}")
		print(f"  Participation ratio:     {participation_ratio:.2f}")
		print(f"    (1 = single state, 256 = uniformly spread)")
		print(f"  Top 10 basis states by probability:")
		print(f"    {'state':>5s}  {'binary':>8s}  {'|amp|²':>10s}  {'in N_g?':>7s}")
		for s in top_idx:
			in_ng = "yes" if s in N_g else "no"
			print(f"    {s:>5d}  {format(s, '08b'):>8s}  {probs[s]:>10.6f}  {in_ng:>7s}")
		print()

	print("Reading guide:")
	print("  Combinatorial (Insight 1-ish):")
	print("    - Gen 0 eigenmode has LOW participation ratio (e.g. 1-3)")
	print("    - Amplitude concentrates on 1-2 basis states")
	print("    - Low magic explained by 'almost a basis state'")
	print()
	print("  Structured cancellation (Insight 2/3):")
	print("    - Gen 0 eigenmode has HIGH participation ratio (e.g. 10+)")
	print("    - Amplitude spread across many basis states")
	print("    - Low magic despite spread → structural Pauli cancellation")
	print("    - Compare against Gen 1 and Gen 2 to see if their structure differs")


def diagnostic_2c_eigenmode_contributions():
	"""Print the weight × magic contribution of every eigenmode to gen-0 ensemble."""
	U = build_walk_operator(delta=DELTA, random_phases=True)
	M2 = (U.conj().T @ U) @ (U.conj().T @ U)
	evals, evecs = np.linalg.eigh(M2)

	N_0 = get_neighborhood(192, radius=2)

	weights = np.array([np.sum(np.abs(evecs[:, j][N_0]) ** 2) for j in range(256)])
	weights /= weights.sum()
	magics = np.array([compute_stabiliser_renyi_entropy(evecs[:, j]) for j in range(256)])
	contributions = weights * magics

	# Sort by contribution and show top 20
	idx = np.argsort(weights)[::-1]
	print(f"\n=== DIAGNOSTIC 2C: Gen-0 contributions (random-phase, single draw) ===")
	print(f"Total ensemble magic: {np.sum(contributions):.4f}")
	print(f"\nTop 20 eigenmodes by weight in N_0:")
	print(f"{'idx':>4s}  {'eval':>8s}  {'weight':>8s}  {'magic':>8s}  {'contrib':>8s}")
	for j in idx[:20]:
		print(f"{j:>4d}  {evals[j]:>8.4f}  {weights[j]:>8.4f}  {magics[j]:>8.4f}  {contributions[j]:>8.4f}")


def diagnostic_2d_radius_1_baseline():
	"""Re-run random-phase baseline at radius=1 instead of radius=2."""
	import sys
	# Same as random-phase baseline but with radius=1 in get_neighborhood calls
	run_spectral_analysis(delta=DELTA, random_phases=False)

if __name__ == "__main__":
	diagnostic_2c_eigenmode_contributions()

'''
if __name__ == "__main__":
	# ... your existing code ...
	diagnostic_2d_radius_1_baseline()



if __name__ == "__main__":
	# ... your existing code ...
	diagnostic_2b_eigenmode_structure()

if __name__ == "__main__":
	# ... your existing code ...
	diagnostic_1_basis_state_magic()


if __name__ == "__main__":
	# Existing runs
	print("=== PROJECT ECHO: NULL BASELINE (δ = 0) ===")
	run_spectral_analysis(delta=0.0)

	print("\n=== PROJECT ECHO: PHYSICAL LATTICE (δ = 2/9) ===")
	run_spectral_analysis(delta=DELTA)

	# (Random-phase baseline removed for brevity in this run; uncomment when needed)
	# print("\n=== PROJECT ECHO: RANDOM PHASE STATISTICAL BASELINE ===")
	# ...

	# New: Experiment A
	results_a = run_experiment_a()

'''
'''
if __name__ == "__main__":
	print("=== PROJECT ECHO: NULL BASELINE (δ = 0) ===")
	run_spectral_analysis(delta=0.0)

	print("\n=== PROJECT ECHO: PHYSICAL LATTICE (δ = 2/9) ===")
	run_spectral_analysis(delta=DELTA)

	print("\n=== PROJECT ECHO: RANDOM PHASE STATISTICAL BASELINE ===")
	entropies = {0: [], 1: [], 2: []}
	for _ in range(1000):
		res = run_spectral_analysis(delta=DELTA, random_phases=True)
		for g in range(3):
			if res[g] is not None:
				entropies[g].append(res[g][2])
	physical_E_g = [3.9939, 4.1491, 4.2016]
	for g in range(3):
		mean_E = np.mean(entropies[g]) if entropies[g] else 0.0
		std_E = np.std(entropies[g])
		z = (physical_E_g[g] - mean_E) / std_E
		print(f"Gen {g}: μ={mean_E:.3f} σ={std_E:.3f} z={z:.1f}")
		print(f"Gen {g} Random Phase Mean E_g: {mean_E:.3f}")
'''