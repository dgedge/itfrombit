import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

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

# ---------------------------------------------------------
# PANEL A: Analytical Band Structure (Bloch Hamiltonian)
# ---------------------------------------------------------
k_vals = np.linspace(-np.pi, np.pi, 200)
bands = np.zeros((4, 200))

for i, k in enumerate(k_vals):
    # Diagonal slice through the Brillouin Zone (kx = ky = k)
    H = np.array([
        [0, 1, np.exp(-1j*k), 1],
        [1, 0, 1, np.exp(1j*k)],
        [np.exp(1j*k), 1, 0, 1],
        [1, np.exp(-1j*k), 1, 0]
    ])
    bands[:, i] = np.sort(np.linalg.eigvalsh(H))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.plot(k_vals, bands[0, :], color='#1f77b4', lw=2.5, label='Slow Branch (Lower)')
ax1.plot(k_vals, bands[1, :], color='#ff7f0e', lw=2.5, label='Fast Branch (Inner)')
ax1.plot(k_vals, bands[2, :], color='#ff7f0e', lw=2.5)
ax1.plot(k_vals, bands[3, :], color='#1f77b4', lw=2.5, label='Slow Branch (Upper)')

ax1.set_title("Analytical Walk Spectrum of 4.8.8 Topology\n(Bifurcated Dispersion Branches)", fontsize=13)
ax1.set_xlabel("Momentum ($k_x = k_y$)", fontsize=12)
ax1.set_ylabel("Energy $E(k)$", fontsize=12)
ax1.set_xticks([-np.pi, 0, np.pi])
ax1.set_xticklabels(['$-\pi$', '$0$', '$\pi$'])
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper right', fontsize=10)

# ---------------------------------------------------------
# PANEL B: Empirical Light Cone (Multi-Velocity Extraction)
# ---------------------------------------------------------
distances = np.array([4, 6, 8, 10, 12])

# Extracted from exact N=3600 sparse continuous-time Krylov evolution
t_sq = np.array([4, 5, 7, 8, 9])
t_hx = np.array([5, 7, 12, 13, 15])
t_488_fast = np.array([5, 7, 9, 13, 14])
t_488_slow = np.array([8, 13, 15, 17, 18])

# Calculate actual velocities (slope = delta_distance / delta_time)
v_sq = linregress(t_sq, distances).slope
v_hx = linregress(t_hx, distances).slope
v_488_f = linregress(t_488_fast, distances).slope
v_488_s = linregress(t_488_slow, distances).slope

ax2.plot(t_sq, distances, 's:', color='gray', ms=7,
         label=f'Square ($z=4$) Baseline [$v={v_sq:.2f}$]')
ax2.plot(t_hx, distances, '^--', color='gray', ms=7,
         label=f'Hex ($z=3$) Control [$v={v_hx:.2f}$]')

ax2.plot(t_488_fast, distances, 'o-', color='#ff7f0e', lw=2.5, ms=8,
         label=f'4.8.8 Fast Wavefront [$v_1={v_488_f:.2f}$]')
ax2.plot(t_488_slow, distances, 'o-', color='#1f77b4', lw=2.5, ms=8,
         label=f'4.8.8 Slow Wavefront [$v_2={v_488_s:.2f}$]')

# Shade the birefringent gap
ax2.fill_betweenx(distances, t_488_fast, t_488_slow, color='indigo', alpha=0.08, label='Lattice Birefringence')

ax2.set_title("Empirical Light Cone: Wavefront Splitting\n(Extracted from N=3600 OTOC dips)", fontsize=13)
ax2.set_xlabel("Arrival Time of Echo Dip ($t$)", fontsize=12)
ax2.set_ylabel("Graph Distance (Shortest-Path Hops)", fontsize=12)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='lower right', fontsize=10)

plt.tight_layout()
plt.savefig("Empirical_Light_Cone.pdf", bbox_inches="tight")
plt.show()