"""
Quantum walk wavefront on 2D lattices.

Renders a continuous-time quantum walk
    i d_t psi = A psi
on one of three lattices (4.8.8, honeycomb, square) starting from a
delta-function at the centre, and saves it as an MP4.

Time evolution is done with scipy.sparse.linalg.expm_multiply (Krylov),
so it is strictly unitary up to numerical precision -- not Trotterised.

Knobs are at the top. Run with: python wavefront.py
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

# ---------------------------------------------------------------------------
# Knobs
# ---------------------------------------------------------------------------

LATTICE   = "4.8.8"          # "4.8.8" | "hex" | "square"
NX, NY    = 24, 24           # unit cells (or grid dimensions for square)
T_MAX     = 18.0             # final time
N_FRAMES  = 220              # number of frames in the animation
FPS       = 25               # video frame rate
COLOURING = "intensity"      # "intensity" (|psi|^2) or "phase" (Re psi)
OUTPUT    = "wavefront.mp4"  # output path

# ---------------------------------------------------------------------------
# Lattice builders -- each returns (vertices, edges)
#   vertices : list of (x, y)
#   edges    : set of (i, j) with i < j
# ---------------------------------------------------------------------------

def build_4_8_8(nx, ny, L=1.0, R=0.32):
    """Truncated square tiling. 4 vertices per unit cell arranged as a
    diamond at each grid point; each diamond connects to its right and
    upper neighbour via single octagon edges."""
    KOFF = [(0, R), (R, 0), (0, -R), (-R, 0)]   # N, E, S, W
    verts, idx = [], {}
    for i in range(nx):
        for j in range(ny):
            for k in range(4):
                idx[(i, j, k)] = len(verts)
                verts.append((i * L + KOFF[k][0], j * L + KOFF[k][1]))
    edges = set()
    def add(a, b):
        if a is None or b is None or a == b: return
        edges.add((min(a, b), max(a, b)))
    for i in range(nx):
        for j in range(ny):
            n_, e_, s_, w_ = (idx[(i, j, k)] for k in range(4))
            # intra-cell square
            add(n_, e_); add(e_, s_); add(s_, w_); add(w_, n_)
            # inter-cell octagon edges
            if i + 1 < nx: add(e_, idx[(i + 1, j, 3)])
            if j + 1 < ny: add(n_, idx[(i, j + 1, 2)])
    return verts, edges


def build_hex(nx, ny):
    """Honeycomb lattice. Two atoms per unit cell, bond length 1."""
    a1 = np.array([1.5, np.sqrt(3) / 2])
    a2 = np.array([1.5, -np.sqrt(3) / 2])
    verts, idx = [], {}
    for i in range(nx):
        for j in range(ny):
            base = i * a1 + j * a2
            idx[(i, j, 0)] = len(verts); verts.append(tuple(base))
            idx[(i, j, 1)] = len(verts); verts.append(tuple(base + np.array([1.0, 0.0])))
    edges = set()
    def add(a, b):
        if a is None or b is None or a == b: return
        edges.add((min(a, b), max(a, b)))
    for i in range(nx):
        for j in range(ny):
            A_ = idx[(i, j, 0)]
            B_ = idx[(i, j, 1)]
            add(A_, B_)                                      # in-cell bond
            if i - 1 >= 0:        add(A_, idx[(i - 1, j, 1)])  # to B of (i-1, j)
            if j - 1 >= 0:        add(A_, idx[(i, j - 1, 1)])  # to B of (i, j-1)
    return verts, edges


def build_square(nx, ny):
    """Plain square lattice."""
    verts, idx = [], {}
    for i in range(nx):
        for j in range(ny):
            idx[(i, j)] = len(verts)
            verts.append((float(i), float(j)))
    edges = set()
    for i in range(nx):
        for j in range(ny):
            v = idx[(i, j)]
            if i + 1 < nx: edges.add((v, idx[(i + 1, j)]))
            if j + 1 < ny: edges.add((v, idx[(i, j + 1)]))
    return verts, edges


BUILDERS = {"4.8.8": build_4_8_8, "hex": build_hex, "square": build_square}

# ---------------------------------------------------------------------------
# Build lattice + adjacency matrix
# ---------------------------------------------------------------------------

verts, edges = BUILDERS[LATTICE](NX, NY)
N = len(verts)
xs = np.array([v[0] for v in verts])
ys = np.array([v[1] for v in verts])
print(f"Lattice {LATTICE}: {N} vertices, {len(edges)} edges")

rows, cols = [], []
for a, b in edges:
    rows += [a, b]; cols += [b, a]
A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N))

# Pick origin nearest to geometric centre
cx, cy = xs.mean(), ys.mean()
origin = int(np.argmin((xs - cx) ** 2 + (ys - cy) ** 2))

# ---------------------------------------------------------------------------
# Time evolution: psi(t) = exp(-i A t) psi0 via Krylov
# ---------------------------------------------------------------------------

psi0 = np.zeros(N, dtype=complex)
psi0[origin] = 1.0

print("Computing time evolution...")
op = (-1j) * A.astype(complex)
snapshots = expm_multiply(op, psi0, start=0.0, stop=T_MAX,
                          num=N_FRAMES, endpoint=True)
times = np.linspace(0.0, T_MAX, N_FRAMES)

# Quick unitarity sanity check
norms = np.linalg.norm(snapshots, axis=1)
print(f"Norm drift: max |1 - ||psi|| | = {abs(1 - norms).max():.2e}")

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

edge_segments = [[(xs[a], ys[a]), (xs[b], ys[b])] for (a, b) in edges]

if COLOURING == "intensity":
    cmap = LinearSegmentedColormap.from_list("wf", [
        (0.00, (0.31, 0.22, 0.71)),
        (0.50, (1.00, 0.39, 0.16)),
        (1.00, (1.00, 0.98, 0.85)),
    ])
else:  # phase: diverging cool/warm
    cmap = plt.get_cmap("RdBu_r")

fig, ax = plt.subplots(figsize=(8, 8), dpi=110, facecolor="#14101f")
ax.set_facecolor("#14101f")
margin = 0.5
ax.set_xlim(xs.min() - margin, xs.max() + margin)
ax.set_ylim(ys.min() - margin, ys.max() + margin)
ax.set_aspect("equal"); ax.axis("off")

ax.add_collection(LineCollection(edge_segments,
                  colors=(0.55, 0.51, 0.78, 0.13), linewidths=0.4))
scat = ax.scatter(xs, ys, s=4, c="#3a3260", zorder=3, edgecolors="none")

t_text = ax.text(0.02, 0.97, "", transform=ax.transAxes, color="#dcd7f0",
                 family="monospace", fontsize=12, va="top")
ax.text(0.98, 0.97, f"{LATTICE}  ·  N = {N}", transform=ax.transAxes,
        color="#dcd7f0", family="monospace", fontsize=12, va="top", ha="right")


def update(frame):
    psi = snapshots[frame]
    amp2 = (psi.real ** 2 + psi.imag ** 2)
    m = max(amp2.max(), 1e-12)

    if COLOURING == "intensity":
        x = np.minimum(1.0, (amp2 / m) ** 0.45)
        sizes = 4 + x * 90
        rgba = cmap(x)
        rgba[:, 3] = np.where(x < 0.04, 0.0, 0.45 + x * 0.50)
    else:  # phase
        intensity = np.minimum(1.0, (amp2 / m) ** 0.45)
        signed = psi.real / max(np.abs(psi.real).max(), 1e-12)
        sizes = 4 + intensity * 90
        rgba = cmap(0.5 + 0.5 * signed)
        rgba[:, 3] = np.where(intensity < 0.04, 0.0, 0.45 + intensity * 0.50)

    scat.set_sizes(sizes)
    scat.set_facecolors(rgba)
    t_text.set_text(f"t = {times[frame]:5.2f}")
    return scat, t_text


print("Encoding MP4...")
anim = animation.FuncAnimation(fig, update, frames=N_FRAMES,
                               interval=1000 // FPS, blit=True)
writer = animation.FFMpegWriter(fps=FPS, bitrate=2600, codec="h264")
anim.save(OUTPUT, writer=writer, savefig_kwargs={"facecolor": "#14101f"})
print(f"Wrote {OUTPUT}")
