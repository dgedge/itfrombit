#!/usr/bin/env python3
"""3D geometric picture: colour = orientation of an oblate square bipyramid, and a
proton = three orthogonal colour-bipyramids meeting at a shared vertex (no overlap).
Renders a two-panel PNG. Run under ~/bin/py13_7 (matplotlib)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

A = 1.0
AX = {"x": np.array([1, 0, 0.]), "y": np.array([0, 1, 0.]), "z": np.array([0, 0, 1.])}
PERP = {"x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y")}
COL = {"x": "#d23b3b", "y": "#2e9e44", "z": "#3b6fd2"}  # red / green / blue = the 3 colours


def bipyramid(P0, axis, a=A):
    """8-face oblate square bipyramid: inner apex at P0, outer apex at P0+a*dir,
    square waist (edge a) at the midpoint, in the plane perpendicular to the axis."""
    d = AX[axis]; u, v = AX[PERP[axis][0]], AX[PERP[axis][1]]
    m = P0 + 0.5 * a * d
    w = [m + 0.5 * a * (su * u + sv * v) for su, sv in [(1, 1), (1, -1), (-1, -1), (-1, 1)]]
    Ain, Aout = P0, P0 + a * d
    faces = [[Ain, w[i], w[(i + 1) % 4]] for i in range(4)]
    faces += [[Aout, w[i], w[(i + 1) % 4]] for i in range(4)]
    return faces


def add_bipyramid(ax, P0, axis, a=A, alpha=0.55):
    pc = Poly3DCollection(bipyramid(P0, axis, a), facecolor=COL[axis], edgecolor="#222",
                          linewidths=0.6, alpha=alpha)
    ax.add_collection3d(pc)


def add_c3_circuit(ax, O=np.zeros(3), t=0.42, R=0.66):
    """Body-diagonal C3 axis + the r->g->b circulation ring with small arrowheads."""
    n = np.array([1, 1, 1.]) / np.sqrt(3)
    e1 = np.array([1, -1, 0.]) / np.sqrt(2)
    e2 = np.cross(n, e1); e2 /= np.linalg.norm(e2)
    cen = O + t * n
    ax.plot(*zip(O - 0.6 * n, O + 1.7 * n), color="#666", lw=1.2, ls=(0, (5, 3)))  # C3 axis
    ang = lambda d: np.arctan2((d - (d @ n) * n) @ e2, (d - (d @ n) * n) @ e1)
    angs = {k: ang(AX[k]) for k in "xyz"}
    pt = lambda a: cen + R * (np.cos(a) * e1 + np.sin(a) * e2)
    for a, b in [("x", "y"), ("y", "z"), ("z", "x")]:
        a0 = angs[a]; d = (angs[b] - a0 + np.pi) % (2 * np.pi) - np.pi
        aa = np.linspace(a0 + 0.14 * d, a0 + 0.80 * d, 24); arc = np.array([pt(x) for x in aa])
        ax.plot(arc[:, 0], arc[:, 1], arc[:, 2], color=COL[a], lw=2.8)
        v = arc[-1] - arc[-2]; v = v / np.linalg.norm(v)
        ax.quiver(*arc[-1], *v, color=COL[a], length=0.17, arrow_length_ratio=0.6, lw=2.8)
    ax.text2D(0.5, 0.66, r"C$_3$ axis (body diagonal):  r $\to$ g $\to$ b",
              transform=ax.transAxes, ha="center", fontsize=10.5, color="#333")


def draw_cube(ax, c=np.zeros(3), s=A, **kw):
    r = s / 2
    pts = np.array([[x, y, z] for x in (-r, r) for y in (-r, r) for z in (-r, r)]) + c
    edges = [(0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)]
    for i, j in edges:
        ax.plot(*zip(pts[i], pts[j]), color="#9aa0a6", lw=0.8, ls=(0, (4, 3)), **kw)


def setup(ax, title):
    O = np.zeros(3)
    for axis in "xyz":
        ax.plot(*zip(O, 1.35 * A * AX[axis]), color=COL[axis], lw=1.4)
        p = 1.5 * A * AX[axis]
        ax.text(*p, f"{axis}", color=COL[axis], fontsize=12, fontweight="bold", ha="center", va="center")
    draw_cube(ax)
    ax.scatter(*O, color="k", s=22, depthshade=False)
    ax.set_xlim(-0.7, 1.4); ax.set_ylim(-0.7, 1.4); ax.set_zlim(-0.7, 1.4)
    ax.set_box_aspect([1, 1, 1]); ax.view_init(elev=20, azim=34)
    ax.set_axis_off(); ax.set_title(title, fontsize=12.5, pad=2)


fig = plt.figure(figsize=(12.6, 6.4))
fig.suptitle("Colour = orientation: how three quarks build a proton  (geometric reading of the framework)",
             fontsize=14, fontweight="bold", y=0.97)

# panel 1: single quark
ax1 = fig.add_subplot(121, projection="3d")
add_bipyramid(ax1, np.zeros(3), "x")
setup(ax1, "one quark  =  one colour  =  one orientation")
ax1.text2D(0.5, 0.02, "one orientation (x / red)\n→ anisotropic → confined",
           transform=ax1.transAxes, ha="center", fontsize=10, color="#333")

# panel 2: baryon = 3 orthogonal colour-bipyramids meeting at the shared vertex
ax2 = fig.add_subplot(122, projection="3d")
for axis in "xyz":
    add_bipyramid(ax2, np.zeros(3), axis, alpha=0.42)
add_c3_circuit(ax2)
setup(ax2, "proton  =  red + green + blue  =  x + y + z")
ax2.text2D(0.5, 0.02, "all three orientations meet at the vertex\n→ isotropic → colour-neutral",
           transform=ax2.transAxes, ha="center", fontsize=10, color="#333")

fig.text(0.5, 0.015,
         "Oblate square bipyramids tile space with NO overlap; the three quark-cells share only the central vertex (and faces), they do not superpose.  "
         "Geometric reading — the rigorous binding (confinement) is the open Yang–Mills problem.",
         ha="center", fontsize=8.6, color="#666")

out = "/private/tmp/claude-501/-Users-davidelliman-tax-bridge/0ee20d99-29dd-449b-870b-32dbe167ea89/scratchpad/baryon3d.png"
fig.savefig(out, dpi=140, bbox_inches="tight"); print("wrote", out)
