#!/usr/bin/env python3
"""Generate a TWO-panel solid/translucent-face TikZ of the baryon geometry, with correct
painter's-algorithm occlusion (faces emitted back-to-front for a fixed orthographic view)."""
import numpy as np

elev, azim = np.radians(22), np.radians(35)
d = np.array([np.cos(elev) * np.cos(azim), np.cos(elev) * np.sin(azim), np.sin(elev)])  # toward camera
u = np.cross([0, 0, 1.], d); u /= np.linalg.norm(u)
v = np.cross(d, u); v /= np.linalg.norm(v)
SC = 3.3
AX = {"x": np.array([1, 0, 0.]), "y": np.array([0, 1, 0.]), "z": np.array([0, 0, 1.])}
PERP = {"x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y")}
COL = {"x": "qr", "y": "qg", "z": "qb"}


def pr(p):  # -> "(sx,sy)" string  and depth
    p = np.asarray(p, float); s = SC * np.array([p @ u, p @ v]); return f"({s[0]:.3f},{s[1]:.3f})", p @ d


def s2(p):
    p = np.asarray(p, float); return SC * np.array([p @ u, p @ v])


def bipyramid_faces(axis, a=1.0, O=np.zeros(3)):
    dd = AX[axis]; e1, e2 = AX[PERP[axis][0]], AX[PERP[axis][1]]
    m = O + 0.5 * a * dd
    w = [m + 0.5 * a * (su * e1 + sv * e2) for su, sv in [(1, 1), (1, -1), (-1, -1), (-1, 1)]]
    Ai, Ao = O, O + a * dd
    F = [[Ai, w[i], w[(i + 1) % 4]] for i in range(4)] + [[Ao, w[i], w[(i + 1) % 4]] for i in range(4)]
    return F


def faces_block(axes):
    items = []
    for ax in axes:
        for f in bipyramid_faces(ax):
            depth = np.mean([p @ d for p in f])
            items.append((depth, COL[ax], f))
    items.sort(key=lambda t: t[0])  # far -> near
    out = []
    for _, c, f in items:
        pts = "--".join(pr(p)[0] for p in f)
        out.append(f"  \\fill[{c},fill opacity=0.50,draw=black!70,line width=0.35pt,line join=round] {pts}--cycle;")
    return "\n".join(out)


def cube_block():
    r = 0.5; cor = {(sx, sy, sz): np.array([sx * r, sy * r, sz * r]) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)}
    E = [((-1, -1, -1), (1, -1, -1)), ((-1, -1, -1), (-1, 1, -1)), ((-1, -1, -1), (-1, -1, 1)),
         ((1, 1, -1), (-1, 1, -1)), ((1, 1, -1), (1, -1, -1)), ((1, 1, -1), (1, 1, 1)),
         ((-1, 1, 1), (1, 1, 1)), ((-1, 1, 1), (-1, -1, 1)), ((-1, 1, 1), (-1, 1, -1)),
         ((1, -1, 1), (-1, -1, 1)), ((1, -1, 1), (1, 1, 1)), ((1, -1, 1), (1, -1, -1))]
    return "\n".join(f"  \\draw[gray!50,dashed,very thin] {pr(cor[a])[0]}--{pr(cor[b])[0]};" for a, b in E)


def axes_block():
    L = [("x", "qr", "$x$ (r)", "north"), ("y", "qg", "$y$ (g)", "west"), ("z", "qb", "$z$ (b)", "west")]
    out = []
    for ax, c, lab, anc in L:
        tip = 1.42 * AX[ax]
        out.append(f"  \\draw[{c},-{{Stealth[length=4pt]}},thick] {pr(np.zeros(3))[0]}--{pr(tip)[0]} node[anchor={anc},{c},font=\\small]{{{lab}}};")
    out.append(f"  \\shade[ball color=black] {pr(np.zeros(3))[0]} circle (0.8pt);")
    return "\n".join(out)


def c3_block():
    n = np.array([1, 1, 1.]) / np.sqrt(3)
    out = [f"  \\draw[gray!70,dashed,thin] {pr(-0.55*n)[0]}--{pr(1.12*n)[0]};"]
    e1 = np.array([1, -1, 0.]) / np.sqrt(2); e2 = np.cross(n, e1); e2 /= np.linalg.norm(e2)
    cen = 0.42 * n; R = 0.66
    ang = lambda dv: np.arctan2((dv - (dv @ n) * n) @ e2, (dv - (dv @ n) * n) @ e1)
    angs = {k: ang(AX[k]) for k in "xyz"}
    for a, b in [("x", "y"), ("y", "z"), ("z", "x")]:
        a0 = angs[a]; dd = (angs[b] - a0 + np.pi) % (2 * np.pi) - np.pi
        aa = np.linspace(a0 + 0.16 * dd, a0 + 0.82 * dd, 16)
        poly = "--".join(pr(cen + R * (np.cos(t) * e1 + np.sin(t) * e2))[0] for t in aa)
        out.append(f"  \\draw[{COL[a]},line width=1.4pt,-{{Stealth[length=5pt,bend]}}] {poly};")
    return "\n".join(out)


def panel(axes, title, caption, with_c3):
    pieces = [cube_block()]
    if with_c3:
        n = np.array([1, 1, 1.]) / np.sqrt(3)
        pieces.append(f"  \\draw[gray!70,dashed,thin] {pr(-0.55*n)[0]}--{pr(1.12*n)[0]};")
    pieces.append(faces_block(axes))
    pieces.append(axes_block())
    if with_c3:
        pieces.append(c3_block())
    # bbox for title/caption placement
    ext = [s2(c * AX[k]) for k in "xyz" for c in (1.42, -0.55)] + [s2(np.array([sx * .5, sy * .5, sz * .5])) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    ext = np.array(ext); cx = (ext[:, 0].min() + ext[:, 0].max()) / 2
    ytop, ybot = ext[:, 1].max(), ext[:, 1].min()
    pieces.append(f"  \\node[font=\\bfseries,align=center] at ({cx:.2f},{ytop+0.7:.2f}) {{{title}}};")
    pieces.append(f"  \\node[align=center,gray!55!black,font=\\small] at ({cx:.2f},{ybot-0.6:.2f}) {{{caption}}};")
    return "\\begin{tikzpicture}[line cap=round]\n" + "\n".join(pieces) + "\n\\end{tikzpicture}"


doc = r"""\documentclass[border=12pt,varwidth=27cm]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta,bending}
\definecolor{qr}{HTML}{D23B3B}\definecolor{qg}{HTML}{2E9E44}\definecolor{qb}{HTML}{3B6FD2}
\begin{document}
\centering
{\large\bfseries Colour $=$ orientation: how three quarks build a proton}\par\smallskip
{\color{gray!60!black}\small geometric reading of the framework (solid / translucent faces)}\par\medskip
""" + panel(["x"], "one quark $=$ one colour $=$ one orientation",
            r"one orientation (x / red)\\ $\to$ anisotropic $\to$ \textbf{confined}", False) + \
    "\\hspace{1.0cm}%\n" + \
    panel(["x", "y", "z"], "proton $=$ r $+$ g $+$ b $=$ x $+$ y $+$ z",
          r"three orientations meet at the vertex (C$_3$: r$\to$g$\to$b)\\ $\to$ isotropic $\to$ \textbf{colour-neutral}", True) + \
    "\n\\end{document}\n"

open("/private/tmp/claude-501/-Users-davidelliman-tax-bridge/0ee20d99-29dd-449b-870b-32dbe167ea89/scratchpad/baryon_solid.tex", "w").write(doc)
print("wrote baryon_solid.tex")
