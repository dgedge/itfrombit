#!/usr/bin/env python3
r"""0nubb 'lobster' plot: m_bb vs the lightest neutrino mass, with the framework's
prediction marked. Saves item87_mbb_lobster.{png,pdf}.

Bands: standard NuFIT-central NO/IO allowed regions (Majorana phases scanned).
Framework point: the pinned neutrino sector (Koide R_nu=1, delta_nu=1/3 -> spectrum
m=(0.79,8.72,50.2) meV, NO) with the discrete sign-pointer Majorana parities ->
m_bb ~ 1.8 meV (sigma=+) / 3.1 meV (sigma=-); envelope [0.7,4.2] meV.
"""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# --- oscillation parameters (NuFIT-central) ---
S12sq, S13sq = 0.303, 0.0222
DM21 = 7.42e-5
DM31_NO = 2.515e-3
DM32_IO = 2.498e-3
C12sq, C13sq = 1 - S12sq, 1 - S13sq
Ue = np.array([C12sq * C13sq, S12sq * C13sq, S13sq])   # |U_ei|^2


def band(mlmin, mlmax, ordering, n=400, nph=61):
    ml = np.logspace(math.log10(mlmin), math.log10(mlmax), n)
    lo, hi = np.zeros(n), np.zeros(n)
    ph = np.linspace(0, math.pi, nph)                  # only relative phases matter for |.|
    a2, a3 = np.meshgrid(ph, ph)
    e2, e3 = np.exp(2j * a2).ravel(), np.exp(2j * a3).ravel()
    for i, m0 in enumerate(ml):
        if ordering == "NO":
            m1, m2, m3 = m0, math.sqrt(m0**2 + DM21), math.sqrt(m0**2 + DM31_NO)
        else:
            m3 = m0; m1 = math.sqrt(m0**2 + DM32_IO); m2 = math.sqrt(m1**2 + DM21)
        vals = np.abs(Ue[0]*m1 + Ue[1]*m2*e2 + Ue[2]*m3*e3)
        lo[i], hi[i] = vals.min(), vals.max()
    return ml, lo, hi


def main():
    fig, ax = plt.subplots(figsize=(7.2, 5.6))

    # NO / IO allowed bands (meV)
    for ordering, color, label in [("IO", "#3b6fb0", "Inverted ordering"),
                                   ("NO", "#3fa34d", "Normal ordering")]:
        ml, lo, hi = band(1e-4, 0.3, ordering)
        ax.fill_between(ml * 1e3, lo * 1e3, hi * 1e3, color=color, alpha=0.28, lw=0)
        ax.plot(ml * 1e3, lo * 1e3, color=color, lw=0.8, alpha=0.7)
        ax.plot(ml * 1e3, hi * 1e3, color=color, lw=0.8, alpha=0.7)

    # experimental sensitivity bands (NME-spread)
    ax.axhspan(28, 122, color="0.5", alpha=0.16, lw=0)
    ax.axhline(28, color="0.4", lw=0.8, ls="--")
    ax.text(1.2e-1, 130, "KamLAND-Zen (2024) exclusion  (NME spread)", fontsize=7.5, color="0.35", va="bottom")
    ax.axhspan(9, 20, color="#b06a3b", alpha=0.14, lw=0)
    ax.text(1.2e-1, 20.5, "next-gen reach: LEGEND-1000 / nEXO (~9-20 meV)", fontsize=7.5, color="#b06a3b", va="bottom")

    # framework prediction: spectrum m1 = 0.79 meV, m_bb = 1.8 (sigma+) and 3.1 (sigma-) meV
    m1_fw = 0.79
    ax.plot([m1_fw, m1_fw], [0.7, 4.2], color="#c0392b", lw=2.4, solid_capstyle="round",
            alpha=0.55, zorder=5)                                   # continuous-phase envelope
    ax.plot(m1_fw, 1.8, marker="*", ms=17, color="#c0392b", mec="k", mew=0.6, zorder=6)
    ax.plot(m1_fw, 3.1, marker="*", ms=13, color="#e08e0b", mec="k", mew=0.6, zorder=6)
    ax.annotate(r"framework:  $m_{\beta\beta}\approx1.8$-$3.1$ meV",
                (m1_fw, 3.1), (1.6, 7.0), fontsize=9, color="#c0392b",
                arrowprops=dict(arrowstyle="-", color="#c0392b", lw=0.8))
    ax.annotate(r"$m_1\!\approx\!0.79$ meV" "\n" r"($\Sigma m_\nu\!\approx\!60$ meV, NO)",
                (m1_fw, 0.7), (0.13, 0.16), fontsize=8, color="#7a2519",
                arrowprops=dict(arrowstyle="-", color="#7a2519", lw=0.8))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e-1, 3e2); ax.set_ylim(0.3, 5e2)
    ax.set_xlabel(r"lightest neutrino mass  $m_{\rm lightest}$  [meV]")
    ax.set_ylabel(r"effective Majorana mass  $m_{\beta\beta}$  [meV]")
    ax.set_title(r"Neutrinoless double-beta decay: framework prediction ($m_{\beta\beta}\!\sim\!2$-$3$ meV, NO)",
                 fontsize=10.5)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

    handles = [
        Patch(facecolor="#3fa34d", alpha=0.35, label="Normal ordering (NuFIT)"),
        Patch(facecolor="#3b6fb0", alpha=0.35, label="Inverted ordering (NuFIT)"),
        Line2D([], [], color="#c0392b", marker="*", ms=13, mec="k", mew=0.5, ls="none",
               label=r"framework $\sigma=+$ ($m_{\beta\beta}\approx1.8$ meV)"),
        Line2D([], [], color="#e08e0b", marker="*", ms=11, mec="k", mew=0.5, ls="none",
               label=r"framework $\sigma=-$ ($\approx3.1$ meV)"),
        Line2D([], [], color="#c0392b", lw=2.4, alpha=0.55,
               label=r"framework phase envelope [0.7, 4.2] meV"),
    ]
    ax.legend(handles=handles, fontsize=8, loc="lower right", framealpha=0.92)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"python_code/item87_mbb_lobster.{ext}", dpi=160, bbox_inches="tight")
    print("saved python_code/item87_mbb_lobster.png and .pdf")


if __name__ == "__main__":
    main()
