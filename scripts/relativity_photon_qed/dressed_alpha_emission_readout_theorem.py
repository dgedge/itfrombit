#!/usr/bin/env python3
r"""Observable-theorem test for the dressed-alpha CPTP route: is the current-link readout the
item-79 emission probability?

THE THEOREM HALF. The CPTP route (dressed_alpha_monitor_web_cptp_audit.py) reads the dressed coupling
off the two Peierls current-link occupations, alpha_eff = (p_link/p_uniform)*alpha0. The bare leg is
already a theorem: item79_unital_channel.py proves (Evans-Frigerio) that the monitored bridge has the
UNIQUE fixed point I/137, so equipartition is derived and alpha0 = P(emission mode) = 1/137 -- a
trace-of-dynamics, not a count. But that licensed readout is the occupation of an EXPLICIT emission mode
(the "+1" 137th state); the dressed route instead reads the CONFINED current-link occupation. They
coincide at equipartition (both 1/137); the dressed (off-equipartition) equivalence is the open gap.

This script tests it decisively. Build the dressed system in the EXPLICIT-emission-mode picture, where
alpha = P(emission mode) is theorem-licensed by item 79:
    D = 137  (136 confined pairs + 1 emission mode),
    H couples the emission mode to the confined sector through the Peierls current portal (strength w),
    monitor = position dephasing on the 136 confined states (the item-79 syndrome extraction),
    web loss = the emission mode leaks to a sink at rate kappa (the CPTP-to-sink escape).
On the dressed quasi-stationary mode read BOTH observables on the SAME state:
    A (licensed):  alpha = P(emission mode) = rho[em,em]          -> delta_A = 1/alpha - 137
    B (the route): alpha = p_link / p_uniform * alpha0            -> delta_B = 1/alpha - 137
If delta_B == delta_A across (w, kappa, rate), the route's confined current-link readout IS the item-79
emission probability (readout LICENSED). If they diverge, the current-link readout is a separate choice.
"""
from __future__ import annotations

import numpy as np

import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_monitor_web_cptp_audit as aud


H2, PAIRS, IDX, _BAS = bw.build_pair_system()
DCONF = len(IDX)                       # 136 confined pair states
I08, I19 = IDX[(0, 8)], IDX[(1, 9)]
T_M = 1.0 / 3.0
A0 = bw.ALPHA0                         # 1/137


def dressed_explicit(w: float, kappa: float, rate: float, eta: int = bw.ETA_PIN):
    """Dressed QSS in the explicit 137-dim emission-mode picture; return both readouts."""
    d = DCONF
    D = d + 1
    em = d
    H = np.zeros((D, D))
    H[:d, :d] = T_M * H2
    port = bw.current_portal(IDX, eta)          # Peierls current portal on (0,8),(1,9)
    H[:d, em] = w * port
    H[em, :d] = w * port
    pconf = np.ones(D); pconf[em] = 0.0         # diagonal of the confined projector

    def matvec(x: np.ndarray) -> np.ndarray:
        rho = x.reshape((D, D))
        out = -1j * (H @ rho - rho @ H)
        # position dephasing on the confined sector: rate*( diag_conf(rho) - 1/2 {P_conf, rho} )
        dm = -0.5 * (pconf[:, None] * rho + rho * pconf[None, :])
        dm[np.arange(d), np.arange(d)] += np.diag(rho)[:d]
        out += rate * dm
        # web loss: emission mode leaks to a sink, -1/2 kappa {|em><em|, rho}
        loss = np.zeros_like(rho)
        loss[em, :] += rho[em, :]
        loss[:, em] += rho[:, em]
        out += -0.5 * kappa * loss
        return out.reshape(-1)

    ev, vec, res, _solver = aud.scipy_dominant(matvec, D * D)
    rho = vec.reshape((D, D))
    rho = 0.5 * (rho + rho.conj().T)
    tr = np.trace(rho).real
    if tr < 0:
        rho = -rho; tr = -tr
    rho /= tr
    diag = np.diag(rho).real
    p_em = diag[em]
    p_link = 0.5 * (diag[I08] + diag[I19])
    delta_A = 1.0 / p_em - 137.0                          # licensed: alpha = P(emission mode)
    alpha_B = (p_link / (1.0 / D)) * A0                   # route: confined current-link occ / uniform
    delta_B = 1.0 / alpha_B - 137.0
    return delta_A, delta_B, p_em, p_link, res


def readout_subchoice_within_route():
    """Caveat-free corroboration: on the route's OWN dressed QSS, how much does delta move
    between defensible confined readouts (2 current links / 4 interface pairs / full K2 cell)?"""
    h = T_M * H2
    s_eta, _ = bw.photon_form_factor(n_grid=200)
    gamma_web, *_ = aud.grouped_web_width_operator(H2, IDX, s_eta)
    labels = aud.pair_basis_labels(PAIRS, IDX)
    cminus, cplus = bw.current_portal(IDX, -1), bw.current_portal(IDX, +1)
    cand = aud.MonitorCandidate("native", "pair_site", 1.0, "", "")
    mv = aud.make_matvec(h, gamma_web, cand, labels, cminus, cplus)
    ev, vec, _res, _ = aud.scipy_dominant(mv, DCONF * DCONF)
    rho = vec.reshape((DCONF, DCONF)); rho = 0.5 * (rho + rho.conj().T)
    tr = np.trace(rho).real
    if tr < 0:
        rho = -rho; tr = -tr
    rho /= tr
    diag = np.diag(rho).real

    def delta_of(indices):
        p = float(np.mean([diag[i] for i in indices]))
        return 1.0 / ((p / (1.0 / DCONF)) * A0) - 137.0

    readouts = {
        "2 current links (0,8),(1,9)": [IDX[(0, 8)], IDX[(1, 9)]],
        "4 interface pairs": [IDX[(0, 8)], IDX[(1, 9)], IDX[(0, 9)], IDX[(1, 8)]],
        "full K2 cell {0,1,8,9}": [IDX[(min(a, b), max(a, b))] for a in (0, 1, 8, 9) for b in (0, 1, 8, 9) if a <= b],
    }
    print("\n[within-route] same dressed QSS (native tick), different confined readouts:")
    print(f"    {'readout':<30} {'delta':>11} {'/target':>9}")
    ds = {}
    for name, ids in readouts.items():
        dd = delta_of(ids)
        ds[name] = dd
        print(f"    {name:<30} {dd:>11.4e} {dd / bw.DELTA_TARGET:>9.3f}")
    vals = list(ds.values())
    print(f"    spread across these 3 defensible readouts: {max(vals) / min(vals):.2f}x in delta")
    return ds


def main() -> None:
    print("OBSERVABLE-THEOREM TEST — is the current-link readout the item-79 emission probability?")
    print("=" * 96)
    print("Both readouts on the SAME dressed QSS (explicit 137-dim emission-mode picture).")
    print("delta_A = licensed (alpha = P(emission mode), item79_unital_channel theorem).")
    print("delta_B = the CPTP route's confined current-link readout.")
    print(f"\n{'w':>5} {'kappa':>8} {'rate':>6} {'delta_A':>11} {'delta_B':>11} {'B/A':>8} {'P(em)':>9} {'p_link':>9}")
    rows = []
    for w in (0.3, 0.5, 1.0):
        for kappa in (1e-3, 3e-3, 1e-2):
            for rate in (1.0,):
                dA, dB, pem, pl, res = dressed_explicit(w, kappa, rate)
                ratio = dB / dA if dA else float("nan")
                rows.append((w, kappa, rate, dA, dB, ratio))
                print(f"{w:>5.1f} {kappa:>8.0e} {rate:>6.2f} {dA:>11.4e} {dB:>11.4e} {ratio:>8.3f} {pem:>9.5f} {pl:>9.5f}")

    ratios = [r[5] for r in rows if np.isfinite(r[5])]
    rmean = float(np.mean(ratios)); rspread = float(np.std(ratios))
    rmin, rmax = min(ratios), max(ratios)
    print("\n" + "=" * 96)
    print("VERDICT (observable theorem):")
    print(f"  delta_B / delta_A = {rmean:.3f} +/- {rspread:.3f}  (range {rmin:.3f}..{rmax:.3f}) across w, kappa.")
    licensed_eq = abs(rmean - 1.0) < 0.05 and rspread < 0.03
    licensed_const = rspread / rmean < 0.05
    if licensed_eq:
        print("  => LICENSED (equal): the confined current-link occupation reproduces the item-79")
        print("     emission probability exactly. The route's readout IS the licensed observable.")
    elif licensed_const:
        print(f"  => LICENSED UP TO A FIXED FACTOR {rmean:.3f}: the current-link readout tracks the")
        print("     emission probability with a constant proportionality (a fixed, derivable rescaling),")
        print("     NOT readout-dominated -- but the route's bare normalisation absorbs this factor.")
    else:
        print("  => NOT LICENSED: delta_B/delta_A varies with the couplings, so the confined")
        print("     current-link occupation is a SEPARATE readout choice, not the item-79 emission")
        print("     probability -- the readout joins t_m/t_web, eta, E_ref, rate as an assignment knob.")
    print("\n  NB scope: this licenses the readout PRINCIPLE (confined current-link occ vs emission")
    print("  probability) in the explicit-emission-mode model; the bare leg (equipartition ->")
    print("  alpha0=1/137) is already proven (item79_unital_channel.py, Evans-Frigerio).")

    readout_subchoice_within_route()


if __name__ == "__main__":
    main()
