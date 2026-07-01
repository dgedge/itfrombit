#!/usr/bin/env python3
r"""Candidate sector-billing theorem for dressed alpha: Ward identity -> current-current correlator.

The internal CPTP route reads dressed alpha off the two-link OCCUPATION (a monitored Born weight).
But measured low-energy QED alpha is billed by the CURRENT/SELF-ENERGY slot: the Ward identity
Z1 = Z2 forces the charge renormalisation to come ENTIRELY from the photon self-energy
Z3 = <J J> (the vacuum polarisation), not from any state occupation -- e_phys = sqrt(Z3) e_bare.
So the suitable sector-billing observable is the substrate's Peierls current-current correlator, not
p_link.

CANDIDATE THEOREM (Ward + optical theorem + Kramers-Kronig):
  Bill dressed alpha by the substrate vacuum polarisation Pi_sub(0):
        1/alpha_phys = (1/alpha0) * (1 - Pi_sub(0)),     delta = -137 * Pi_sub(0).
  Pi_sub is the J-resolvent on the bridge spectrum, Pi_sub(w) = g_v^2 sum_g |<E_g|J>|^2/(w - w_g),
  g_v^2 = alpha0 * G_portal^2 the established emission vertex.  Its IMAGINARY part is the web
  emission spectrum (the route's Gamma_g, optical theorem); its REAL part at w=0 -- which bills alpha
  -- is the Kramers-Kronig transform, the dispersion sum  Re Pi_sub(0) = -g_v^2 sum_g w_g-weight/w_g.
  This bills alpha by the self-energy slot exactly as QED requires, is parameter-free given the
  spectrum + the established vertex, and uses the SAME spectrum the route already computes -- but
  combined DISPERSIVELY (Re part), not as an occupation.

The script computes:
  [1] the current normalisation (f-sum) and the Ward structure of |J> (does it carry net energy /
      connect the ground state?);
  [2] Pi_sub(0) bare (matter loop) and web-dressed (KK from the emission spectrum) -> delta;
  [3] comparison to the target 0.036 and to the occupation readout (0.9956 x target);
  [4] the honest normalisation caveats that decide whether this is SUITABLE or merely right-order.
"""
from __future__ import annotations

import numpy as np

import dressed_alpha_bridge_web_open_system as bw
import bridge_web_lindblad_keldysh_poles as poles


H2, PAIRS, IDX, _BAS = bw.build_pair_system()
T_M = 1.0 / 3.0
A0, G = bw.ALPHA0, bw.G_PORTAL
DELTA_TARGET = bw.DELTA_TARGET

EVALS_UNIT, EVECS = np.linalg.eigh(H2)
EVALS = T_M * EVALS_UNIT
E_REF = float(EVALS.min())
OMEGA = EVALS - E_REF
JVEC = bw.current_portal(IDX, bw.ETA_PIN)
AMPS = EVECS.T @ JVEC                      # <E_g | J>
GROUPS = poles.eigen_groups(EVALS)
S_ETA, _OM = bw.photon_form_factor(n_grid=240)


def main() -> None:
    print("CANDIDATE SECTOR-BILLING THEOREM — Ward identity: alpha billed by <J J>, not the occupation")
    print("=" * 96)

    # [1] current normalisation + Ward structure
    total_w = float(np.sum(AMPS ** 2))
    jhj = float(JVEC @ H2 @ JVEC) * T_M
    gs_weight = float(np.sum(AMPS[OMEGA < 1e-6] ** 2))
    print("\n[1] current normalisation + Ward structure of the Peierls portal |J>")
    print(f"    sum_g |<E_g|J>|^2            = {total_w:.6f}   (=1: J normalised — the f-sum)")
    print(f"    <J|H|J>                      = {jhj:.3e}    (=0: the current carries no net energy)")
    print(f"    weight at omega~0 (ground)   = {gs_weight:.3e}    (small: J connects vacuum->excited)")

    # [2] vacuum polarisation Pi_sub(0); g_v^2 = alpha0 G^2 (the established emission vertex)
    bare = web = emis = 0.0
    for g in GROUPS:
        w = float(np.sum(AMPS[g] ** 2))
        og = float(np.mean(OMEGA[g]))
        if og < 1e-6:
            continue
        se = S_ETA(og, bw.ETA_PIN)
        bare += w / og
        web += w * se / og
        emis += w * se
    gv2 = A0 * G * G
    delta_bare = 137.0 * gv2 * bare          # delta = -137 Pi(0), Pi(0) = -g_v^2 sum w/omega
    delta_web = 137.0 * gv2 * web
    gamma_esc = 2 * np.pi * gv2 * emis       # cross-check: the route's Im-part escape ~1.57e-3
    print("\n[2] substrate vacuum polarisation  Pi_sub(0)  (Ward billing: alpha <- <J J>)")
    print(f"    bare dispersion   sum_g w/omega         = {bare:.4f}")
    print(f"    web-dressed (KK)  sum_g w*S_eta/omega   = {web:.4f}")
    print(f"    g_v^2 = alpha0*G^2                       = {gv2:.6e}")
    print(f"    -> delta (bare matter loop)             = {delta_bare:.4e}   /target {delta_bare / DELTA_TARGET:.3f}")
    print(f"    -> delta (web-dressed, KK)              = {delta_web:.4e}   /target {delta_web / DELTA_TARGET:.3f}")
    print(f"    [x-check] Gamma_esc = 2pi g_v^2 sum wS  = {gamma_esc:.4e}   (route Gamma_esc 1.57e-3)")

    # [3] comparison of billing observables
    print("\n[3] billing observables compared (all on the same substrate, x target)")
    print(f"    {'observable':<42} {'delta/target':>12}  {'slot':<22}")
    print(f"    {'occupation p_link (route)':<42} {0.996:>12.3f}  {'service Born weight':<22}")
    print(f"    {'Ward <JJ> bare matter loop':<42} {delta_bare / DELTA_TARGET:>12.3f}  {'self-energy (dispersive)':<22}")
    print(f"    {'Ward <JJ> web-dressed (KK)':<42} {delta_web / DELTA_TARGET:>12.3f}  {'self-energy (dispersive)':<22}")

    # [4] verdict
    bare_ratio = delta_bare / DELTA_TARGET
    web_ratio = delta_web / DELTA_TARGET
    print("\n" + "=" * 96)
    print("VERDICT — is the Ward current-current correlator a suitable sector-billing candidate?")
    print("  PRINCIPLED: yes. Z1=Z2 makes charge renormalisation = the photon self-energy Z3 = <JJ>,")
    print("  so billing alpha by the substrate vacuum polarisation (not the occupation) is the")
    print("  QED-correct slot. The Peierls J is the natural current; Im Pi = the web emission spectrum")
    print("  the route already has; Re Pi(0) = its Kramers-Kronig transform (the dispersion sum).")
    print(f"  MAGNITUDE: the self-energy billing gives delta/target = {bare_ratio:.2f} (bare loop) and")
    print(f"  {web_ratio:.2f} (web-dressed). It is the RIGHT ORDER and brackets/approaches the target,")
    print("  but it does NOT reproduce the occupation route's 0.996 — it is a DIFFERENT number, as it")
    print("  must be (occupation != self-energy). So a Ward sector-billing theorem RELOCATES the")
    print("  magnitude off 0.996 onto the dispersion value.")
    print("  CAVEATS (what keeps it a candidate, not a closure): (a) the overall normalisation has")
    print("  O(1) freedom (the KK subtraction point, the 2pi/2 in Im Pi=Gamma/2, the renormalisation")
    print("  scale) that the bare-vs-web spread (>2x) exposes; (b) Ward REQUIRES |J> to be the exactly")
    print("  conserved lattice current (continuity), which is the item-115 Ward-compatible-current")
    print("  question, not yet closed on this bridge. SUITABLE in form and slot; not yet a derivation")
    print("  of the magnitude until the normalisation and the conserved-current identity are pinned.")
    print(f"\n  CONSISTENCY (decisive): the web-dressed Ward self-energy delta = {delta_web:.4f} lands inside")
    print("  the independently-established charge-weighted 1-loop vacuum polarisation (DRIFT 1429:")
    print("  Sum Q_f^2 = 8-16 -> delta 0.009-0.019). So billing the PHYSICAL self-energy slot reproduces")
    print("  standard QED's 1-loop magnitude AND its 2-4x undershoot of 0.036 — it CONFIRMS the magnitude")
    print("  is not closed (per 1429, no loop order reaches 0.036), rather than rescuing the occupation's")
    print("  0.996. The occupation (service) slot and the self-energy (physical) slot are different, and")
    print("  the physical one gives the known undershoot. No billing theorem makes the physical slot")
    print("  equal 0.996; 137.036 stays the Part-12 fit.")


if __name__ == "__main__":
    main()
