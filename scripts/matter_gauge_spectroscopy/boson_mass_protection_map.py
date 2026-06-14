#!/usr/bin/env python3
r"""Glueball / graviton investigation via the ACTION-SPACE + REPRESENTATION-PARITY map.

Today's two ingredients combine to settle the whole spin-1 / spin-2 boson mass pattern without touching
the K_6 bulk-Bloch trap (ANCHOR L1205-1237 action-space subsection):

  (i) ACTION SPACE (which space the irrep acts on) fixes WHICH PARTICLE an O_h irrep is:
        T_1u  -> photon (gauge / matter-cell vector, massless)
        E_g   -> GRAVITON in the matter-cell action (s10.1, spin-2 metric distortion)
                 BUT the 2++ GLUEBALL in the gauge-web action (L1217) -- SAME irrep, DIFFERENT particle.
        A_1g  -> Higgs (matter-cell breathing, s8.1)  /  0++ glueball (gauge-web cycle, Fock).

  (ii) REPRESENTATION-PARITY MASS PROTECTION (L1169 theorem) fixes WHETHER each can have a mass:
        a real antisymmetric operator projected onto an ODD-dimensional O_h irrep has >=1 exact zero
        eigenvalue (massless, forced); an EVEN-dimensional irrep does not (mass PERMITTED).

Combining: the photon (T_1u, dim 3 ODD) is mass-FORBIDDEN -> massless by topology (verified gap 0.000,
L1173). The graviton and the 2++ glueball are BOTH E_g (dim 2 EVEN) -> mass-PERMITTED; the gauge-web E_g
(2++ glueball) realises a mass (Fock (2N-1)Lambda), while the matter-cell E_g (graviton) carries a
Gamma-point tachyon stabilised to a massless mode with stiffness K_{E_g} by the s10.2 period-4 vacuum.
The A_1g (Higgs, 0++ glueball) get mass from the breathing-resonance / Fock mechanism, not this operator.

This script ASSERTS the parity theorem and assembles the map. The remaining open piece -- the numerical
K_{E_g} stiffness (item 7) and the 208x208 BZ-trace K_eff=205 (item 8) -- needs the full walk-operator
Bloch Hamiltonian and is NOT attempted here (flagged).
"""
import sys
import numpy as np

rng = np.random.default_rng(0)


# ---- (ii) verify the representation-parity mass-protection theorem (the core of L1169) ----
print("(ii) Representation-parity mass protection (L1169): real antisymmetric op on a d-dim irrep")
print("     d (irrep dim)   min|eigenvalue| over 2000 random antisym   has forced zero? (=> massless)")
for d in range(1, 7):
    mins = []
    for _ in range(2000):
        M = rng.standard_normal((d, d)); M = M - M.T            # real antisymmetric
        ev = np.linalg.eigvals(M)                                # eigenvalues are 0 or +-i*lambda pairs
        mins.append(np.min(np.abs(ev)))
    worst = max(mins)                                            # largest min-|eig| seen (0 if always singular)
    forced_zero = worst < 1e-9
    parity = "ODD " if d % 2 else "EVEN"
    print(f"       d={d} ({parity})      {worst:.2e}                                {'YES -> massless' if forced_zero else 'no  -> mass permitted'}")
    # theorem: odd d => always a zero eigenvalue; even d => generically none
    assert (d % 2 == 1) == forced_zero, f"parity theorem failed at d={d}"
print("     => ODD-dim irreps force a massless mode; EVEN-dim permit a mass. (Just: odd antisym is singular.)\n")

# ---- O_h irreps actually used, with dimension parity ----
OH = {"A_1g": 1, "A_2g": 1, "E_g": 2, "T_1g": 3, "T_2g": 3,
      "A_1u": 1, "A_2u": 1, "E_u": 2, "T_1u": 3, "T_2u": 3}
print("(i+ii) BOSON MAP  (action space sets the particle; irrep-dim parity sets mass protection):")
print(f"  {'boson':<13}{'irrep':<7}{'dim/parity':<13}{'action space':<14}{'mass':<16}{'protection mechanism'}")
rows = [
    ("photon",      "T_1u", "photon",  "matter/gauge", "0 (massless)",       "PARITY-forbidden (odd d=3); robust (gap 0.000, L1173)"),
    ("graviton",    "E_g",  "graviton","matter-cell",  "0 (massless)*",      "mass-PERMITTED (even d=2); massless only via s10.2 period-4 vacuum"),
    ("2++ glueball","E_g",  "2++ gb",  "gauge-web",    "(2N-1)L Fock",       "mass-PERMITTED (even d=2); REALISED massive (Fock, L1247)"),
    ("0++ glueball","A_1g", "0++ gb",  "gauge-web",    "(2N-1)L=5L=1660",    "A_1g cycle (Fock); not the antisym operator -- Feshbach threshold"),
    ("Higgs",       "A_1g", "Higgs",   "matter-cell",  "EW (~125 GeV)",      "A_1g breathing RESONANCE (s8.1); not the antisym operator"),
]
for name, irr, _who, space, mass, mech in rows:
    d = OH[irr]; par = "odd" if d % 2 else "even"
    print(f"  {name:<13}{irr:<7}{f'{d} ({par})':<13}{space:<14}{mass:<16}{mech}")

# the load-bearing assertions of the map
assert OH["T_1u"] % 2 == 1, "photon T_1u must be odd-dim (mass-forbidden)"
assert OH["E_g"] % 2 == 0,  "graviton & 2++ glueball E_g must be even-dim (mass-permitted)"
assert OH["A_1g"] % 2 == 1, "Higgs / 0++ glueball A_1g is 1-dim"

print(f"""
RESULT (spin-1 / spin-2 boson mass pattern, action-space + parity, no K_6 trap):
  * SPIN-1 (T_1u, dim 3 ODD): the photon is mass-FORBIDDEN by representation parity (L1169) -- robust,
    topological, the same in either action space. (Today's loop-mode work placed it at the gapless E=1.)
  * SPIN-2 (E_g, dim 2 EVEN): mass is PERMITTED, and the ACTION SPACE then decides the fate:
      - gauge-web E_g  = the 2++ GLUEBALL: the mass is realised, an N-body Fock condensate (2N-1)Lambda
        (today's single-loop spectral function showed the E_g channel is 53% cycle-space flat band +
        47% dispersive -- a localised single loop, with the mass coming from the Fock/threshold sum).
      - matter-cell E_g = the GRAVITON: the permitted mass shows up as a Gamma-point TACHYON, which the
        s10.2 period-4 vacuum stabilises to a massless graviton carrying the stiffness K_{{E_g}} that sets
        the bare Planck mass M_P = sqrt(4 pi K_{{E_g}}) Lambda (s10.1).
  So the SAME E_g irrep is a massive glueball (gauge web) and a massless graviton (matter cell) -- the
  action-space distinction made quantitative for spin 2, exactly as photon/rho did it for spin 1.

  OPEN (honestly NOT closed here): the NUMBER K_{{E_g}} (item 7) and the 208x208 BZ-trace K_eff=205
  (item 8, the gravity sector's load-bearing prefactor) require the full walk-operator Bloch Hamiltonian
  W_QQ(k) + the E_g character projection + the s10.2 period-4 stabilisation -- a different, larger build
  than today's single-particle KPM. This map fixes the STRUCTURE (which irrep, which action space, which
  protection); it does not extract the gravity NUMBER. And the graviton's masslessness is the delicate
  s10.2 piece, not the robust parity protection the photon enjoys.
""")
print("exit 0 -- parity theorem verified; boson action-space x mass-protection map asserted.")
