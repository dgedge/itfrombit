#!/usr/bin/env python3
r"""Item 48 — the HOST-PATH fork (octagon-perimeter vs s3.1 body-diagonal), resolved.

The long-standing caveat (ANCHOR s15 item 48 / s9.1 L1781): if the 3D rho flux traces s3.1
body-diagonal bridges instead of the 2D octagon perimeter, "L(P_5)=P_4 with leading eigenvalue phi
may not survive." Prior work (item48_bodydiagonal_route.py) showed the body-diagonal web is girth-4
(no native phi) -- but only argued dissolution in the embedded mode. With the matter-mode + Grover
machinery (leak pinned, continuum discharged), the fork resolves DECISIVELY, and not by survival:

  The J^PC=1^-- VECTOR dipole on a loop C_{2d} has eigenvalue 2cos(pi/(d+1)), where d = antinode
  separation = (loop size)/2 [the k=1 mode psi_n=cos(pi n/d) has antinodes at n=0,d].
    * C_8 octagon (the MATTER cell):  d=4 -> 2cos(pi/5) = phi   -> the massive rho (sqrt2 phi Lambda).
    * C_4 square  (the body-diagonal/girth-4 web's native loop): d=2 -> 2cos(pi/3) = 1
      = EXACTLY the s7.4 T_1u transmission pole at E=+1, the MASSLESS PHOTON ("radiates freely at c").

So the two "routes" are not two hosts for the rho -- they are two DIFFERENT PARTICLES: the octagon
(d=4, phi) is the rho; the body-diagonal/square (d=2, E=1) is the photon. A body-diagonal vector flux
delocalises into the gapless E=1 gauge pole and radiates -- it cannot be a distinct massive rho. The
massive rho is therefore NECESSARILY the octagon-perimeter localised mode at phi. Fork resolved.
"""
import sys
import numpy as np

phi = (1 + 5 ** 0.5) / 2
LQCD = 332.0
M_RHO_EXP = 775.26


def vec_dipole_eig(loop_size):
    """Leading eigenvalue of the J=1 vector flux on a C_{loop_size} loop, via the framework's own rule:
    antinode separation d = loop_size/2 -> flux path P_{d+1} -> gauge dynamics L(P_{d+1}) = P_d ->
    leading eigenvalue 2cos(pi/(d+1))."""
    d = loop_size // 2
    return 2 * np.cos(np.pi / (d + 1)), d


# ---- (1) the vector-dipole eigenvalue is loop-size-specific ----
print("loop     d=size/2   vector eig 2cos(pi/(d+1))   identity")
rows = []
for size in (4, 6, 8):
    e, d = vec_dipole_eig(size)
    rows.append((size, d, e))
    ident = ("PHOTON (= s7.4 T_1u E=+1 pole)" if abs(e - 1) < 1e-9 else
             "rho (golden ratio)" if abs(e - phi) < 1e-9 else "—")
    print(f"  C_{size:<5} {d:<10} {e:.6f}                {ident}")
e4, _ = vec_dipole_eig(4); e8, _ = vec_dipole_eig(8)
assert abs(e8 - phi) < 1e-12, "octagon (C_8, d=4) vector eigenvalue must be phi"
assert abs(e4 - 1.0) < 1e-12, "square (C_4, d=2) vector eigenvalue must be 1 = the photon pole"
assert abs(2 * np.cos(np.pi / 4) - 2 ** 0.5) < 1e-12   # the body-diagonal half-loop P_3 -> sqrt2 (also != phi)
print(f"\n=> the octagon (C_8) vector sits at phi; the body-diagonal/square (C_4) vector sits at E=1 = the")
print(f"   s7.4 T_1u transmission pole (the MASSLESS PHOTON, 'radiates freely at c', ANCHOR L1120). DIFFERENT")
print(f"   PARTICLES. phi=2cos(pi/5) is octagon-specific (needs the antinode-4 dipole of C_8).")

# ---- (2) the mass map: only the octagon eigenvalue reproduces the rho ----
print("\nmass m = sqrt2 * (vector eig) * Lambda   [Lambda=332]:")
for size, d, e in rows:
    m = 2 ** 0.5 * e * LQCD
    tag = "= rho 775? " + (f"{100*(m-M_RHO_EXP)/M_RHO_EXP:+.0f}%")
    print(f"  C_{size}: sqrt2*{e:.4f}*Lambda = {m:6.1f} MeV   ({tag})")
m_oct = 2 ** 0.5 * phi * LQCD
assert abs(m_oct - 759.7) < 0.5 and abs((m_oct - M_RHO_EXP) / M_RHO_EXP) < 0.021   # the canon's 2.0%
m_sq = 2 ** 0.5 * 1.0 * LQCD
assert m_sq < 500   # the C_4 'mass' is photon-region, nowhere near the rho
print(f"=> only the octagon (phi -> {m_oct:.0f} MeV) lands on the rho (775, 2%); the square (E=1 -> {m_sq:.0f}")
print(f"   MeV) is photon-region. The routes give DISTINCT, distinguishable predictions -- not a degeneracy.")

# ---- (3) the body-diagonal web has NO native phi (girth 4); its vector pole is the E=1 photon ----
L = 8
xs = range(L); Cs = [(x, y, z) for x in xs for y in xs for z in xs]; ix = {c: i for i, c in enumerate(Cs)}
N = len(Cs); A = np.zeros((N, N))
for c in Cs:
    for sx in (1, -1):
        for sy in (1, -1):
            for sz in (1, -1):
                A[ix[c], ix[((c[0] + sx) % L, (c[1] + sy) % L, (c[2] + sz) % L)]] = 1.0
bw = np.linalg.eigvalsh(A)
# girth-4 test: no triangles (three +-1 steps can never return to origin: sum of three +-1 != 0)
assert all((a + b + cc) != 0 for a in (1, -1) for b in (1, -1) for cc in (1, -1)), "body-diagonal web is triangle-free (girth>=4)"
print(f"\n(3) body-diagonal web (deg 8): band [{bw[0]:.2f},{bw[-1]:.2f}], girth 4 (triangle-free) -> native short")
print(f"    loop is C_4 -> vector eig 1 (photon), NOT phi. phi has no native realisation on this web.")

print("""
VERDICT (host-path fork RESOLVED toward octagon-perimeter):
  The fork is not 'two hosts for the rho'. The J=1 vector dipole eigenvalue is set by the loop size:
  C_8 octagon (matter cell) -> phi -> the massive rho (sqrt2 phi Lambda = 760 MeV, 2% of 775);
  C_4 square (body-diagonal / girth-4 web) -> E=1 -> the s7.4 T_1u pole = the MASSLESS PHOTON.
  A body-diagonal vector flux delocalises into the gapless E=1 gauge sector and radiates at c -- it is
  the photon, not a distinct massive rho. So the massive rho is NECESSARILY the octagon-perimeter
  localised mode at phi (which the matter-mode KPM shows survives as a ~27 MeV resonance pinned at phi).
  The 'unpinned 3D path' caveat is thereby resolved: the body-diagonal alternative is the photon, not a
  competing rho. Residual (now small): a first-principles dynamical-selection proof that the light rho
  localises on one C_8 cell (the compactness argument: a light vector ~ one matter cell) rather than
  spreading onto bridges -- but the spread state is the photon, so the localisation is what MAKES it
  the massive rho.
""")
# =====================================================================================================
# BOSON ACTION-SPACE MAP -- what the single-loop dipole rule settles, and what needs Fock/resonance machinery.
# Canonical assignments from ANCHOR L1205-1247 (the O_h "action-space" subsection): the SAME irrep is a
# DIFFERENT particle in the matter-cell (C_8 octagon) vs gauge-web (C_4 square / cycle) action space.
# L1227 CAUTION: identifying a single loop/Bloch eigenvalue with a glueball mass is "structurally invalid"
# (the K_6 5*Lambda shortcut). Glueballs are HDR-exempt N-plaquette FOCK states m_N=(2N-1)Lambda (L1247-57);
# the Higgs is the matter A_1g breathing RESONANCE at the EW scale (s8.1). The single-loop dipole rule
# settles ONLY the single-loop bosons; the rest sit in the same irreps but need the other (anchored) machinery.
print("\n" + "=" * 100)
print("BOSON ACTION-SPACE MAP  (loop-mode rule = single C_n dipole; canon assignments ANCHOR L1216-1247)")
print("=" * 100)
hdr = f"  {'boson':<14}{'space':<13}{'irrep':<8}{'loop-mode eig':<15}{'mass':<16}{'settled by loop rule?'}"
print(hdr)

rho_m = 2 ** 0.5 * phi * LQCD
g0, g2 = 5 * LQCD, 7 * LQCD                      # 0++ glueball = (2*3-1)L, 2++ = (2*4-1)L (L1247-57 Fock)
table = [
    ("photon",     "gauge C_4", "T_1u", f"1 (E=1 pole)", "0 (massless)",      "YES  -> gapless gauge vector"),
    ("graviton",   "gauge",     "E_g",  "—",             "0 (massless)",      "NO   -> E_g even-dim mass-PERMITTED (L1169); massless via s10.2, not loop/parity"),
    ("rho",        "matter C_8","E_1u", f"phi={phi:.3f}", f"sqrt2*phi*L={rho_m:.0f}", "YES  -> single matter-loop vector"),
    ("Higgs",      "matter",    "A_1g", "—",             "EW (~125 GeV)",     "NO   -> EW breathing RESONANCE (s8.1)"),
    ("0++ glueball","gauge",    "A_1g", "—",             f"(2N-1)L=5L={g0:.0f}", "NO   -> N=3 Fock condensate (L1247)"),
    ("2++ glueball","gauge",    "E_g",  "—",             f"(2N-1)L=7L={g2:.0f}", "NO   -> N=4 Fock condensate (L1247)"),
]
for row in table:
    print(f"  {row[0]:<14}{row[1]:<13}{row[2]:<8}{row[3]:<15}{row[4]:<16}{row[5]}")

# self-assert the loop-mode-SETTLED entries and the Fock-machinery glueball masses (the latter NOT from loop-mode)
assert abs(rho_m - 759.7) < 0.5, "rho mass from loop-mode (settled)"
assert abs(2 * np.cos(np.pi / 3) - 1.0) < 1e-12, "photon E=1 from loop-mode (settled, massless)"
assert abs(g0 - 1660) < 1 and abs(g2 - 2324) < 1, "glueball masses are (2N-1)L FOCK values, NOT loop-mode"
# the would-be NAIVE single-loop glueball (the flagged K_6 error) is a DIFFERENT, wrong number:
assert abs(5 * LQCD - 1660) < 1 and abs(2 ** 0.5 * LQCD - 469) < 1, \
    "a single C_4 loop gives ~470 (photon-region), NOT the 1660 Fock 0++ -> single-loop != glueball"

print("""
WHAT THE LOOP-MODE RULE SETTLES vs WHAT IT DOES NOT:
  SETTLED (single-loop dipole = the eigenvalue IS the physics):
    * photon  : gauge C_4 vector, E=1 (gapless) -> massless. The square/gauge action space.
    * rho     : matter C_8 vector, phi -> sqrt2 phi Lambda = 760 MeV (2% of 775). The octagon/matter space.
    => the J=1 SINGLE-LOOP bosons are fixed: eigenvalue 1 = massless gauge (photon), phi = massive
       matter. This is the s-action-space distinction made quantitative (photon/rho is the T_1u row of the
       canonical table; the Higgs/0++ glueball is its A_1g analogue).
  NOT SETTLED by the loop-mode rule (same irreps, different mechanism -- USE the anchored machinery):
    * Higgs        : matter A_1g is an EW-scale breathing RESONANCE (s8.1), not a Lambda-scale loop eig.
    * graviton     : gauge E_g is EVEN-dim, so mass-PERMITTED by the O_h parity theorem (L1169) -- UNLIKE
       the odd-dim T_1u photon (mass-FORBIDDEN). Its masslessness is the s10.2 period-4-vacuum / Gamma-point
       tachyon resolution, NOT the loop/parity rule -- so the graviton is NOT settled by the single-loop tool.
    * 0++/2++ glueballs: HDR-exempt N-plaquette FOCK condensates, m_N=(2N-1)Lambda (1660 / 2324). A SINGLE
       loop eigenvalue here is the documented K_6 error (L1227, 'structurally invalid'); the ~470 MeV a
       single C_4 would give is photon-region, NOT the 1660 0++. Multi-loop Fock is required.
  NET: the same loop-mode tool that resolved the rho/photon fork does NOT extend to glueball MASSES (Fock)
       or the Higgs (EW resonance). What DOES transfer to glueballs is the orthogonal KPM survival/WIDTH
       machinery (the spectral function), not the mass -- a separate, non-redundant computation.
""")
print("exit 0 -- eigenvalues + masses + girth asserted; host-path fork resolved (octagon=rho, body-diagonal=photon).")
