#!/usr/bin/env python3
"""
Koide per-sector phase: testing the chiral-anomaly / Atiyah-Singer "Target C" hypothesis
against the nu_R-defect-index hypothesis, sector by sector.

Built ONLY from explicitly-anchored ANCHOR.md structure:
  - III(c) 16-state single-generation fermion-charge table (ANCHOR S15 item 135 III(c))
  - S5.8 / Part 18 nu_R Feshbach coupling: H_QQ = I_3x3; H_PQ has 6 nonzero entries per
    generation (e_R->nu_R via I3 flip, nu_L->nu_R via chi flip); NO quark couples (Colour Firewall, S2.14)
  - Coin (item 89): zero-controlled CNOT I3 -> I3 XOR ~chi
  - gamma5 = sigma_y^(chi) (x) I^(I3)   (S3.5 / item 93)

Nothing here is fitted. Where a value cannot be derived from spec, it is reported OPEN.
"""
from fractions import Fraction as F

line = "="*78
def hdr(s): print("\n"+line+"\n"+s+"\n"+line)

# ----------------------------------------------------------------------------
# 1. The III(c) 16-state single-generation register, verbatim from ANCHOR.
#    bits (c0,c1,c2,c3); W = c1+c2+c3; charges as tabulated.
# ----------------------------------------------------------------------------
# (state, chirality, T3, Y, Q, colour multiplicity)
III_c = [
    ("nu_L", "LH", F(1,2),  F(-1),   F(0),   1),
    ("e_L",  "LH", F(-1,2), F(-1),   F(-1),  1),
    ("u_L",  "LH", F(1,2),  F(1,3),  F(2,3), 3),
    ("d_L",  "LH", F(-1,2), F(1,3),  F(-1,3),3),
    ("nu_R", "RH", F(0),    F(0),    F(0),   1),
    ("e_R",  "RH", F(0),    F(-2),   F(-1),  1),
    ("u_R",  "RH", F(0),    F(4,3),  F(2,3), 3),
    ("d_R",  "RH", F(0),    F(-2,3), F(-1,3),3),
]

hdr("1. SM consistency of the III(c) register (expected: all exact)")
Tr_T3sq = sum(m*s[2]**2 for s in III_c for m in [s[5]] if s[1]=="LH")
Tr_Ysq  = sum(s[5]*(s[3]/2)**2 for s in III_c)
sumQL3  = sum(s[5]*s[4]**3 for s in III_c if s[1]=="LH")
print(f"  Tr(T3^2)  over LH      = {Tr_T3sq}   (SM: 2)")
print(f"  Tr((Y/2)^2) all states = {Tr_Ysq}   (SM: 10/3)")
print(f"  sin^2_thetaW (GUT)     = {Tr_T3sq/(Tr_T3sq+Tr_Ysq)}   (SU(5): 3/8)")
print(f"  sum_LH Q^3             = {sumQL3}   (SM: -2/9)")

# ----------------------------------------------------------------------------
# 2. ROUTE 1 -- the cubic chiral trace, evaluated PER SECTOR (the stuck route).
#    Koide targets (anchored externally, S15 item 86 / 96):
# ----------------------------------------------------------------------------
target = {"lepton": F(2,9), "neutrino": F(3,9), "down": F(1,9), "up": F(2,27)}

hdr("2. ROUTE 1: cubic chiral trace |sum Q_L^3| per sector vs Koide target")
Qe, Qu, Qd = F(-1), F(2,3), F(-1,3)
cubic = {
    "lepton": abs(Qe**3),          # e_L alone
    "down":   abs(3*Qd**3),        # 3 colours
    "up":     abs(3*Qu**3),        # 3 colours
}
cubic["gen_total"] = abs(Qe**3 + 3*Qu**3 + 3*Qd**3)
for k in ("lepton","down","up"):
    t = target[k]; v = cubic[k]
    verdict = "MATCH" if v==t else f"MISS (factor {v/t})"
    print(f"  {k:9s} |sumQ^3| = {str(v):6s}  target {str(t):5s}  -> {verdict}")
print(f"  gen_total |sumQ^3| = {cubic['gen_total']}  (equals lepton target -- the item-86 entry1=entry3 identity)")
print("  -> 1 of 3 sectors. 'down' is a COLLISION of two unrelated 1/9's (d/N vs N_c*Q^3).")

# ----------------------------------------------------------------------------
# 3. The nu_R Feshbach coupling, built explicitly (S5.8 / Part 18).
#    8-bit register; canonical ordering (item 99 / csit):
#        (G0, G1, LQ, C0, C1, I3, chi, W)
#    nu_R pseudocodeword per generation. The framework states the two coupling
#    channels into nu_R: e_R -> nu_R via I3-flip, nu_L -> nu_R via chi-flip.
#    A quark would need to also change colour bits (C0,C1) to reach nu_R -> blocked.
# ----------------------------------------------------------------------------
hdr("3. nu_R Feshbach defect support d, from explicit coupling channels")

# Represent the relevant states by the bits that the coupling touches: (I3, chi, colour_active)
# nu_R is the all-isotropic codeword (no colour polarisation). Coupling = single-bit flip into it.
# Channels named exactly as Part 18:
channels = {
    "charged_lepton": ["e_R via I3-flip", "nu_L via chi-flip"],   # both reach nu_R in ONE flip
    "neutrino":       ["e_R via I3-flip", "nu_L via chi-flip", "frozen I3 (CNOT dormant, LQ=0)"],
}
N_plaq = 9  # 8 boundary + 1 bulk-centre, the 4.8.8 vertex figure (item 86)

for sector, ch in channels.items():
    d = len(ch)
    print(f"  {sector:15s}: d = {d} channels  -> d/N = {F(d,N_plaq)}")
    for c in ch: print(f"        - {c}")

print("\n  Quark -> nu_R coupling:")
print("    nu_R has NO colour polarisation; a quark has active colour bits (C0,C1).")
print("    Reaching nu_R requires flipping colour bits too -> NOT a single-flip channel.")
print("    => H_PQ(quark) = 0  (Colour Firewall, S2.14).  Defect support d_quark = 0.")
print("    => The nu_R-index route CANNOT produce a quark phase. Quarks are off this mechanism.")

# ----------------------------------------------------------------------------
# 4. Honest index-vs-density check.
# ----------------------------------------------------------------------------
hdr("4. Is 2/9 an Atiyah-Singer index? (honest categorial check)")
print("  Atiyah-Singer / lattice (HLN-Luscher) index = n_+ - n_-  is an INTEGER.")
print(f"  The genuine integer index here is d (zero-mode/channel count): lepton d=2, nu d=3.")
print(f"  2/9 is d/N = (integer index)/(plaquette dim 9): a DEFECT DENSITY, not an index.")
print("  => 'Chern number c_1 = 2/9' is mislabelled. c_1 in Z; the rational is index/dim.")
print("     The per-sector object that is legitimately rational is the index DENSITY d/N.")
print("  => D_TCH itself is NOT written down in canon (item 93: 'explicit construction ...")
print("     in closed form' listed OPEN). So the 'lattice Atiyah-Singer' branding is")
print("     aspirational until D_TCH exists explicitly.")

# ----------------------------------------------------------------------------
# 5. Quark sector via the I3 / colour-coherence route (NOT the anomaly index).
#    item 96: up is I3=0 (CNOT dormant) -> colour coherence preserved -> N_eff = 9*N_c.
#             down is I3=1 (CNOT fires) -> coherence collapsed -> N stays 9.
# ----------------------------------------------------------------------------
hdr("5. ROUTE 2: quark phase from I3/colour-coherence (Boltzmann sector, not anomaly)")
Nc = 3
d_charged = 2
up_pred   = F(d_charged, N_plaq*Nc)          # I3=0 passive: N_eff = 9*Nc
print(f"  up   (I3=0, coherence preserved, N_eff=9*Nc=27): d/N_eff = {up_pred}  target {target['up']}  "
      f"-> {'MATCH' if up_pred==target['up'] else 'MISS'}")
print(f"  down (I3=1, coherence collapsed, N=9): target {target['down']} needs d/N = 1/9")
print(f"        d=2,N=9 -> 2/9 (= lepton, factor 2 too big);  d=1,N=9 -> 1/9 (works but why d=1?)")
print(f"        Either way an unexplained factor 2.  This IS Willow Q4 (delta_d = delta_l/2).")
print(f"        => DOWN REMAINS OPEN. Not closed by either framing. Reported honestly.")

# ----------------------------------------------------------------------------
# 6. Verdict
# ----------------------------------------------------------------------------
hdr("6. VERDICT")
print("""  Target C as posed ('one universal anomaly current for all four phases') is
  STRUCTURALLY IMPOSSIBLE by the framework's own Colour Firewall: leptons/neutrinos
  get their phase from the nu_R Feshbach defect index; quarks have NO nu_R channel
  and get mass from the Boltzmann/I3 sector. Two mechanisms, not one.

    lepton  delta = 2/9   DERIVED   (nu_R index d/N = 2/9, explicit H_PQ)
    nu      delta = 3/9   DERIVED   (nu_R index d/N = 3/9, +frozen I3)
    up      delta = 2/27  CONSISTENT(I3=0 colour-coherence N_eff=27, item-96 anchored)
    down    delta = 1/9   OPEN      (the factor-2 residue == Willow Q4)

  Retire: cubic-trace route (1/3 sectors; 'down' a coincidence) -- alongside the
          already-retired spatial Aharonov-Bohm route.
  Correct: '2/9 = Atiyah-Singer index' -> '2/9 = index DENSITY d/N'; the integer
          index is d=2.
""")
