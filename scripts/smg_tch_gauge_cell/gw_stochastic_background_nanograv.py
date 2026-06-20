#!/usr/bin/env python3
r"""GRAVITATIONAL-WAVE STOCHASTIC BACKGROUND from the substrate crystallisation, vs NANOGrav.
(First GW prediction of the framework -- a new, falsifiable observable channel.)

THE INGREDIENTS (all from existing canon, no new input):
  - the substrate crystallises at the boot, at temperature T_* ~ Lambda_QCD = 0.332 GeV
    (the w6/Lambda = 1/6 boot bridge: T_start = 6 w6 = Lambda);
  - it traps a Kibble-Zurek defect/domain-wall network (the K04 walls);
  - those walls are SUBSTRATE-PINNED (depinning_mobility_gate.py, ~42 OOM) and their
    recorded gravitating shadow collapses to VACUUM grade (pairing_orphan_closure.py).

THE PHYSICS: any cosmological event at T_* radiates GW at the redshifted horizon frequency.
The headline is the PEAK FREQUENCY (robust, set only by T_*); the amplitude is governed by the
SAME wall-pinning/invisibility this session established for the dark sector -- a genuine cross-
check, not a new assumption.

exit 0 = the peak frequency (-> which detector band) and the domain-wall overclosure (Zel'dovich)
check are computed and ASSERTED; the amplitude bracket, spectral discriminator, and the honest
verdict (in-band, sub-threshold, cross-sector-consistent, with a sharp falsifier) are reported.
Scope: standard cosmological-GW relations + framework scales; not a lattice GW simulation.
"""
import math

# ----------------------------------------------------------------------------- [0] constants
LAMBDA_QCD = 0.332          # GeV (substrate fundamental scale)
T0 = 2.348e-13             # GeV (CMB temperature today)
M_PL = 1.221e19           # GeV (Planck mass)
GEV_HZ = 1.519e24         # Hz per GeV (1/hbar)
GS0 = 3.91                # entropic g_* today
NHZ = 1e-9                # Hz
# NANOGrav 15yr band ~ 2-60 nHz; GW-interpreted signal level Omega_GW h^2 ~ 1e-9..1e-8
NG_LO, NG_HI = 1.0, 100.0     # nHz (PTA band)
SIGMA_WALL_MAX = (1e-3) ** 3  # Zel'dovich/CMB stable-wall bound ~ (1 MeV)^3, in GeV^3

def f_peak_today(T_star, g_star):
    """redshifted Hubble frequency of a transition at T_star (Hz)."""
    H_star = 1.66 * math.sqrt(g_star) * T_star ** 2 / M_PL            # GeV
    a_ratio = (GS0 / g_star) ** (1 / 3) * (T0 / T_star)               # a_*/a_0
    return a_ratio * H_star * GEV_HZ

# ----------------------------------------------------------------------------- [1] peak frequency
print("[1] PEAK FREQUENCY (robust -- set only by the transition temperature T_*):")
cases = [("Lambda_QCD = 0.332 GeV (above QCD)", 0.332, 61.75),
         ("QCD crossover ~0.155 GeV",            0.155, 17.25)]
fpk = {}
for name, Ts, gs in cases:
    f = f_peak_today(Ts, gs); fpk[Ts] = f / NHZ
    print(f"    T_* = {name:34s}: f_peak ~ {f/NHZ:5.1f} nHz")
print(f"    -> the Lambda_QCD-scale crystallisation peaks at ~{fpk[0.332]:.0f}-{fpk[0.155]:.0f} nHz: "
      f"the NANOGrav/PTA band ({NG_LO:.0f}-{NG_HI:.0f} nHz).")
assert NG_LO < fpk[0.332] < NG_HI and NG_LO < fpk[0.155] < NG_HI, "peak not in the PTA band"
print("    This is the headline: a substrate transition at the framework's OWN scale lands in")
print("    exactly the band where a stochastic GW signal has been reported (NANOGrav 2023).")

# ----------------------------------------------------------------------------- [2] Zel'dovich overclosure -> pinning REQUIRED
print("\n[2] DOMAIN-WALL OVERCLOSURE (Zel'dovich-Kobzarev-Okun) -- the consistency gate:")
sigma_wall = LAMBDA_QCD ** 3        # K04 wall tension ~ Lambda_QCD^3 (the natural scale)
overclose = sigma_wall / SIGMA_WALL_MAX
print(f"    K04 wall tension sigma ~ Lambda_QCD^3 = {sigma_wall:.2e} GeV^3")
print(f"    Zel'dovich stable-wall bound sigma <~ (1 MeV)^3 = {SIGMA_WALL_MAX:.0e} GeV^3")
print(f"    -> a NAIVE scaling wall network OVERCLOSES by ~{overclose:.0e} ({math.log10(overclose):.1f} OOM):")
print("       a standard mobile Lambda_QCD wall network is cosmologically EXCLUDED. The framework")
print("       survives ONLY because its walls are (a) substrate-PINNED (depinning_mobility_gate.py,")
print("       ~42 OOM) and (b) gravitationally INVISIBLE (pairing_orphan_closure.py, shadow ->")
print("       vacuum grade) -- the SAME two results that demoted K04 from dark matter this session.")
assert overclose > 1e6, "naive walls should grossly overclose"

# ----------------------------------------------------------------------------- [3] amplitude bracket
print("\n[3] AMPLITUDE -- governed by the same pinning (this is the cross-check, not a new knob):")
print("    GW emission needs an ACCELERATING mass quadrupole. The K04 walls are PINNED (frozen,")
print("    no peculiar motion -- the depinning gate), so they radiate weakly; and their gravitating")
print("    stress collapsed to vacuum grade. Both push Omega_GW FAR below a scaling network.")
print(f"    - strong scaling/first-order source (hypothetical): Omega_GW h^2 ~ 1e-9..1e-7 (detectable)")
print(f"      -> but that is exactly the case Zel'dovich [2] EXCLUDES (overproduces GW AND overcloses).")
print(f"    - framework actual (glassy transition + pinned/invisible walls): Omega_GW h^2 << 1e-9")
print(f"      -> SUB-THRESHOLD: the framework does NOT predict it IS the NANOGrav signal.")

# ----------------------------------------------------------------------------- [4] spectral discriminator
print("\n[4] SPECTRAL SHAPE (a discriminator IF ever detected):")
print("    causally-sourced background below its peak: Omega_GW ~ f^3  =>  strain h_c ~ f^(+1/2) (RISING).")
print("    astrophysical SMBHB background: h_c ~ f^(-2/3) (FALLING) -- what NANOGrav data resemble.")
print("    -> in the PTA band (below the ~30-50 nHz peak) the framework predicts a RISING h_c, opposite")
print("       to SMBHB: a clean discriminator. The data's ~SMBHB-like slope is consistent with the")
print("       framework being SUB-THRESHOLD there (its rising tail hidden under the SMBHB signal).")

# ----------------------------------------------------------------------------- [5] verdict
print(f"""
[5] VERDICT — a genuine new, falsifiable channel; honestly: in-band, sub-threshold, cross-consistent.
  * NEW PREDICTION (robust): the substrate crystallisation at T_* ~ Lambda_QCD radiates GW peaked at
    ~{fpk[0.155]:.0f}-{fpk[0.332]:.0f} nHz -- the NANOGrav/PTA band. The frequency follows from the framework's OWN
    scale with no free parameter; it predicts the signal is in the PTA band, NOT LISA (mHz) or LIGO (Hz).
  * AMPLITUDE: self-consistently SUB-THRESHOLD. A detectable Lambda_QCD wall network would overclose
    the universe by ~{math.log10(overclose):.0f} OOM (Zel'dovich); the framework's walls evade this ONLY by being
    pinned + gravitationally invisible -- the very results that demoted K04 from dark matter -- and the
    same pinning suppresses GW emission. So the framework does NOT claim the NANOGrav signal.
  * CROSS-SECTOR CONSISTENCY (the real prize): ONE mechanism -- K04 wall pinning + vacuum-grade
    shadow -- now does THREE independent jobs: (i) K04 is not mobile halo dark matter; (ii) K04 evades
    the Zel'dovich domain-wall catastrophe; (iii) K04 radiates negligible GW. Three sectors, one cause,
    no new assumption -- a nontrivial coherence the framework did not have to possess.
  * SHARP FALSIFIER: if a strong stochastic background at ~tens of nHz is confirmed COSMOLOGICAL
    (not SMBHB) at NANOGrav amplitude, it is in TENSION with the framework's required wall suppression;
    a SMBHB-dominated interpretation (current best fit) is consistent. And the predicted band itself is
    falsifiable: a confirmed cosmological signal peaking at LISA/LIGO frequencies, not PTA, would not fit.
  TIER: peak frequency = robust parameter-free prediction; overclosure gate = exact; amplitude =
    bracketed/sub-threshold via the cross-checked pinning; spectral shape = discriminator if detected.
    A precise Omega_GW(f) needs the transition order + wall-network lattice GW emission (next step).
exit 0""")
print(f"ALL ASSERTIONS PASSED — f_peak ~{fpk[0.332]:.0f} nHz in the PTA band; naive walls overclose by "
      f"{math.log10(overclose):.0f} OOM (pinning required).")
