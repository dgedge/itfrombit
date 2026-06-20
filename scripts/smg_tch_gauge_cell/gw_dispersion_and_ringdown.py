#!/usr/bin/env python3
r"""GW DISPERSION + DISCRETE-AREA RINGDOWN — the two sibling GW predictions (item 152 covered
only the stochastic background). Honest question: are they observable, or null?

[A] GW DISPERSION (graviton = shear phonon of the Z3 lattice): same lattice as the photon
    (item 150), so inversion forbids linear LV and the first anisotropy is O(k^4),
    Delta v/v ~ (a0 k)^2/12, with the same Brillouin cutoff ~Lambda_QCD. The photon's
    3e-18 was for OPTICAL wavelengths; GW wavelengths are ~10^12 longer, so the SAME law
    gives a far smaller number. Compute it and compare to the GW170817 speed-of-gravity bound.
[B] DISCRETE-AREA RINGDOWN (S=A/4 + the 55/8 ledger): the horizon area is quantized
    (Bekenstein-Mukhanov, dA = 4 ln2 ell_P^2 per minimal record), giving quantized QNMs /
    possible echoes. Compute the area-quantum number for a stellar BH, the QNM shift, and the
    (conditional) echo delay.

exit 0 = both effects computed and compared to bounds; the verdict (real but observationally
null / consistency-grade, with one conditional echo caveat) is reported honestly. Scope:
order-of-magnitude GW phenomenology + the framework's lattice/area scales.
"""
import math

# ----------------------------------------------------------------------------- [0] constants
C = 2.99792458e8           # m/s
HBAR = 1.054571817e-34     # J s
G = 6.674e-11              # SI
M_SUN = 1.989e30           # kg
A0 = 0.594e-15             # lattice spacing (m)
ELL_P = 1.616e-35          # Planck length (m)
LAMBDA_QCD_HZ = 0.332 * 1.602e-10 / (2 * math.pi * HBAR)   # 0.332 GeV in J / (2 pi hbar) -> Hz
GW170817_BOUND = 1e-15     # |c_g/c - 1| bound from the GW+GRB time delay (representative)
print(f"[0] a0={A0:.2e} m ; ell_P={ELL_P:.2e} m ; graviton BZ cutoff ~ {LAMBDA_QCD_HZ:.1e} Hz")

# ----------------------------------------------------------------------------- [A] GW dispersion
def dv_over_v(f_hz):
    k = 2 * math.pi * f_hz / C       # GW wavenumber
    return (A0 * k) ** 2 / 12.0      # O(k^4) anisotropy, same law as the photon
print("\n[A] GW DISPERSION (graviton = shear phonon; same O(k^4) law as the photon):")
for label, f in (("LIGO 100 Hz", 1e2), ("LISA 1 mHz", 1e-3), ("PTA 30 nHz", 3e-8)):
    print(f"    {label:14s}: Delta v/v = {dv_over_v(f):.2e}")
dv_ligo = dv_over_v(1e2)
print(f"    GW170817 speed-of-gravity bound |c_g/c-1| < ~{GW170817_BOUND:.0e}")
print(f"    -> framework at LIGO band is {GW170817_BOUND/dv_ligo:.0e}x BELOW the bound: passes by "
      f"~{math.log10(GW170817_BOUND/dv_ligo):.0f} OOM.")
assert dv_ligo < GW170817_BOUND
# frequency at which the effect would reach the bound:
f_reach = math.sqrt(12 * GW170817_BOUND) / (2 * math.pi * A0 / C)
print(f"    Delta v/v reaches {GW170817_BOUND:.0e} only at f ~ {f_reach:.1e} Hz (UV-wavelength gravitons,")
print(f"    ~{f_reach/1e9:.0e} GHz) -- no such GW source exists. Cutoff at ~{LAMBDA_QCD_HZ:.0e} Hz is unreachable.")
print("    VERDICT: real (tensor analogue of the photon dispersion) but UNOBSERVABLE -- GW wavelengths")
print("    (~km) are ~10^12x longer than optical, so the same law gives ~10^-43 not ~3e-18. Null/consistent.")

# ----------------------------------------------------------------------------- [B] discrete-area ringdown
print("\n[B] DISCRETE-AREA RINGDOWN (S=A/4 + 55/8 ledger -> Bekenstein-Mukhanov area quantum):")
M = 30 * M_SUN                                   # typical LIGO remnant
R_bh = 2 * G * M / C ** 2                        # Schwarzschild radius
A_bh = 4 * math.pi * R_bh ** 2
dA = 4 * math.log(2) * ELL_P ** 2                # area quantum (minimal record dS=ln2)
N_quanta = A_bh / dA
qnm_shift = (ELL_P / R_bh) ** 2                  # fractional QNM deviation ~ (ell_P/R)^2
echo_delay = (R_bh / C) * math.log(R_bh / ELL_P) # IF the horizon reflects (Planck-membrane echo)
print(f"    stellar BH (30 Msun): R={R_bh/1e3:.0f} km, area A={A_bh:.1e} m^2")
print(f"    area quantum dA = 4 ln2 ell_P^2 = {dA:.1e} m^2  ->  N_quanta = A/dA = {N_quanta:.1e}")
print(f"    fractional QNM shift from quantization ~ (ell_P/R)^2 = {qnm_shift:.1e}  -> NEGLIGIBLE")
assert qnm_shift < 1e-50
print(f"    -> the horizon is ~10^{math.log10(N_quanta):.0f} quanta: discreteness is invisible in the")
print("       continuous QNM spectrum. Default prediction = STANDARD GR ringdown (consistent with LIGO).")
print(f"    CONDITIONAL echo: IF the QEC horizon reflects coherently (a Planck-membrane, NOT established")
print(f"    -- the discreteness is sub-horizon), echoes would appear at dt ~ {echo_delay*1e3:.0f} ms after")
print("       ringdown (LIGO band). This is the only potentially-observable handle, and it is conditional.")

# ----------------------------------------------------------------------------- [C] verdict
print(f"""
[C] VERDICT — both sibling predictions are REAL but observationally NULL (the answer to 'did we
    address them': now yes, and they do not yield discriminating tests):
  * GW DISPERSION: the graviton inherits the photon's lattice dispersion (no linear LV, O(k^4)
    anisotropy, ~{LAMBDA_QCD_HZ:.0e} Hz cutoff), but at GW wavelengths Delta v/v ~ {dv_ligo:.0e} (LIGO band) --
    ~{math.log10(GW170817_BOUND/dv_ligo):.0f} OOM below the GW170817 bound. The framework passes the speed-of-gravity
    test trivially; there is no observable GW dispersion (the 3e-18 was optical-wavelength-specific).
  * DISCRETE-AREA RINGDOWN: the area IS quantized (Bekenstein-Mukhanov, dA=4 ln2 ell_P^2), but a
    stellar horizon holds ~10^{math.log10(N_quanta):.0f} quanta, so the QNM shift ~(ell_P/R)^2 ~ {qnm_shift:.0e} is
    unmeasurable. Default = standard GR ringdown. A ~{echo_delay*1e3:.0f} ms echo is the ONLY possible signature,
    and only IF the QEC horizon reflects coherently -- not established (the discreteness is sub-horizon).
  NET: the framework's DISCRETENESS is invisible in the GW propagation/ringdown sector -- it passes
    GW170817 and LIGO ringdown by 28-76 OOM (a strong consistency, not an overclaim). So the GW sector
    has exactly ONE useful channel: the STOCHASTIC BACKGROUND (item 152, in the NANOGrav/PTA band).
    The honest correction to the earlier 'two sibling GW predictions': real, but not tests. The one
    live handle is the conditional ringdown echo, gated on a reflecting-QEC-horizon theorem.
exit 0""")
print(f"ALL ASSERTIONS PASSED — GW dispersion {dv_ligo:.0e} (28 OOM below GW170817); ringdown shift "
      f"{qnm_shift:.0e} (negligible); one channel (stochastic background) survives.")
