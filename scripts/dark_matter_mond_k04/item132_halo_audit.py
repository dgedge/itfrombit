#!/usr/bin/env python3
r"""ITEM 132 (dark halo EoS / cusp-core / BTFR) AUDIT -- self-asserting (exit 0 = every quoted number checked).
Companion to DRIFT M14. Canon item 132 (ANCHOR L3357-3373, anchored 2026-05-27):
  L3367: isotropic w_eff=-4/15 fluid gives a cored halo via "negative pressure -> outward support".
  L3369: BTFR v^4=a0 G M_b because "the halo topological tension is strictly proportional to M_b".
  L3370: a0 = c H0 / 2pi.
This shows: (1) L3367 fails the Jeans sign test; (2) L3369's LINEAR tension gives the WRONG BTFR power;
(3) the radial-string fix that works needs the "nonlinear flux law", which IS deep-MOND (so MOND is ASSUMED
unless the R4 action is derived; companion item132_p3_r4_action.py narrows this by forcing p=3 from
scale covariance once a local R4 action is granted, and item132_r4_local_action_lift.py reduces that
action premise to the R4 line-density law.  Companion item132_r4_line_density_dynamics.py shows
generic stable dynamics gives lambda_R4=(2/3)chi_R4|g|/a0; companion
item132_chi_unit_poisson.py gives chi_R4=1 under the Poisson QEC rate-matching lemma);
(4) a0=cH0/2pi is the a0<->horizon coincidence;
(5) the cored-profile MOND interpolation matches the RAR no better than a zero-parameter MOND function (§16.3).
"""
import numpy as np
c=2.998e8; G=6.674e-11; Mpc=3.086e22; Msun=1.989e30
ok=True
def check(cond,msg):
    global ok; ok = ok and bool(cond); print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")

print("[1] Isotropic w_eff<0 FAILS the Jeans sign test (refutes L3367 'negative pressure -> outward support'):")
# hydrostatic dp/dr=-rho GM/r^2 with p=w rho c^2 (w const) => drho/dr = -rho GM/(r^2 w c^2); c_s^2=dp/drho=w c^2
w=-4/15
check(-np.sign(w)>0, f"w_eff={w:.3f}<0 => drho/dr>0: density RISES outward (unphysical for a bound halo)")
check(np.sign(w)<0,  f"c_s^2 = w_eff c^2 < 0: imaginary sound speed -> dynamically unstable")

print("\n[2] Canon L3369 'tension proportional to M_b' gives the WRONG BTFR power:")
# flat rotation from rho=A/r^2 gives v^2=4 pi G A.  BTFR is v^4 ∝ M_b.
Mb=np.array([1e9,1e10,1e11])*Msun
p_lin =np.polyfit(np.log(Mb),np.log((4*np.pi*G*Mb)**2),1)[0]        # A∝M_b
p_sqrt=np.polyfit(np.log(Mb),np.log((4*np.pi*G*np.sqrt(Mb))**2),1)[0] # A∝sqrt(M_b)
check(abs(p_lin -2)<1e-9, f"LINEAR A∝M_b  => v^4 ∝ M_b^{p_lin:.2f}  (canon's mechanism -> WRONG, BTFR needs exponent 1)")
check(abs(p_sqrt-1)<1e-9, f"SQRT  A∝√M_b => v^4 ∝ M_b^{p_sqrt:.2f}  (BTFR) -- requires a NONLINEAR amplitude law")

print("\n[3] The 'nonlinear R4 Gauss law'  ∮(g^2/a0)dA = 4πG M_b  IS deep-MOND  (g^2 = a0 g_N):")
H0=67.4*1e3/Mpc; a0=c*H0/(2*np.pi)
r=np.logspace(18,22,50); Mb1=1e11*Msun
gN=G*Mb1/r**2; g_nl=np.sqrt(a0*G*Mb1)/r                  # spherical soln of the nonlinear flux law
check(np.allclose(g_nl**2, a0*gN), "g_4^2 = a0 g_N exactly == the deep-MOND field equation")
check(np.allclose((g_nl*r)**2, a0*G*Mb1), "=> v_inf^4 = a0 G M_b (BTFR) is the TEXTBOOK deep-MOND consequence, not a derivation")

print("\n[4] a0 = c H0 / 2pi  (the a0<->horizon coincidence), vs Milgrom 1.2e-10 m/s^2:")
check(abs(a0-1.04e-10)<0.03e-10, f"a0 = {a0:.3e} = {a0/1.2e-10*100:.0f}% of Milgrom at H0=67.4 (6-13% low over H0=67-73)")

print("\n[5] §16.3: the cored-profile nu(y) RAR match is NOT discriminating:")
y=np.logspace(-4,2,400)
nu_user =1+(1/np.sqrt(y))*(1-(np.sqrt(y)/3)*np.arctan(3/np.sqrt(y)))  # cored profile, r_c=r_M/3
nu_simple=0.5+np.sqrt(0.25+1/y)                                       # zero-param 'simple' MOND fn
nu_RAR  =1/(1-np.exp(-np.sqrt(y)))                                    # McGaugh 2016 empirical RAR
d_user  =np.max(np.abs(nu_user/nu_RAR -1))*100
d_simple=np.max(np.abs(nu_simple/nu_RAR-1))*100
check(1.5<d_user<2.5, f"nu_user vs RAR = {d_user:.1f}%  (inside the RAR's ~10% intrinsic scatter)")
check(d_simple<6,     f"zero-param 'simple' nu vs RAR = {d_simple:.1f}%  -> a ~2% match singles out nothing")
check(abs(nu_user[-1]-1)<0.01 and nu_user[0]>50, "nu_user limits: ->1 (Newtonian), ->1/√y (deep-MOND) == it IS a MOND fn")

print("\n"+"="*94)
print("VERDICT (item 132): the radial-string reframing FIXES canon's two errors (isotropic-Jeans L3367,")
print("linear-tension BTFR L3369) and supplies an explicit v_c(r) -- genuine improvements. BUT the BTFR/flat-")
print("curve result ASSUMES MOND: the load-bearing 'nonlinear flux law' = deep-MOND (g^2=a0 g_N), and BTFR is")
print("its textbook consequence. Companions item132_p3_r4_action.py and item132_r4_local_action_lift.py")
print("narrow the substrate target: p=3 is forced; item132_r4_line_density_dynamics.py shows generic")
print("stable dynamics gives lambda_R4=(2/3) chi_R4 |g|/a0. item132_chi_unit_poisson.py gives")
print("chi_R4=1 if the boundary-QEC line ledger has matched creation/erasure rates and Poisson")
print("count support. Full closure still requires deriving those lemmas intrinsically, plus")
print("a core-boundary axiom (the 1/3) and the a0<->horizon coincidence. Tier:")
print("cusp-core/flat-curve mechanism CONSISTENT WITH MOND; MOND-from-substrate is the open problem. NOT")
print("a parameter-free derivation.")
print("="*94)
assert ok, "an item-132 audit check FAILED"
print("exit 0 -- two canon errors confirmed; BTFR/MOND assumed-not-derived; reframing is a real mechanism fix.")
