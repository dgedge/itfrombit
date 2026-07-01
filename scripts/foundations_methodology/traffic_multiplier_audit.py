#!/usr/bin/env python3
r"""TRAFFIC-MULTIPLIER AUDIT (R1 follow-up) — the parameter ledger for the record-action measure.

Full record-action form (the +Sum_a lambda_a Q_a constraint term, restored):

    A_record = A_traffic - (1/2) dS_record + i Phi_rec + Sum_a lambda_a Q_a .

R1 factors the traffic sector as  A_traffic = Omega_clock * T_rel[gamma].  The predictive criterion
is whether  A_traffic = Omega_sat * T_code[gamma]  with T_code built from finite code counts / orbit
incidences / stabiliser constraints / graph degrees / scheduler capacities -- NOT fitted reals.
This audit does NOT attempt any new derivation; it LISTS every open coefficient, CLASSIFIES each
from the *current* canon status (verified, cited), and COUNTS survivors so the ledger inequality
  K_traffic + K_phase << N_constants  can be evaluated, with the constraint/source columns kept
SEPARATE from the continuous-multiplier count.

STATUS classes:
  forced       -- reduced to code counts / group invariants / a closed theorem. NOT a multiplier.
  conditional  -- code-structured but forced ONLY under a named, not-yet-discharged premise.
  free         -- a fitted/continuous real coefficient or class-2 candidate (the danger).
  premise      -- a BINARY structural question (does a charge/zero-mode exist?), not continuous.

SECTORS / how each counts (per the corrected taxonomy):
  clock        -- Omega, the one absolute activity scale (saturation target).         K_clock
  traffic-core -- cosmological/QEC service relative activities.                        K_rel
  traffic-EW   -- electroweak stiffness ratios (item 55) -- a SECOND absolute scale.   K_rel
  traffic-nuc  -- nuclear binding stiffnesses (item 113) -- the nuclear scale.         K_rel
  holonomy     -- Im A_record: CP phases / signs.                                      K_phase
  constraint   -- Q_a existence (binary; M15 dust). NOT a continuous multiplier.       K_Q
  source       -- if a charge exists, what populates it? multiplier ONLY if fitted.    K_source(free)

Every 'where' is a canon location verified 2026-06-19 (DRIFT line / script).
"""

# coef, sector, status, where (verified citation), note
LEDGER = [
    # ---- the absolute clock (saturation leg) -------------------------------------------------
    ("Omega: H0 = Lambda/N_lock (cosmological endpoint)", "clock", "conditional",
     "DRIFT 'frontier snapshot' D + cosmological_selector_lock_theorem.py",
     "selector/service-span lock at conditional-theorem grade, -0.06%; this IS Omega_sat for 1 member"),

    # ---- traffic-core: cosmological / QEC service relative activities -------------------------
    ("K04 w4/w6 = 2 (orbit normal-ordering)", "traffic-core", "conditional",
     "DRIFT K04 orbit-deriv (l.99) k04_w4_w6_orbit_derivation.py",
     "w4/w6=2 RECOVERED from Q3 loop-orbit edge-incidence (w4=1/2, w6=1/4); conditional on the "
     "edge-ledger normal-ordering premise. (w6/Lambda + cooling law still open; SIGNS forced separately, l.95)"),
    ("CC rho_Lambda coefficient (generation-vertex loop)", "traffic-core", "conditional",
     "DRIFT l.3910-3922 cc_generation_vertex_item115_loop.py + cc_clause24_proof.py",
     "CONVERTED model-grade: C_loop=(1/3)<1/|sin k|>_BZ=0.30356 (one-insertion RS shift; alpha_0^1 billing "
     "FORCED); rho_Lambda=rho_obs at 0.10sigma; -0.39% residual ill-posed. Balance-form A_req=1.219=sqrt(3/2) "
     "REJECTED (cross-sector gravity C_loop; rate-multiplicity, not a determinant) -- handoff_action_derivation.py"),
    ("gravity coefficient (alpha^1 / C_loop=3/2)", "traffic-core", "conditional",
     "DRIFT G7 l.224 item119_jump_channel.py + g7_constructibility_audit.py",
     "alpha^2/K_eff Feshbach trace NOT constructible (4 invented ops + E-tuning, G7) AND structurally "
     "impossible (item119); SUPERSEDED by the alpha^1 erasure route M_P^2=(2/3)alpha Lambda^3 R_dS, "
     "C_loop=3/2 from canonical jump operators, conditional on item 79. Absolute M_P scale = horizon-input (selector lock)"),
    ("HBC scalar amplitude A_nu = (3/4) alpha0^4", "traffic-core", "conditional",
     "DRIFT l.686 item131_hbc_channel_whitening_closure.py + hbc_amplitude_status.py",
     "CONDITIONALLY closed (2026-06-15): F_eff=1, Pi_k closed, S_j(k=aH)=1 (local ledgers), latch "
     "lambda_shell=C_F via queue balance -> A_nu=(3/4)alpha0^4=2.13e-9. Residual: 3/4=1/C_F selection "
     "(T8), early H_* absolute scale (item 42), locality premise. R_HBC=psi-nu is the derived gauge-inv field"),
    ("BH flux coefficient 10/27", "traffic-core", "conditional",
     "DRIFT B5 (l.1989) bh_flux_scheduler_stencil_theorem.py",
     "10/27 = (9 Moore face + 1 latch)/27 IFF the item-120 Landauer-Moore alphabet is universal for 3D erasure"),
    ("HBC saturated-printer latch lambda_shell = C_F", "traffic-core", "conditional",
     "DRIFT snapshot C 2016",
     "VALUE is the Casimir C_F (group-forced); the marginality/stop mechanism is the open dynamical axiom"),
    ("item-126 T3 absolute-eta normalisation", "traffic-core", "conditional",
     "DRIFT item126 T3 + snapshot 2021, scheduler_alpha_composition.py",
     "c_4=alpha_0 DERIVED, s_1=ln(8x137) DERIVED; residual = the O(10) photon-bookkeeping convention only"),

    # ---- traffic-EW: electroweak stiffness ratios (a SECOND absolute scale, §16.2) ------------
    ("Higgs M_H/M_Z: matter-cell stiffnesses (55d, <=4)", "traffic-EW", "free",
     "DRIFT item55 (l.2059) item55_matter_cell_dynamical_matrix.py",
     "12->=4 central-force constants from Lambda_QCD UNPINNED; breathing mode 266x below 125 GeV"),
    ("Higgs M_H/M_Z: bridge k_shear, k_mix", "traffic-EW", "free",
     "DRIFT item55 (l.2061) item55_gauge_bridge_mz_ratio.py",
     "M_H, M_Z proven to share NO parameter -> ratio irreducibly cross-space, not parameter-free"),
    ("Higgs M_H/M_Z: Cauchy-Born k_t/k_l + coupling g", "traffic-EW", "free",
     "DRIFT item55f-ii (l.2065) item55f_cauchy_born_homogenization.py",
     "isotropy is automatic/geometric so it does NOT fix k_t/k_l (Poisson 0.83->0); stays free"),

    # ---- traffic-nuc: nuclear SEMF binding stiffnesses (the nuclear scale, §16.2) -------------
    ("nuclear a_V (bulk binding eps~5.25 MeV/bond)", "traffic-nuc", "free",
     "DRIFT item113 a_V (l.2073) item113_aV_tch_coordination.py",
     "z=6 coordination geometric; eps unpinned (Lambda/64 is a forking-path coincidence, not adopted)"),
    ("nuclear a_S (surface)", "traffic-nuc", "free",
     "DRIFT item113 (l.2069/2077) item113_discrete_liquid_drop.py", "many-body, not pairwise; form only"),
    ("nuclear a_A (asymmetry)", "traffic-nuc", "free",
     "DRIFT item113 (l.2077)", "I3/chi bipartite-matching combinatorics; half-closed"),
    ("nuclear a_P (pairing)", "traffic-nuc", "free",
     "DRIFT item113 (l.2071)", "no value computed"),

    # ---- holonomy: Im A_record (CP phases / signs) -------------------------------------------
    ("quark CKM CP phase + Jarlskog sign", "holonomy", "conditional",
     "DRIFT baryon-sign (l.2037) ckm_walk_signed_template.py",
     "J>0 (correct) from the e^{ik pi/4} walk orientation; needs the physical basis (RG/item 88) selected"),
    ("delta_nu (PMNS CP / leptogenesis / baryon SIGN)", "holonomy", "free",
     "DRIFT item87 (l.2055) item87_MR_derivation_attempt.py",
     "EXACT no-go: nu_e=0x00 is an inert walk eigenstate -> delta_nu is a GEOMETRIC PRIMITIVE, not walk-derived"),

    # ---- constraint: Q_a existence (BINARY premise -- NOT a continuous multiplier) ------------
    ("CMB pressureless dust reservoir omega~0.096 (does Q_a exist?)", "constraint", "premise",
     "DRIFT snapshot E 2025 (M15)",
     "BINARY: does the measure admit a conserved pressureless zero-mode? all native escapes closed; "
     "live route = R4 zero-mode or external AeST import; 3rd-peak falsification pressure. K_Q, not K_rel"),

    # ---- source: if the charge exists, what populates it? (multiplier ONLY if fitted) ---------
    ("nu_R + dust ABSOLUTE abundance (alpha_0/208 boot source)", "source", "conditional",
     "DRIFT snapshot 2020, item123_nuR_absolute_density_boot_qec.py",
     "code-STRUCTURED form (alpha_0 derived; 208 = invalid-horizon subspace count, cf B1) but the source "
     "LAW is the boot-QEC residual; RELATIVE 4:1 split is forced (below). NOT fitted -> not a free multiplier"),

    # ---- a sample of what IS forced (resolved -> context; NOT counted) ------------------------
    ("3/14 baryon branching", "traffic-core", "forced",
     "DRIFT item126 channel ledger, item126_channel_ledger.py", "AG(3,2) hyperplane count; 1/14 measure derived"),
    ("c_4 = alpha_0 (item-126 composition)", "traffic-core", "forced",
     "DRIFT item126 T3, scheduler_alpha_composition.py", "Gamma_4=alpha_0^5 Lambda, derived not free"),
    ("s_1 = ln(8x137) record unit", "traffic-core", "forced",
     "ANCHOR §16.4, record_alphabet_derivation.py", "address x channel, pre-registered, -0.1 sigma"),
    ("K04 signs (lambda,w4,w6 > 0)", "traffic-core", "forced",
     "DRIFT K04 (l.95)", "stabiliser-loop minimality"),
    ("nu_R/dust RELATIVE abundance 4:1", "source", "forced",
     "DRIFT snapshot 2011, item123_r4_zero_mode_abundance_ratio.py", "directed-edge R4 incidence; the split is forced"),
    ("nuclear a_C (Coulomb)", "traffic-nuc", "forced",
     "DRIFT item113 (l.2071)", "(3/5) alpha hbar c / r0, r0=2a0, 2.2%"),
]

# ----------------------------------- tally -----------------------------------
def count(sector=None, status=None):
    return sum(1 for c, sec, st, *_ in LEDGER
               if (sector is None or sec == sector) and (status is None or st == status))

free_core = count("traffic-core", "free")
free_ew   = count("traffic-EW", "free")
free_nuc  = count("traffic-nuc", "free")
K_rel     = free_core + free_ew + free_nuc          # surviving FREE continuous traffic multipliers
cond_traffic = (count("traffic-core", "conditional") + count("traffic-EW", "conditional")
                + count("traffic-nuc", "conditional"))

K_clock   = 1                                        # Omega (conditionally pinned by the selector lock)
K_phase   = count("holonomy", "free")               # delta_nu
cond_holo = count("holonomy", "conditional")
K_Q       = count("constraint")                      # binary charge-existence premises (M15)
src_free  = count("source", "free")                 # source laws that are FITTED (would be multipliers)
src_cond  = count("source", "conditional")          # code-structured-but-open source laws (not multipliers yet)
K_traffic = K_clock + K_rel

# ----------------------------------- assertions ------------------------------
assert all(len(row) == 5 and row[3] for row in LEDGER), "every ledger entry must carry a citation"
assert count() == len(LEDGER) >= 20, "ledger must enumerate the open coefficient set (sanity floor)"
assert count("constraint", "premise") == K_Q and K_Q >= 1, "M15 is a BINARY charge premise (K_Q), not in K_rel"
assert src_free == 0, "no source law is currently a FITTED multiplier (alpha_0/208 is code-structured-conditional)"
assert free_ew + free_nuc > free_core, \
    "free multipliers are dominated by the EW+nuclear absolute-scale bolus, not the service core"
assert free_core == 0, "native coefficient ledger EXHAUSTED: zero free (w4/w6 + CC-loop + HBC A_nu + gravity alpha^1 all conditional; alpha^2/K_eff dead)"
assert cond_traffic >= 3, "several core coefficients are conditional (one premise away), not free"
assert K_clock == 1, "saturation collapses the absolute scale to ONE conditionally-pinned clock"

# ----------------------------------- report ----------------------------------
bar = "=" * 92
print(bar)
print("TRAFFIC-MULTIPLIER AUDIT (R1 parameter ledger) — counts from verified canon, 2026-06-19")
print(bar)
order = ["clock", "traffic-core", "traffic-EW", "traffic-nuc", "holonomy", "constraint", "source"]
mark = {"free": "FREE", "conditional": "cond", "forced": " ok ", "premise": "PREM"}
for sec in order:
    print(f"\n[{sec}]")
    for c, _, st, where, note in [r for r in LEDGER if r[1] == sec]:
        print(f"  {mark[st]}  {c}")
print(f"\n{bar}")
print("LEDGER TALLY  (constraint + source kept SEPARATE from the continuous-multiplier count)")
print(f"  K_clock  Omega                 : {K_clock}   (conditionally pinned: H0=Lambda/N_lock, selector lock)")
print(f"  K_rel    free traffic ratios   : {K_rel}  = core {free_core} + EW {free_ew} + nuclear {free_nuc}")
print(f"           (+ {cond_traffic} conditional: 10/27, latch value, item-126 T3)")
print(f"  K_phase  free holonomy         : {K_phase}   (delta_nu, geometric primitive) [+{cond_holo} cond: CKM sign]")
print(f"  K_Q      charge-existence      : {K_Q}   (BINARY premise: CMB dust; NOT a continuous knob)")
print(f"  K_source fitted source laws    : {src_free}   (+{src_cond} conditional: alpha_0/208, code-structured)")
print(f"  => K_traffic = K_clock + K_rel : {K_traffic}")
print(f"""
{bar}
VERDICT (audit, exit 0):
  With M15 and the abundance correctly removed from K_rel, the ledger SPLITS — and the split
  is the finding:

  * Saturation leg on track. Omega is not free — the cosmological-selector/service-span lock
    pins the H0 endpoint at conditional-theorem grade (-0.06%). "Saturated bandwidth fixes the
    one native clock" is REALISED for one cluster member.  The later two-anchor audit rejects
    the stronger claim that this one clock serves the EW/top second anchor.

  * Native cosmological/QEC core EXHAUSTED: {free_core} free multipliers (NONE remain) +
    {cond_traffic} conditional (w4/w6, 10/27, latch, T3, CC-loop, HBC A_nu, gravity alpha^1),
    against a dozen-plus forced wins
    (alpha, 3/14, c_4=alpha_0, s_1, Omega_Lambda, w_0, OZI, clock, 4:1, ...). K_rel(core) < N_core
    by ~2-3x — predictive-LEANING, margin-dependent on the conditionals discharging.

  * M15 is correctly OUT of K_rel: it is a BINARY Q_a-existence premise (K_Q=1), and its abundance
    is a code-structured-but-open source law (alpha_0/208), NOT a fitted multiplier. Neither
    inflates the continuous count. (If the dust abundance later proves fitted, src_free rises.)

  * The DANGER is localised and pre-existing: K_rel is DOMINATED by the EW ({free_ew}) + nuclear
    ({free_nuc}) stiffness bolus — a SECOND absolute scale (§16.2), not service-history relative
    activities. Counting will not reduce them; they need a second saturated anchor.

  NET: R1 is predictive on its native domain, does NOT worsen the matter-scale walls, and does not
  solve them. The native traffic-core free column is now EMPTY: alpha^2/K_eff is closed-negative /
  superseded by the alpha^1 erasure route, while the surviving single native phase frontier is
  delta_nu (holonomy). The continuous danger is the EW+nuclear second-scale bolus, plus one binary
  premise (M15 dust) that the framework flags as its most exposed sector.
{bar}""")
print(f"exit 0 -- ledger: K_clock={K_clock}, K_rel={K_rel} (core {free_core}), K_phase={K_phase}, "
      f"K_Q={K_Q} (binary), K_source(fitted)={src_free}; danger = the EW+nuclear second-scale bolus.")
