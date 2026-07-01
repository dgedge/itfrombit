#!/usr/bin/env python3
r"""ITEM 144 — close the flagged gap: tie the arrow-of-time entropy-production RATE to the
substrate error rate, replacing the toy placeholder C_REC=0.5 bits/cell/tick in
`entropy_arrow_monotonicity.py` with a DERIVED value.

The toy (entropy_arrow_monotonicity.py) shows S(0)=0 and dS/dt>=0, but leaves the *rate*
C_REC as a free knob. Both factors that fix it are already canon:

  rate of record-writing events   = the per-cell-per-tick leak/error rate  = alpha   [item 79]
                                    (alpha = Tr_non-unit[W|K2]; the fine-structure constant
                                     in its substrate role as the QEC leak probability "p")
  information per record event     = log2(address alphabet)                = log2(8) = 3 bits
                                    [record_content_from_syndrome.py, Lemma 1: a single-site
                                     flip toggles EXACTLY the 3 checks incident to its vertex,
                                     and the 8 incident triples are 8 distinct 12-bit syndromes
                                     -> the increment is a vertex label, alphabet 8, DERIVED,
                                     verified over all 256x8 flips]

  =>  dS/dt | per active cell  =  alpha * log2(8)  =  3 alpha   bits / cell / tick.

So the QEC error probability the gap asks for IS alpha (item 79) = the fine-structure constant,
and the arrow of time advances at 3*alpha bits/cell/tick.

TWO LEDGERS (the distinction that resolves the apparent rate ambiguity):
  * RECORD entropy S  -- STOCHASTIC syndrome increments at rate alpha, each a 3-bit address
    record. This accumulates as the irreversible archive = the arrow of time (item 144).
  * WASTE HEAT  rho_Lambda -- the DETERMINISTIC every-tick R4 housekeeping erasure (the I3=1
    parity flip S:t->~t of section 5.2, erased every tick). A deterministic flip carries no
    Shannon uncertainty, so it produces Landauer HEAT (-> dark-energy exhaust, item 123), not
    archive entropy. Different rate (~1/tick), different ledger.
The arrow-of-time rate is therefore the ERROR rate alpha, not the housekeeping rate.

This script DERIVES the rate and self-asserts the arithmetic; it does not re-run the toy
(monotonicity is independent of the value of C_REC). Honest scope is printed in the verdict.
"""
import math

# ---- canon inputs (cited, not re-derived here) ----
ALPHA      = 1.0 / 137.036          # item 79: substrate leak rate per cell per tick == fine-structure const
A_ADDRESS  = 8                      # record_content_from_syndrome.py Lemma 1: syndrome-increment address alphabet
C_REC_TOY  = 0.5                    # entropy_arrow_monotonicity.py: the FREE placeholder we are replacing

# ---- the derived intensive entropy-production rate ----
bits_per_event = math.log2(A_ADDRESS)               # = 3 exactly
dS_dt_bits     = ALPHA * bits_per_event              # bits / active-cell / tick   (= 3 alpha)
dS_dt_nats     = ALPHA * math.log(A_ADDRESS)         # k_B (nats) / active-cell / tick (= alpha*ln8)

# ---- the shared-alpha thread (one substrate constant, three record-channel outputs) ----
eta_baryon   = (3.0/14.0) * ALPHA**4                 # item 126 magnitude: baryon-to-photon ~ alpha^4
eta_over_arrow = eta_baryon / dS_dt_bits             # = (1/14) alpha^3  -- pure number, shared-alpha by-product
expect_ratio = (1.0/14.0) * ALPHA**3

print("="*86)
print("ITEM 144 — entropy-production rate from the substrate error rate (closing the C_REC gap)")
print("="*86)
print(f"  leak / record-event rate           alpha            = {ALPHA:.6e}  per cell per tick  [item 79]")
print(f"  record content per increment       log2(8)          = {bits_per_event:.3f} bits           [record_content_from_syndrome.py L1]")
print(f"  --------------------------------------------------------------------------------------")
print(f"  DERIVED  dS/dt | per active cell  = alpha*log2(8)   = {dS_dt_bits:.6f} bits/cell/tick (= 3 alpha)")
print(f"                                    = alpha*ln 8      = {dS_dt_nats:.6f} k_B /cell/tick")
print(f"  toy placeholder (now superseded)   C_REC            = {C_REC_TOY} bits/cell/tick  (free knob -> {C_REC_TOY/dS_dt_bits:.0f}x too large)")
print()
print(f"  shared-alpha thread:  arrow = 3*alpha = {dS_dt_bits:.4e} bits/cell/tick")
print(f"                        baryon = (3/14)alpha^4 = {eta_baryon:.4e}   (item 126 magnitude)")
print(f"                        eta / (dS/dt) = (1/14)alpha^3 = {eta_over_arrow:.4e}  (one alpha, two outputs)")

# ---------------------------------- self-assertions ----------------------------------
assert abs(bits_per_event - 3.0) < 1e-12,                       "address alphabet 8 must give exactly 3 bits"
assert abs(dS_dt_bits - 3.0*ALPHA) < 1e-15,                     "dS/dt must equal 3*alpha bits/cell/tick"
assert abs(dS_dt_nats - dS_dt_bits*math.log(2)) < 1e-15,        "bits<->nats conversion must be exact (3 bits = ln8 nats)"
assert 0.0 < dS_dt_bits < C_REC_TOY,                           "derived rate must be positive and below the toy placeholder"
assert abs(eta_over_arrow - expect_ratio) < 1e-30,             "eta/(dS/dt) must reduce to (1/14)alpha^3 (shared-alpha)"
assert 6.0e-10 < eta_baryon < 6.2e-10,                         "baryon magnitude sanity (item 126, ~6.08e-10)"

print(f"""
{"="*86}
VERDICT (item 144 rate, DERIVED):
  The arrow-of-time entropy-production rate is no longer a free knob:
     dS/dt = alpha * log2(8) = 3*alpha = {dS_dt_bits:.5f} bits per active cell per tick.
  Both factors are canon: the event rate is the leak rate alpha (item 79 = fine-structure
  constant), and the per-event content is the DERIVED address alphabet 8 (=3 bits,
  record_content_from_syndrome.py L1). The same alpha threads the baryon asymmetry
  ((3/14)alpha^4, item 126) and the dark-energy exhaust (item 123) -- one substrate constant,
  several record-channel outputs.

  TIER: derived-conditional (intensive rate). Solid: the FORM (rate x content), alpha [item 79],
  the alphabet-8 content [verified all 256x8]. OPEN/HONEST SCOPE:
    * the EXTENSIVE cosmic budget dS/dt_total = 3*alpha * N_active(t) needs the crystallisation
      / expansion history N_active(t) -- NOT supplied here;
    * 3*alpha is the FULL-record value; sub-threshold corrections that heal without committing
      an address record give an effective coefficient <= 3, so dS/dt = c*alpha with c in (0,3];
    * this fixes the RATE only; monotonicity / S(0)=0 are unchanged (entropy_arrow_monotonicity.py).
  This is the RECORD ledger (rate alpha); the deterministic every-tick housekeeping erasure is a
  SEPARATE ledger (Landauer heat -> rho_Lambda, item 123), not archive entropy.
{"="*86}""")
print("exit 0 -- dS/dt = 3*alpha bits/cell/tick DERIVED from canon (item 79 rate x alphabet-8 content); C_REC gap closed.")
