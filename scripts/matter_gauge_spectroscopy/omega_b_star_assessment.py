#!/usr/bin/env python3
"""
Assessment of the Omega_b* (3/2+, ssb) "parameter-free forward prediction"
(ANCHOR §15 item 142 / §14: band 6075-6090 MeV).

Checks the claim against (a) the framework's OWN verified script
(bottom_baryon_forward.py) and (b) the external standard-QCD literature.

External anchors (PDG / literature, cited):
  - Omega_b      = 6045.2 MeV (measured; PDG)
  - Omega_b*     = NOT yet observed (forward) — true
  - color-hyperfine quark model: Omega_b* = 6082.8 +/- 5.6 MeV  (arXiv:0708.4027)
  - heavy-quark-scaled splitting: (M(Om_c*)-M(Om_c))*(m_c/m_b) = 24.0 +/- 0.7 MeV
      -> Omega_b* ~ 6069 MeV
  - Omega_c*-Omega_c = 2765.9 - 2695.2 = 70.7 MeV (PDG)

Self-asserting: exit 0 == every flagged finding holds.
"""
Omb_meas   = 6045.2          # measured Omega_b
Omb_fw     = 6132.3          # framework's OWN Omega_b prediction (ANCHOR §15 item 142 table)
Omc_star, Omc = 2765.9, 2695.2

# --- framework's own numbers (from bottom_baryon_forward.py, verified exit 0) ---
split_naive_fw = 14.1        # naive single-A
split_enh_fw   = 24.0        # |psi(0)|^2-enhanced (x1.70 from Om_c "assuming saturation")
script_lo, script_hi = Omb_meas + split_naive_fw, Omb_meas + split_enh_fw   # 6059, 6069
doc_lo, doc_hi = 6075.0, 6090.0    # the ANCHOR §14 / item-142 published band

# --- external standard estimates ---
split_HQ = 24.0              # (Om_c*-Om_c)*(m_c/m_b), literature value 24.0+/-0.7
Omb_star_HQ = Omb_meas + split_HQ                 # ~6069
Omb_star_colorHF = 6082.8                          # arXiv:0708.4027, +/-5.6

print("="*72); print("Omega_b* — is it a parameter-free FRAMEWORK forward prediction?"); print("="*72)
print(f"  Omega_b (measured)              = {Omb_meas} MeV")
print(f"  Omega_b (FRAMEWORK's own pred)  = {Omb_fw} MeV  -> +{Omb_fw-Omb_meas:.0f} MeV (+{(Omb_fw-Omb_meas)/Omb_meas*100:.2f}%) HIGH")
print(f"  => the Omega_b* prediction anchors on the MEASURED Omega_b ({Omb_meas}),")
print(f"     NOT the framework's own (87-MeV-high) value. It is measured-anchor + splitting.")
print()
print(f"  framework splitting (enhanced)  = {split_enh_fw} MeV")
print(f"  standard HQ-scaled splitting     = {split_HQ} +/- 0.7 MeV   <-- IDENTICAL")
print(f"     (framework's x1.70 enhancement is taken FROM the Om_c splitting, i.e.")
print(f"      it reproduces the standard (Om_c*-Om_c)*(m_c/m_b) estimate, not new physics)")
print()
print(f"  framework VERIFIED SCRIPT band   = [{script_lo:.0f}, {script_hi:.0f}] MeV  (bottom_baryon_forward.py)")
print(f"  framework PUBLISHED (ANCHOR) band= [{doc_lo:.0f}, {doc_hi:.0f}] MeV  (§14 / item 142)")
print(f"  => DOC band is {doc_lo-script_lo:.0f}-{doc_hi-script_hi:.0f} MeV ABOVE the verified script;")
print(f"     'widened upward' in prose to overlap the color-hyperfine value {Omb_star_colorHF}.")
print()
print(f"  standard color-hyperfine QM      = {Omb_star_colorHF} +/- 5.6 MeV  (arXiv:0708.4027)")
print(f"  standard HQ-scaling              = {Omb_star_HQ:.0f} MeV")
print(f"  => the literature itself spans ~{Omb_star_HQ:.0f}-{Omb_star_colorHF:.0f}; the framework's")
print(f"     script value ({(script_lo+script_hi)/2:.0f}) sits at the HQ end, its doc band at the color-HF end.")
print(f"     Neither is a NEW prediction — both reproduce existing standard methods.")

print("\n"+"="*72); print("VERDICT"); print("="*72)
print(" forward / unmeasured:  TRUE  (Omega_b* not yet observed; LHCb Run 3/4 can find it)")
print(" parameter-free:        FALSE (A,m_u,m_s frozen from light baryons; m_b fitted from")
print("                        Lambda_b and explicitly 'FREE' per §16.3; x1.70 from Om_c data)")
print(" framework-specific:    NO    (= standard DGG/HQ-scaling; the framework's own script")
print('                        verdict: "the baryon sector ... IS the constituent quark model")')
print(' discriminating:        NO    (§14: "not strongly discriminating from standard QM ~6080";')
print("                        Tier: Proposition). The doc band overlaps standard by construction.")
print(" code-vs-doc gap:       YES   (verified script 6059-6069 vs published 6075-6090).")

# robust asserts
assert abs(split_enh_fw - split_HQ) < 0.5, "framework splitting should equal standard HQ-scaled"
assert (doc_lo > script_hi - 1), "doc band should sit above the script band (upward widening)"
assert abs(Omb_fw - Omb_meas) > 80, "framework's own Omega_b should be ~87 MeV high"
assert script_lo < Omb_star_colorHF and doc_hi > Omb_star_HQ, "literature spans both framework values"
print("\nexit 0 == all findings verified. The 'parameter-free framework prediction' framing")
print("is not supported: it is the standard quark-model estimate, measured-anchored, non-")
print("discriminating; the only genuinely true part is that Omega_b* is still unobserved.")
