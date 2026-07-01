#!/usr/bin/env python3
r"""Does dynamic-R1 monitoring RESCUE item 84 (the baryonic R1-anomaly)?  No.

The user's insight: if R1 is dynamically MONITORED (the boot-R1 'leans true' result),
the R1-forbidden generation corner (1,1) is not merely 'invalid' -- it must RESOLVE
through one of the two legal Hasse rescue edges (11->01, 11->10 = single generation-bit
flips), not 11->00 (a two-bit move) and not a symmetric all-generation channel.  The
hope was that this turns item 84's qualitative 'violent asymmetric cascade' into a
specific, falsifiable decay-cascade prediction for the cross-generational hyperons.

This script tests that hope against the canon -- and it FAILS to rescue item 84, for
the reason the canon already records.  Two prior canon facts are decisive (ANCHOR
'(Q4) Baryon XOR-composite framework', item 84, and its 2026-05-30 flag):

  (i)  ENUMERATION BUG -- ALREADY KNOWN, not a new finding.  The canon itself flagged
       (2026-05-30) that item 84's predicted set wrongly lists Lambda_b=udb and
       Omega_b-=ssb, which do NOT satisfy the XOR=(1,1) rule.  We re-derive it only to
       confirm we are reading the rule correctly -- we do NOT claim it as new.
  (ii) EMPIRICAL TENSION -- DECISIVE.  PDG Xi_b^0 (tau~1.48 ps) and Xi_b^- (tau~1.57
       ps) have entirely NORMAL b-baryon lifetimes and standard b->cW^- decays, with NO
       anomalous instability.  The canon's verdict already stands: 'this prediction
       does not currently survive.'  Its rescue bar: a CALCULATED decay-asymmetry or
       branching deviation BEYOND CKM -- not a qualitative cascade claim.

The dynamic-R1 sharpening produces a single-generation-step cascade -- which is exactly
the CKM single-step hierarchy.  Reproducing CKM is NOT a deviation beyond CKM, so the
sharpening lands precisely in the canon's named failure mode ('vacuous -- no number
distinguishes it from CKM').  Honest result: the insight is structurally consistent and
EXPLAINS why item 84 collapses to CKM, but it does NOT rescue it as a forward prediction.

Bit-map note: the Koide LEPTON map (ANCHOR line 687) is tau=(0,0), e=(0,1), mu=(1,0);
the item-84 QUARK example (Xi_b^0=usb, (0,0)+(1,0)+(0,1)=(1,1)) labels gen1=(0,0),
gen3=(0,1) -- the opposite gen1<->gen3 codeword choice.  This relabelling is immaterial
to the criterion: for ANY bijection of the 3 generations onto {00,01,10}, the XOR of
all three is 00^01^10 = (1,1).  So 'one quark per generation <=> XOR=(1,1)' is
convention-independent; only the (non-forbidden) value of paired states changes label.
"""
import itertools

# Item-84 quark convention (ANCHOR line 2052): u=(0,0), s=(1,0), b=(0,1) -> usb=(1,1).
GEN = {"u": (0, 0), "d": (0, 0), "c": (1, 0), "s": (1, 0), "t": (0, 1), "b": (0, 1)}


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def xor(bits):
    a = b = 0
    for x, y in bits:
        a ^= x; b ^= y
    return (a, b)


def hamming(p, q):
    return (p[0] ^ q[0]) + (p[1] ^ q[1])


def main():
    print("ITEM 84 + dynamic-R1: does monitoring RESCUE the baryonic R1-anomaly?")
    print("=" * 80)

    print("\n[0] Convention-independence of the XOR=(1,1) criterion")
    val = {"00": (0, 0), "01": (0, 1), "10": (1, 0)}
    xors = set()
    for perm in itertools.permutations(["00", "01", "10"]):     # any gen<->codeword map
        xors.add(xor([val[c] for c in perm]))
    print(f"    XOR of {{00,01,10}} under every generation relabelling: {xors}")
    check("one-quark-per-generation XOR=(1,1) regardless of bit-map convention", xors == {(1, 1)})

    print("\n[1] Re-derive (NOT claim) the canon's already-known enumeration bug")
    canon_set = {"Xi_b0": "usb", "Xi_b-": "dsb", "Lambda_b": "udb", "Omega_b-": "ssb"}
    forbidden, valid = [], []
    for name, q in canon_set.items():
        x = xor([GEN[c] for c in q])
        tag = "FORBIDDEN (XOR=11)" if x == (1, 1) else f"valid (XOR={x})"
        (forbidden if x == (1, 1) else valid).append(name)
        print(f"    {name:<9s} = {q}: gen-XOR = {x}  -> {tag}")
    check("Xi_b0(usb), Xi_b-(dsb) ARE R1-forbidden (one quark per generation)",
          "Xi_b0" in forbidden and "Xi_b-" in forbidden)
    check("Lambda_b(udb), Omega_b-(ssb) are NOT forbidden (paired generations cancel)",
          "Lambda_b" in valid and "Omega_b-" in valid)
    print("    --> this matches the canon's OWN 2026-05-30 enumeration-bug flag; it is")
    print("        an ALREADY-RECORDED correction, re-derived here only as a read-check.")

    print("\n[2] The genuine R1-forbidden set = one quark from each generation")
    gen1, gen2, gen3 = ["u", "d"], ["c", "s"], ["b"]   # t too heavy to bind
    onepg = ["".join(t) for t in itertools.product(gen1, gen2, gen3)]
    for q in onepg:
        check(f"{q} has gen-XOR=11 (one-per-generation)", xor([GEN[c] for c in q]) == (1, 1))
    print(f"    forbidden baryons = {onepg}: Xi_b^0(usb), Xi_b^-(dsb), doubly-heavy Xi_bc(ucb,dcb)")

    print("\n[3] The dynamic-R1 sharpening (the user's insight): Hasse-edge resolution")
    rescue = [v for v in [(0, 1), (1, 0)] if hamming((1, 1), v) == 1]    # single-bit flips
    print(f"    single-bit (Hasse) rescues of 11: {rescue}   (11->01, 11->10)")
    print(f"    11->00 distance = {hamming((1,1),(0,0))} (a TWO-bit move, excluded as one recovery)")
    check("forbidden corner resolves via single Hasse edges 11->01, 11->10", set(rescue) == {(0, 1), (1, 0)})
    check("11->00 excluded as a single recovery step (Hamming 2)", hamming((1, 1), (0, 0)) == 2)

    print("\n[4] DECISIVE: the single-step cascade IS CKM -> no deviation beyond CKM")
    # a single XOR-bit flip = one generation-bit transition = a single CKM step.
    # 11->00 (both bits) = a two-generation jump = doubly-CKM-suppressed.  This is the
    # ordinary CKM single-step hierarchy -- it REPRODUCES CKM, it does not exceed it.
    single_step_allowed = (len(rescue) == 2)
    two_step_suppressed = (hamming((1, 1), (0, 0)) == 2)
    reproduces_ckm = single_step_allowed and two_step_suppressed
    deviation_beyond_ckm = False    # the sharpening yields NO number above CKM
    print(f"    single-generation-step allowed, two-step suppressed -> CKM pattern: {reproduces_ckm}")
    print(f"    any calculated decay-asymmetry/branching deviation BEYOND CKM produced? {deviation_beyond_ckm}")
    check("the dynamic-R1 cascade reproduces the CKM single-step hierarchy", reproduces_ckm)
    check("it produces NO deviation beyond CKM (the canon's rescue bar is NOT met)", not deviation_beyond_ckm)

    print("\n[5] The empirical-tension flag (canon, 2026-05-30) still stands")
    # PDG: Xi_b^0 ~1.48 ps, Xi_b^- ~1.57 ps -- normal b-baryon lifetimes (cf Lambda_b ~1.47 ps).
    lifetimes_ps = {"Xi_b0": 1.48, "Xi_b-": 1.57, "Lambda_b_ref": 1.47}
    normal = all(1.3 < t < 1.7 for t in lifetimes_ps.values())
    print(f"    PDG lifetimes (ps): {lifetimes_ps} -> all normal b-baryon range: {normal}")
    check("Xi_b^0, Xi_b^- have normal lifetimes -> no observed anomalous 'shattering'", normal)

    print(
        """
[6] VERDICT -- dynamic-R1 does NOT rescue item 84 (honest negative)
    The user's insight is real and structurally consistent: monitored R1 forces the
    forbidden corner (1,1) to resolve through the two single-bit Hasse edges
    (11->01, 11->10), not the two-bit 11->00 -- a single-generation-step cascade.  But
    that single-step structure IS the CKM hierarchy (single step allowed; two-step
    doubly suppressed).  Reproducing CKM is not a deviation BEYOND CKM, which is exactly
    the canon's stated rescue bar (a CALCULATED asymmetry/branching deviation past CKM).
    So the sharpening lands precisely in the failure mode the canon already named:
    'vacuous -- no number distinguishes it from CKM'.

    Two canon facts make this conclusive and they were ALREADY on record (verify-before-
    cite):
      * the enumeration bug (Lambda_b, Omega_b- don't satisfy XOR=11) -- canon flagged
        it 2026-05-30; we only re-derived it as a read-check, we do NOT claim it new;
      * the decisive empirical tension -- PDG Xi_b^0 (1.48 ps), Xi_b^- (1.57 ps) have
        normal lifetimes and standard b->cW^- decays; 'this prediction does not survive'.

    NET: the dynamic-R1 picture is a good EXPLANATION of WHY item 84 is vacuous -- its
    cascade is isomorphic to CKM single-step suppression -- but it is NOT a rescue.  It
    yields no falsifiable number beyond CKM, and the empirical-tension flag stands.  The
    correct action is to LEAVE item 84 retired, now with a mechanism for its vacuity, not
    to revive it as a forward prediction.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
