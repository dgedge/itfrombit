#!/usr/bin/env python3
r"""fermion_q4_base_channels.py -- promote the Q4 candidate's BASE d=2 from asserted to DERIVED.

The Q4 charge-neutrality candidate (DRIFT M9 cluster; fermion_q4_charge_selector.py) reads the
four-sector defect numerators as d = 2 (base) + 1 (on the smaller-|Q| doublet member). Promotion
condition (1) was: NAME the base d=2 on quarks as the firewall-free analogue of the lepton's two
nu_R channels (e_R->nu_R via I3-flip, nu_L->nu_R via chi-flip) -- i.e. an explicit count of the two
single-flip channels at the RH QUARK singlet, since the nu_R-Feshbach index is firewalled on LQ=1.

THIS SCRIPT enumerates the reduced 8-bit register (R2 chi=W imposed, generation fixed) and counts,
for each sector's RH singlet, the NON-COLOUR single-flip channels {LQ, I3, chi} that land on a
physical (R3,R4-allowed) state. The claim to verify:
  * lepton ref = nu_R (the R4-sterile up-type RH lepton): exactly 2 channels (e_R via I3, nu_L via chi);
  * quark refs = u_R, d_R (physical RH quark singlets):   exactly 2 channels each, and BOTH targets
    are quarks (LQ=1) -> FIREWALL-FREE (no nu_R / lepton involved);
and the MECHANISM that forces 2: the LQ-flip at any RH singlet is always R3-forbidden (LQ flips
lepton<->quark, which demands simultaneously toggling colour from 00<->coloured), so only the
I3-partner and chi-partner survive. Colour-flips are excluded by construction -- they are the
N_eff = N*N_c colour multiplicity of the candidate, not defect channels.

Convention note: I3=1 := up-type, I3=0 := down-type; chi=1 := RH, chi=0 := LH. The CHANNEL COUNT is
convention-independent (it depends on WHICH bit differs, not its value); only the labels move.

RESULT: base d=2 for all sectors, forced by R3; firewall-free for quarks. Promotion condition (1) MET.
Still open: condition (2) the §16.3 null-check on the |Q| +1-selector; the N_eff=27 colour-coherence
(item 96, separate); the down-R anomaly sqrt(12/5). So the candidate moves: base d=2 DERIVED, the
+1 selector remains the open promotion gate.
"""
from itertools import product

# reduced register bit names (R2 chi=W folded into chi; generation fixed)
# state = (LQ, C0, C1, I3, chi)
NONCOLOUR = {"LQ": 0, "I3": 3, "chi": 4}  # the three non-colour flip axes (indices into the tuple)


def coloured(s):
    return (s[1], s[2]) != (0, 0)


def passes_R3(s):
    # R3: LQ=0  <=>  colourless (C0,C1)=(0,0)
    return (s[0] == 0) == (not coloured(s))


def is_nuR(s):
    # the R4-sterile state: up-type (I3=1) RH (chi=1) colourless lepton (LQ=0, C=00)
    return s == (0, 0, 0, 1, 1)


def physical(s):
    # R1 trivial (generation fixed); R3; R4 removes nu_R from the ACTIVE codespace
    return passes_R3(s) and not is_nuR(s)


def name(s):
    LQ, C0, C1, I3, chi = s
    sect = "quark" if LQ == 1 else "lepton"
    up = "up" if I3 == 1 else "dn"
    hand = "R" if chi == 1 else "L"
    col = {(0, 1): "r", (1, 0): "g", (1, 1): "b", (0, 0): "-"}[(C0, C1)]
    if LQ == 0:
        base = {("up", "L"): "nu_L", ("up", "R"): "nu_R*", ("dn", "L"): "e_L", ("dn", "R"): "e_R"}[(up, hand)]
        return base
    base = {("up", "L"): "u_L", ("up", "R"): "u_R", ("dn", "L"): "d_L", ("dn", "R"): "d_R"}[(up, hand)]
    return f"{base}({col})"


def flip(s, axis):
    s = list(s)
    s[NONCOLOUR[axis]] ^= 1
    return tuple(s)


def base_channels(ref):
    """Non-colour single-flip channels from ref that land on a PHYSICAL state."""
    out = []
    for axis in ("LQ", "I3", "chi"):
        t = flip(ref, axis)
        out.append((axis, t, name(t), physical(t)))
    return out


def main():
    print("=== Q4 promotion: BASE d=2 = non-colour single-flip channels at the RH singlet ===\n")

    # sanity: enumerate the physical codespace (expect 15 = 16 - sterile nu_R)
    allstates = [s for s in product((0, 1), repeat=5)]
    phys = [s for s in allstates if physical(s)]
    nleps = sum(1 for s in phys if s[0] == 0)
    nqks = sum(1 for s in phys if s[0] == 1)
    print(f"  physical gen-1 states: {len(phys)}  (leptons {nleps}: nu_L,e_L,e_R; quarks {nqks}: u/d x L/R x 3col)")
    print(f"  sterile (R4-removed) reference: nu_R* = {name((0,0,0,1,1))} at (LQ,C0,C1,I3,chi)=(0,0,0,1,1)\n")

    refs = [
        ("lepton sector", (0, 0, 0, 1, 1), "nu_R (sterile up-type RH lepton)"),
        ("up sector",     (1, 0, 1, 1, 1), "u_R (red)"),
        ("down sector",   (1, 0, 1, 0, 1), "d_R (red)"),
    ]
    print(f"  {'sector':14s} {'RH singlet':10s}  channels (axis -> target : physical?)")
    base_ok = True
    firewall_ok = True
    for label, ref, desc in refs:
        chans = base_channels(ref)
        live = [c for c in chans if c[3]]
        d = len(live)
        base_ok &= (d == 2)
        # firewall: for quark refs, all live targets must be quarks (LQ=1)
        if ref[0] == 1:
            firewall_ok &= all(t[0] == 1 for _, t, _, ok in [(c[0], c[1], c[2], c[3]) for c in live] for ok in [True])
            firewall_ok &= all(c[1][0] == 1 for c in live)
        chan_str = "  ".join(f"{ax}->{nm}:{'Y' if ok else 'n(R3-forbid)'}" for ax, _, nm, ok in chans)
        print(f"  {label:14s} {name(ref):10s}  d={d}  [{chan_str}]")
    print()

    print("  [mechanism] why d is forced to 2 at every RH singlet:")
    print("    - LQ-flip: toggles lepton<->quark; but R3 ties LQ=0 <-> colourless, so flipping LQ alone")
    print("      makes a coloured lepton or a colourless quark -> R3-FORBIDDEN. Never a channel.")
    print("    - I3-flip: -> the isospin partner (up<->down) RH singlet. Always physical.")
    print("    - chi-flip: -> the LH partner (RH<->LH). Always physical.")
    print("    => exactly 2 surviving channels {I3-partner, chi-partner}: base d=2, FORCED by R3.")
    print("    Colour-flips are excluded (they are the N_eff=N*N_c colour multiplicity, not defects).")

    print("\n[verdict] BASE d=2 DERIVED (promotion condition 1 MET):")
    print("  Every sector's RH singlet has exactly TWO non-colour single-flip channels (the I3-partner")
    print("  and the chi-partner); the lepton<->quark (LQ) flip is always R3-forbidden. For leptons these")
    print("  are e_R->nu_R (I3) and nu_L->nu_R (chi) -- the canonical d=2. For quarks they are u_R<->d_R")
    print("  (I3) and u_R<->u_L / d_R<->d_L (chi) -- both QUARK-internal, so FIREWALL-FREE (no nu_R, the")
    print("  firewalled lepton index is never invoked). So the candidate's base d=2 is no longer asserted:")
    print("  it is forced by R3, identical in structure across leptons and quarks, quark-native.")
    print("  REMAINING for full Q4 closure: (2) the §16.3 null-check on the |Q| (+1) selector; the")
    print("  N_eff=27 colour-coherence (item 96, separate); the down-sector R=sqrt(12/5) anomaly.")

    assert len(phys) == 15 and nleps == 3 and nqks == 12, "codespace must be 15 = 3 leptons + 12 quarks"
    assert base_ok, "every RH singlet must have exactly 2 non-colour single-flip channels (base d=2)"
    assert firewall_ok, "quark-sector channels must be quark-internal (firewall-free)"
    # the LQ-flip must be the R3-forbidden one at each ref
    for _, ref, _ in refs:
        lqflip = flip(ref, "LQ")
        assert not physical(lqflip), "the LQ-flip at an RH singlet must be R3-forbidden"
    print("\nGATES PASSED -- 15-state codespace; base d=2 forced by R3 at every RH singlet; quark channels")
    print("firewall-free; the LQ-flip is the R3-forbidden one. Promotion condition (1) MET.")
    print("exit 0")


if __name__ == "__main__":
    main()
