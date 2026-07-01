#!/usr/bin/env python3
r"""Item 87 -- the user's "crack in the wall": move the closed R1 edge-pair record
from the R4 sterile repair (electroweak, wrong register -> wall) to a BOOT
generation-register stabilisation.

Why it could crack the wall: the R4 wall was register-specific. R4 reads
{LQ,I_3,chi} (electroweak) and so can never see generations. But R1 -- the
generation rule G0*G1 != 1 -- reads the GENERATION register {G0,G1}. So a
generation-edge record must come from an R1/boot process, not R4. And the valid
R1-adjacency among the three allowed generations is exactly the path that gives
B B^T.

The honest catch (from canon): R1 is documented as a STATIC selection rule
(ANCHOR 359/378: "truncates to 3 valid generations", "excludes the 4th
generation"), and the boot SOURCE is the S3-singlet (line 745). So a *dynamical*
boot generation-recovery with a Stinespring environment is NOT documented -- it is
the new primitive the proposal must posit.  This script separates the part that is
right (register + geometry) from the part that is new physics (boot-as-recovery).
"""
import numpy as np

# the documented checks and which register each READS
R4_READS = {"LQ", "I_3", "chi"}     # electroweak/lepton (line 790) -> the R4 wall
R1_READS = {"G0", "G1"}             # generation register (G0*G1 != 1) -> the crack
GEN_BITS = {"G0", "G1"}

VALID = {"tau": (0, 0), "e": (0, 1), "mu": (1, 0)}   # 3 allowed generations
FORBIDDEN = (1, 1)                                   # excluded by R1


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def hamming(a, b):
    return (a[0] ^ b[0]) + (a[1] ^ b[1])


def eplane_ellipticity(C):
    u1 = np.array([1.0, -1.0, 0.0]) / np.sqrt(2)
    u2 = np.array([1.0, 1.0, -2.0]) / np.sqrt(6)
    U = np.column_stack([u1, u2])
    lo, hi = np.sort(np.linalg.eigvalsh(U.T @ C @ U))
    return 0.0 if abs(hi + lo) < 1e-15 else (hi - lo) / (hi + lo)


def main():
    print("ITEM 87 -- boot/R1 crack: right register, correct geometry, one open premise")
    print("=" * 82)

    print("\n[1] THE CRACK: R1 reads the GENERATION register; R4 reads electroweak")
    print(f"    R4 reads {sorted(R4_READS)} (electroweak) -> cannot see generations (the wall)")
    print(f"    R1 reads {sorted(R1_READS)} (generation)   -> CAN host a generation-edge record")
    check("R4 cannot access the generation register (disjoint) -- the documented wall", R4_READS.isdisjoint(GEN_BITS))
    check("R1 DOES read the generation register -- the crack: right channel for the record", R1_READS == GEN_BITS)

    print("\n[2] GEOMETRY: valid R1-adjacency among {tau,e,mu} is the path -> B B^T")
    names = ["e", "tau", "mu"]  # order rows as e,tau,mu
    idx = {n: i for i, n in enumerate(names)}
    # build the valid-adjacency incidence: edges = distance-1 pairs of VALID generations
    edges = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if hamming(VALID[a], VALID[b]) == 1:
                edges.append((a, b))
    print(f"    valid distance-1 edges among generations: {edges}")
    B = np.zeros((3, len(edges)))
    for j, (a, b) in enumerate(edges):
        B[idx[a], j] = 1; B[idx[b], j] = 1
    K = B @ B.T
    print(f"    B B^T =\n{K.astype(int)}")
    check("valid R1-adjacency is the path e-tau-mu (tau central, 2 edges)", edges == [("e", "tau"), ("tau", "mu")])
    check("its covariance is exactly K_R1 = B B^T", np.allclose(K, [[1, 1, 0], [1, 2, 1], [0, 1, 1]]))
    check("E-plane ellipticity 1/2 (-> d/N with the isotropic completion) = delta", abs(eplane_ellipticity(K) - 0.5) < 1e-9)
    # the forbidden-corner edges {11-e, 11-mu} are a DIFFERENT object (error-correction, not valid-adjacency)
    fc = [n for n in names if hamming(VALID[n], FORBIDDEN) == 1]
    check("the proposal needs VALID-adjacency (path), not the forbidden-corner edges (e,mu only)", set(fc) == {"e", "mu"})

    print("\n[3] THE CATCH: R1 is documented STATIC; no dynamical boot generation-recovery")
    print("    ANCHOR 359/378: R1 'truncates to 3 valid generations' / 'excludes the 4th' --")
    print("    a static codespace selection, NOT a recovery with a Stinespring environment.")
    print("    Line 745: the boot SOURCE is the S3-singlet (generation-blind).")
    r1_is_static_selection = True       # documented (lines 359/378)
    boot_source_generation_blind = True # documented (line 745)
    check("R1 is documented as a static selection rule, not a dynamical recovery", r1_is_static_selection)
    check("the boot source is the S3-singlet (generation-blind)", boot_source_generation_blind)

    print(
        """
[4] VERDICT -- a genuine crack, not the R4 wall, but one new premise remains
    The crack is real and well-aimed: unlike R4 (electroweak, register-fatal), the
    R1/generation channel reads the RIGHT register, and the valid R1-adjacency among
    the three allowed generations is EXACTLY the path e-tau-mu -> B B^T (delta) with
    the oriented part giving Phi. So the geometry and the register both work, and the
    "wrong register" wall does NOT apply here.

    The remaining premise is softer than the R4 wall but still NEW PHYSICS: canon
    documents R1 as a STATIC selection rule (truncation), and the boot source as the
    S3-singlet -- so there is no documented *dynamical* boot generation-register
    recovery with a Stinespring environment to write the edge records. The proposal
    must posit exactly that:

        booting the R1-valid generation register is a RECOVERY whose environment
        records each defect's valid (distance-1) R1-adjacency -- closed + oriented --
        while the traced physical channel stays generation-blind.

    This is the best-placed candidate so far: right register, correct geometry,
    unifies delta + Phi, and NOT register-walled.  It is not closed -- the single
    open premise is whether generation-register stabilisation at boot is dynamical
    (a recovery) rather than a static codespace cut.  If the framework's boot is an
    error-correcting establishment of the generation code (plausible -- booting IS a
    process), the record geometry is forced and delta=d/N + Phi-sigma close together.
    If the generation register is only ever a static selection (no boot recovery
    event), the record has no home and delta stays phenomenological.
    The wall is cracked from "wrong register (fatal)" to "is boot dynamical? (open)".
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
