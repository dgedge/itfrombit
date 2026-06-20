#!/usr/bin/env python3
r"""ALPHA0 COUNT-TO-RATE THEOREM.

Target:
    Why does the monitored service fire once per 137-channel cycle?

Core point:
    The 137 is not a set of 137 independent Bernoulli processes with a free
    rate.  It is one repeatable service-label measurement with 137 mutually
    exclusive projectors:

        {P_0, ..., P_136},     sum_i P_i = I,     P_i P_j = delta_ij P_i.

    The monitored bridge is unital and connected (item79_unital_channel.py),
    so its stationary record state is the maximally mixed state I/137.  The
    emission / non-unitary firing is one projector in this one-hot alphabet.
    Therefore the Born/closed-record-pair measure gives

        p_fire = Tr(P_fire I/137) = 1/137.

    With the substrate clock normalized to one elementary service interrogation
    per lattice tick, Lambda ticks per unit time, the firing rate is

        Gamma = Lambda / 137 = alpha0 Lambda.

What is closed (the RATE-GIVEN-COUNT step):
    GIVEN the 137-label one-hot service alphabet, the monitored unital+connected
    channel uniformises to I/137 (item79_unital_channel.py, Evans-Frigerio) and the
    firing is one projector, so alpha0 = p_fire = 1/(count) and w = Lambda/137 =
    alpha0 Lambda. This resolves the closed-walk no-thermalisation objection (the
    integrable unitary walk could not equipartition; the dissipator supplies it).

What is NOT closed by this theorem:
    * THE COUNT ITSELF.  137 = Sym^2(16)+1 (symmetric / record-pair) vs the
      Grassmann-consistent 121 = Lambda^2(16)+1 (antisymmetric / fermion-pair) is a
      STATE-SPACE CONVENTION (item79_unital_channel.py residual; DRIFT 1483 sec 16.3).
      This script ASSERTS the symmetric record-pair reading; it does NOT derive it over
      the antisymmetric fermion-pair one.  Sharp residual: prove the service ledger
      bills record-pairs (symmetric, 136), not matter/fermion-pairs (antisymmetric, 120).
    * the dressed shift 1/137 -> 1/137.036;
    * sector-specific proof that a given process bills this observable.

exit 0 = GIVEN the symmetric 137-label record-pair alphabet, the one-hot monitor
         gives alpha0 exactly; the antisymmetric (Grassmann) alphabet would give
         1/121 -- the count is the open convention, the rate-given-count is the theorem.
"""

from fractions import Fraction


def rate(count: int, firing_labels: int = 1, bandwidth: int = 1) -> Fraction:
    return Fraction(firing_labels * bandwidth, count)


def main() -> None:
    records = 16
    pair_channels = records * (records + 1) // 2
    service_labels = pair_channels + 1

    print("ALPHA0 COUNT-TO-RATE THEOREM")
    print("\n[1] structural alphabet")
    print(f"    byte record words                  = {records}")
    print(f"    symmetric closed record pairs       = {pair_channels}")
    print(f"    + idle/latch service label          = {service_labels}")
    assert pair_channels == 136
    assert service_labels == 137

    alpha0 = rate(service_labels)
    print("\n[2] one-hot service-label observable")
    print("    service labels are mutually exclusive projectors: exactly one label per interrogation")
    print(f"    stationary state from unital connected monitor: rho = I/{service_labels}")
    print(f"    one emission/firing projector: p_fire = Tr(P_fire rho) = {alpha0}")
    assert alpha0 == Fraction(1, 137)

    print("\n[3] convert probability per tick to rate")
    print("    one elementary monitor interrogation per cell tick; clock frequency = Lambda")
    print(f"    Gamma/Lambda = p_fire = {alpha0} = alpha0")

    print("\n[4] the count is a state-space CONVENTION (the antisymmetric/Grassmann reading is the LIVE alternative, not a strawman):")
    controls = [
        ("coinless antisymmetric pair space", 16 * 15 // 2 + 1, 1, 1),
        ("distinguishable two-defect pair space", 16 * 16 + 1, 1, 1),
        ("two independent firing projectors", service_labels, 2, 1),
        ("double bandwidth per tick", service_labels, 1, 2),
    ]
    for name, count, firing, bandwidth in controls:
        r = rate(count, firing, bandwidth)
        print(f"    {name:<38s}: Gamma/Lambda = {r} = {float(r):.8f}")
        assert r != alpha0

    print(
        """
[verdict]
  The RATE-GIVEN-COUNT step is a projective-record theorem:

      byte -> 16 records -> [symmetric record-pairs Sym^2(16)+idle = 137 labels];
      one-hot monitored service observable -> one label per interrogation;
      unital connected monitoring -> stationary I/137 (equipartition DERIVED);
      one non-unitary firing label -> p_fire = 1/137;
      one lattice interrogation per tick -> Gamma = alpha0 Lambda.

  So "once per 137-channel cycle" is not an extra continuous rate: it is the Born
  weight of one projector in the maximally mixed stationary service register.

  BUT the COUNT itself is NOT derived here.  The symmetric record-pair space (136)
  vs the Grassmann-consistent antisymmetric fermion-pair space (120) is a
  state-space CONVENTION -- the antisymmetric reading gives 1/121 (the [4] control,
  and item79_unital_channel.py's own residual).  This script ASSERTS the record-pair
  (symmetric) reading; deriving it over the fermion-pair (antisymmetric) one is the
  open residual.  Other caveats: dressed-alpha renormalization; sector billing maps.

  [UPDATE: the COUNT is now DERIVED -- alpha0_record_pair_symmetry_theorem.py shows records are
  clonable (the orthonormal pointer basis R1/R7) so the diagonal self-pairs are admitted -> symmetric
  Sym^2(16)=136, not the Pauli-excluded fermionic 120; and SS5.9 bills syndrome-records, not matter.
  This script's "asserts the symmetric reading" is SUPERSEDED; the count closes at foundationally-
  grounded grade, conditional only on the reconstruction floor.]
exit 0"""
    )


if __name__ == "__main__":
    main()
