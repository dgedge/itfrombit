#!/usr/bin/env python3
r"""Item 87 -- close the smaller-|Q| / sterility selector.

Question
--------
The Dynamic-R1 contact-count gate reduced the Koide-magnitude theorem to one
remaining selector:

    why does the extra +1 defect contact sit on {nu, d}, not {nu, u}?

The old frozen-I3 rule selects {nu, u}.  It is the wrong object for quarks
because it tracks the up-type reference.  The candidate that matches the data is
"smaller |Q|" / sterility, but as a naked value rule that was only a proposition.

Closure attempted here
----------------------
Use the post-EWSB monitored charge algebra.  The two electroweak records
(T3, Y) are not separately observed at long distance; the unbroken monitor reads

    Q = T3 + Y.

The closed-pair service credit for fusing the two records into the single
unbroken U(1)_EM record is the reduction in quadratic monitored load:

    C_close = (T3^2 + Y^2) - (T3 + Y)^2 = -2 T3 Y.

Thus a binary closed cancellation contact exists exactly when T3 and Y have
opposite signs.  In each weak doublet that is also exactly the member with the
smaller |Q|:

    lepton doublet:  nu  has T3=+1/2, Y=-1/2  -> C_close>0;  e does not.
    quark doublet:   d  has T3=-1/2, Y=+1/6  -> C_close>0;  u does not.

The extra +1 is therefore not a continuous charge-weighted coefficient.  It is
the presence/absence of one closed electroweak cancellation-pair record.  That
is the record-action/QEC type needed by the R1 contact ledger.

Verdict
-------
The selector closes at finite algebra grade, conditional only on the already
adopted premises that (i) electric charge is the unbroken monitored generator
Q=T3+Y, and (ii) contact counts are closed-record-pair counts, not continuous
charge weights.  Under those premises the {nu,d} selector is forced.
"""

from __future__ import annotations

from fractions import Fraction as F


def check(name: str, cond: bool) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


MEMBERS = {
    "nu": dict(doublet="lepton", T3=F(1, 2), Y=F(-1, 2), target_delta=F(3, 9), coloured=False),
    "e":  dict(doublet="lepton", T3=F(-1, 2), Y=F(-1, 2), target_delta=F(2, 9), coloured=False),
    "u":  dict(doublet="quark",  T3=F(1, 2), Y=F(1, 6),  target_delta=F(2, 27), coloured=True),
    "d":  dict(doublet="quark",  T3=F(-1, 2), Y=F(1, 6), target_delta=F(1, 9), coloured=True),
}

DOUBLETS = {"lepton": ("nu", "e"), "quark": ("u", "d")}


def q(m: str) -> F:
    return MEMBERS[m]["T3"] + MEMBERS[m]["Y"]


def closed_pair_credit(m: str) -> F:
    t3, y = MEMBERS[m]["T3"], MEMBERS[m]["Y"]
    return t3 * t3 + y * y - (t3 + y) * (t3 + y)


def plus_one_closed_pair(m: str) -> int:
    return int(closed_pair_credit(m) > 0)


def plus_one_frozen_i3(m: str) -> int:
    return int(MEMBERS[m]["T3"] > 0)


def n_eff(m: str) -> int:
    return 27 if MEMBERS[m]["coloured"] else 9


def delta(m: str, selector) -> F:
    return F(2 + selector(m), n_eff(m))


def main() -> None:
    print("ITEM 87 -- smaller-|Q| / sterility selector from closed EM-pair service")
    print("=" * 84)

    print("\n[1] The unbroken charge monitor produces a binary closed-pair test")
    print("    C_close = (T3^2 + Y^2) - (T3 + Y)^2 = -2 T3 Y")
    print(f"    {'member':6s} {'T3':>6s} {'Y':>6s} {'Q':>6s} {'C_close':>9s} {'+1?':>4s}")
    for m in ("nu", "e", "u", "d"):
        t3, y = MEMBERS[m]["T3"], MEMBERS[m]["Y"]
        c = closed_pair_credit(m)
        check(f"{m}: algebraic identity C_close=-2 T3 Y",
              c == -2 * t3 * y)
        print(f"    {m:6s} {str(t3):>6s} {str(y):>6s} {str(q(m)):>6s} {str(c):>9s} {plus_one_closed_pair(m):4d}")

    selected = {m for m in MEMBERS if plus_one_closed_pair(m)}
    check("closed-pair cancellation selects exactly {nu,d}", selected == {"nu", "d"})

    print("\n[2] It is equivalent to the smaller-|Q| member in each weak doublet")
    for dbl, pair in DOUBLETS.items():
        min_q = min(pair, key=lambda x: abs(q(x)))
        positive = [m for m in pair if plus_one_closed_pair(m)]
        print(f"    {dbl:7s}: smaller |Q| = {min_q}; positive C_close = {positive}")
        check(f"{dbl}: positive closed-pair credit iff smaller |Q|", positive == [min_q])

    print("\n[3] It differs from frozen-I3 exactly where the old route failed")
    for label, selector in (
        ("closed-pair EM cancellation", plus_one_closed_pair),
        ("frozen-I3 / up-type reference", plus_one_frozen_i3),
    ):
        print(f"\n    Selector: {label}")
        ok_all = True
        chosen = []
        for m in ("e", "nu", "u", "d"):
            pred = delta(m, selector)
            ok = pred == MEMBERS[m]["target_delta"]
            ok_all &= ok
            if selector(m):
                chosen.append(m)
            print(f"      {m:2s}: A=2+{selector(m)}={2 + selector(m)}, "
                  f"N_eff={n_eff(m):2d}, delta={str(pred):>4s}, "
                  f"target={str(MEMBERS[m]['target_delta']):>4s} {'OK' if ok else 'MISS'}")
        print(f"      +1 set = {{{', '.join(chosen)}}}; all sectors match = {ok_all}")

    check("closed-pair EM selector reproduces all four Koide twists",
          all(delta(m, plus_one_closed_pair) == MEMBERS[m]["target_delta"] for m in MEMBERS))
    check("frozen-I3 selector swaps the quark phases",
          delta("u", plus_one_frozen_i3) == F(1, 9)
          and delta("d", plus_one_frozen_i3) == F(2, 27))

    print("\n[4] CPT and convention checks")
    print("    Under CPT, T3 and Y both flip sign, so T3*Y and C_close are invariant.")
    for m in MEMBERS:
        t3, y = MEMBERS[m]["T3"], MEMBERS[m]["Y"]
        c_anti = (-t3) * (-t3) + (-y) * (-y) - (-(t3 + y)) * (-(t3 + y))
        check(f"{m}: closed-pair selector is CPT invariant", c_anti == closed_pair_credit(m))

    print(
        """
[5] VERDICT -- selector closed at finite charge-monitor grade
    The +1 is the existence of one closed electroweak cancellation-pair record,
    not a continuous charge-weighted multiplier.  The post-EWSB monitored
    service algebra supplies the binary test:

        C_close = (T3^2 + Y^2) - Q^2 = -2 T3 Y > 0.

    This condition is forced by the already-canonical charge map Q=T3+Y and
    selects exactly {nu,d}.  It is equivalent to "smaller |Q|" inside each weak
    doublet, but its derivation is stronger: it says why the smaller-charge
    member gets a record contact, namely because the T3 and Y records close a
    cancellation pair under the unbroken EM monitor.

    The old frozen-I3 rule is now identified as the wrong billing map for this
    purpose: it follows the up-type reference and therefore selects {nu,u},
    swapping the quark phases.  The physical R1-contact selector follows the
    closed EM-pair service record and selects {nu,d}.

    Therefore the remaining +1 selector in the defect-conditioned R1 contact
    theorem is conditionally closed under two already-used framework premises:
      (i) Q=T3+Y is the unbroken monitored charge record;
      (ii) contact multiplicities count closed record-pairs, not continuous
           charge magnitudes.

    Combining this with the previous gate gives the absolute contact ledger
    (e,nu,u,d)=(2,3,2,3) over N_eff=(9,9,27,27).  The Dynamic-R1 mechanism now
    has a sector-native route to A_s=d_s: base d=2 from R3 non-colour channels,
    plus one EM-closed-pair contact exactly on {nu,d}.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
