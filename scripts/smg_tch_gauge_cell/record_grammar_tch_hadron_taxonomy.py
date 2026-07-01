#!/usr/bin/env python3
r"""Closed-shell TCH hadron taxonomy from centre triality.

Purpose
-------
This is the taxonomy companion to the closed-shell Gauss/triality certificate.

The previous script proves the local selection rule:

    a source on a compact TCH colour shell is a gauge-invariant record
    iff its total SU(3) centre triality is zero.

This script applies that rule to the standard hadron families.  It deliberately
does not compute masses, string tensions, widths, or spectroscopy.  It checks
only the exact kinematic question:

    which colour-source words can exist as records on a closed shell?

Result
------
Allowed:
    glueballs / pure Wilson loops       triality 0
    mesons          q qbar              triality 0
    baryons         q q q               triality 0
    antibaryons     qbar qbar qbar      triality 0
    tetraquarks     q q qbar qbar       triality 0
    pentaquarks     q q q q qbar        triality 0
    dibaryons       q^6                 triality 0
    hybrids         q qbar + gluons     triality 0

Forbidden:
    any net-colour word with triality 1 or 2, for example q, qq, q qbar q.

Boundary
--------
This is an exact colour-neutrality / centre-selection result.  It reproduces
the QCD hadron taxonomy, including exotics, but it does not predict which
allowed exotics are dynamically bound or narrow.  Dynamics belongs to the
Wilson-action/bulk rung.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from record_grammar_tch_colour_singlet_gauss import (
    TRIALITY,
    center_projector_value,
    in_gauss_image,
)
from record_grammar_tch_glued_surface_selector import build_glued_complex


@dataclass(frozen=True)
class HadronWord:
    name: str
    quarks: int
    antiquarks: int
    expected_allowed: bool
    comment: str

    @property
    def triality(self) -> int:
        return (self.quarks - self.antiquarks) % TRIALITY

    @property
    def total_partons(self) -> int:
        return self.quarks + self.antiquarks


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<70s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<70s} value={value}")
    if not value:
        raise AssertionError(name)


def assert_false(name: str, value: bool) -> None:
    print(f"  {name:<70s} value={value}")
    if value:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<70s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def source_vector(complex_, word: HadronWord) -> np.ndarray:
    if word.total_partons > complex_.vertices:
        raise AssertionError("not enough vertices for source word")
    charges = np.zeros(complex_.vertices, dtype=int)
    cursor = 0
    for _ in range(word.quarks):
        charges[cursor] = (charges[cursor] + 1) % TRIALITY
        cursor += 1
    for _ in range(word.antiquarks):
        charges[cursor] = (charges[cursor] + 2) % TRIALITY
        cursor += 1
    return charges


def expected_by_center(word: HadronWord) -> bool:
    return word.triality == 0


def word_label(word: HadronWord) -> str:
    return f"{word.name} (q={word.quarks}, qbar={word.antiquarks})"


def named_words() -> list[HadronWord]:
    return [
        HadronWord("glueball / closed Wilson loop", 0, 0, True, "pure gauge, centre-neutral"),
        HadronWord("meson", 1, 1, True, "q qbar"),
        HadronWord("baryon / proton colour sector", 3, 0, True, "qqq epsilon"),
        HadronWord("antibaryon", 0, 3, True, "qbar qbar qbar epsilon"),
        HadronWord("tetraquark", 2, 2, True, "qq qbar qbar"),
        HadronWord("pentaquark", 4, 1, True, "qqqq qbar"),
        HadronWord("hexaquark / dibaryon", 6, 0, True, "q^6"),
        HadronWord("hybrid meson", 1, 1, True, "q qbar plus centre-neutral gluonic loop"),
        HadronWord("lone fundamental", 1, 0, False, "net triality 1"),
        HadronWord("diquark alone", 2, 0, False, "net triality 2"),
        HadronWord("triquark q qbar q", 2, 1, False, "net triality 1"),
        HadronWord("four-quark nonneutral word", 3, 1, False, "net triality 2"),
        HadronWord("five-quark nonneutral word", 3, 2, False, "net triality 1"),
    ]


def main() -> None:
    print("Closed-shell TCH hadron taxonomy from triality")
    print("=" * 96)
    complex_ = build_glued_complex(2)

    print("\n[1] Named hadron families")
    for word in named_words():
        charges = source_vector(complex_, word)
        total_triality = int(np.sum(charges) % TRIALITY)
        gauss_allowed = in_gauss_image(complex_, charges)
        center_allowed = expected_by_center(word)
        projector = center_projector_value(total_triality)

        assert_equal(f"{word_label(word)} total triality", total_triality, word.triality)
        assert_true(f"{word_label(word)} centre rule matches declared status", center_allowed == word.expected_allowed)
        assert_true(f"{word_label(word)} Gauss rule matches centre rule", gauss_allowed == center_allowed)
        if center_allowed:
            assert_less(f"{word_label(word)} centre projector minus 1", abs(projector - 1.0), 1e-12)
            print(f"    allowed:   {word.comment}")
        else:
            assert_less(f"{word_label(word)} centre projector", abs(projector), 1e-12)
            print(f"    forbidden: {word.comment}")

    print("\n[2] Exhaustive small-word check: allowed iff q - qbar = 0 mod 3")
    checked = 0
    allowed = 0
    forbidden = 0
    for quarks in range(0, 7):
        for antiquarks in range(0, 7):
            word = HadronWord("word", quarks, antiquarks, (quarks - antiquarks) % TRIALITY == 0, "")
            charges = source_vector(complex_, word)
            gauss_allowed = in_gauss_image(complex_, charges)
            center_allowed = expected_by_center(word)
            assert_true(f"q={quarks}, qbar={antiquarks} Gauss=center", gauss_allowed == center_allowed)
            checked += 1
            if center_allowed:
                allowed += 1
            else:
                forbidden += 1

    assert_equal("small words checked", checked, 49)
    assert_equal("small words allowed", allowed, 17)
    assert_equal("small words forbidden", forbidden, 32)

    print(
        """
VERDICT:
  PASS.  The closed TCH colour shell reproduces the QCD hadron taxonomy as a
  centre-triality selection rule.  Mesons, baryons, antibaryons, tetraquarks,
  pentaquarks, dibaryons, glueballs, and hybrids are all allowed exactly when
  their total triality is zero.  Lone quarks and all net-colour words are
  forbidden exactly.

  This is kinematic confinement only: it says which colour-source records can
  exist.  It does not decide which allowed exotic channels are dynamically
  bound, narrow, or light.  That requires the bulk Wilson-action dynamics.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
