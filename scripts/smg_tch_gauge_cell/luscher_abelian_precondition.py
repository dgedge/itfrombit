#!/usr/bin/env python3
"""Luscher abelian-sector SOLVABILITY PRECONDITION for the framework's 16-fermion content.

Self-asserting. Computes the full Standard-Model gauge + gravitational + Witten anomaly audit
over one generation (the 16 = SO(10) chiral spinor, item 116 / §2.11, incl. nu_R), in exact
rational arithmetic. All channels vanish => Luscher's 1999 abelian U(1)_Y chiral-gauge measure
is GLOBALLY constructible. This is the *precondition* invoked by the DRIFT K4 Luscher-measure
route (DRIFT line 909).

SCOPE (read before citing): this is the abelian solvability PRECONDITION only, NOT the chiral
measure construction. The non-abelian SU(2)_L measure phase is the open piece (DRIFT K4
§909/§911/§915 — closed-negative at the continuum, with zero geometric leverage by the homotopy
theorem). So a green exit here is a re-derivation / extension of the §2.11 anomaly content
(`smg_code_projection.py`: sum Q = 0, sum Q^2 = 16; §911/§919: sum Y = 0, sum Y^3 = 0) plus the
mixed / SU(3)^3 / Witten channels made explicit. It does NOT advance items 13/143.

Convention: all fields as left-handed Weyl (right-handed fields conjugated); Q = T3 + Y with
Y(Q_L) = 1/6.
"""
from fractions import Fraction as F

# (name, Y, SU(2) dim, signed SU(3): +3 triplet / -3 anti-triplet / +1 singlet)
FERMIONS = [
    ("Q_L",    F(1, 6),  2, +3),   # quark doublet
    ("u_R^c",  F(-2, 3), 1, -3),   # up antiquark
    ("d_R^c",  F(1, 3),  1, -3),   # down antiquark
    ("L_L",    F(-1, 2), 2, +1),   # lepton doublet
    ("e_R^c",  F(1),     1, +1),   # positron
    ("nu_R^c", F(0),     1, +1),   # right-handed neutrino (the framework's 16th)
]


def mult(f):
    """Total Weyl components = SU(2) dim x |SU(3) dim|."""
    return f[2] * abs(f[3])


def anomalies(fermions):
    return {
        "n_weyl":  sum(mult(f) for f in fermions),
        "grav":    sum(mult(f) * f[1] for f in fermions),                        # U(1)_Y - grav^2
        "cubic":   sum(mult(f) * f[1] ** 3 for f in fermions),                   # U(1)_Y^3
        "mixed2":  sum(abs(f[3]) * f[1] for f in fermions if f[2] == 2),         # U(1)_Y - SU(2)^2
        "mixed3":  sum(f[2] * f[1] for f in fermions if abs(f[3]) == 3),         # U(1)_Y - SU(3)^2
        "su3cube": sum(f[2] * (1 if f[3] > 0 else -1)                            # SU(3)^3 (vector-like?)
                       for f in fermions if abs(f[3]) == 3),
        "n_doub":  sum(abs(f[3]) for f in fermions if f[2] == 2),                # Witten SU(2): even?
    }


def main():
    a = anomalies(FERMIONS)
    print(f"Weyl components: {a['n_weyl']}")
    for f in FERMIONS:
        print(f"  {f[0]:6s} Y={str(f[1]):>5s}  SU2={f[2]} SU3={f[3]:+d}  (x{mult(f)})")
    print(f"  sum Y        (grav)       = {a['grav']}")
    print(f"  sum Y^3      (U(1)^3)     = {a['cubic']}")
    print(f"  sum Y|SU(2)  (U1-SU2^2)   = {a['mixed2']}")
    print(f"  sum Y|SU(3)  (U1-SU3^2)   = {a['mixed3']}")
    print(f"  SU(3)^3 signed            = {a['su3cube']}")
    print(f"  # SU(2) doublets (Witten) = {a['n_doub']}")

    assert a["n_weyl"] == 16, a["n_weyl"]
    assert a["grav"] == 0, a["grav"]
    assert a["cubic"] == 0, a["cubic"]
    assert a["mixed2"] == 0, a["mixed2"]
    assert a["mixed3"] == 0, a["mixed3"]
    assert a["su3cube"] == 0, a["su3cube"]
    assert a["n_doub"] % 2 == 0, a["n_doub"]

    print("\nPASS: U(1)^3, U(1)-grav^2, U(1)-SU(2)^2, U(1)-SU(3)^2, SU(3)^3, Witten all cancel over the 16.")
    print("=> Luscher-1999 abelian U(1)_Y chiral measure is GLOBALLY constructible (precondition MET).")
    print("SCOPE: precondition only; the non-abelian SU(2)_L measure phase (the open piece) is untouched.")


if __name__ == "__main__":
    main()
