#!/usr/bin/env python3
"""
Run the complete fermion doubling analysis.

Executes all modules in sequence and generates all figures.

Usage:
    python run_all.py
"""

import sys
import os

# Ensure we're in the right directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

SEPARATOR = "=" * 70


def run_section(title, func):
    print(f"\n{SEPARATOR}")
    print(f"  {title}")
    print(f"{SEPARATOR}\n")
    func()


def section_1():
    """Dirac algebra verification."""
    from dirac_matrices import verify_clifford
    verify_clifford()


def section_2():
    """Naive square lattice — fermion doubling baseline."""
    import naive_square
    # Re-run __main__ logic
    import importlib
    importlib.reload(naive_square)
    exec(open("naive_square.py").read())


def section_3():
    """Wilson fermions — standard fix."""
    exec(open("wilson_fermions.py").read())


def section_4():
    """4.8.8 lattice — physical NNN coupling analysis."""
    exec(open("lattice_488.py").read())


def section_5():
    """Systematic BZ scan."""
    exec(open("dispersion_scan.py").read())


def section_6():
    """Band structure and dispersion plots."""
    from band_structure import plot_band_comparison, plot_dispersion_surfaces
    plot_band_comparison("bands_comparison.png")
    plot_dispersion_surfaces("dispersion_surface.png")


def main():
    print(SEPARATOR)
    print("  FERMION DOUBLING ON THE 4.8.8 LATTICE")
    print("  Companion analysis for 'The Holographic Circlette'")
    print("  D.J. Shorten, Neuro-Symbolic Ltd, 2026")
    print(SEPARATOR)

    run_section("1. DIRAC ALGEBRA VERIFICATION", section_1)
    run_section("2. NAIVE SQUARE LATTICE (BASELINE)", section_2)
    run_section("3. WILSON FERMIONS (STANDARD FIX)", section_3)
    run_section("4. 4.8.8 LATTICE (PHYSICAL NNN COUPLING)", section_4)
    run_section("5. SYSTEMATIC BRILLOUIN ZONE SCAN", section_5)
    run_section("6. GENERATING PLOTS", section_6)

    print(f"\n{SEPARATOR}")
    print("  ANALYSIS COMPLETE")
    print(SEPARATOR)
    print()
    print("Key findings:")
    print("  1. Naive square lattice: 4 Dirac points (1 physical + 3 doublers)")
    print("  2. Wilson mass (by hand): all doublers gapped ✓")
    print("  3. 4.8.8 physical NNN (α₁α₂ × sin·sin): doublers PERSIST ✗")
    print("  4. Resolution: CNOT coin = dynamical Z₂ chiral mixing")
    print("     → Nielsen-Ninomiya inapplicable to discrete chirality")
    print()
    print("Generated files:")
    print("  bands_comparison.png   — Band structures along Γ→X→M→Γ")
    print("  dispersion_surface.png — |E|_min across full BZ")
    print()


if __name__ == "__main__":
    main()
