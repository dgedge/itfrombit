#!/usr/bin/env python3
r"""R4 equation-of-state resolution and CMB completion gate.

This script supersedes the active-EoS reading in item123_cmb_completion.py.

Question
--------
Canon carried an apparent contradiction:

* R4 emitted as Landauer exhaust has w = +1/3.
* An older Item-132/string-gas reading said w = -1/3.

The current Item-132/MOND chain no longer uses the string-gas/negative-pressure
route.  It uses scheduler-clocked Poisson line-current records.  Therefore the
old -1/3 statement is not a second active cosmological R4 phase; it is a retired
branch.

The CMB gate is separate: the third acoustic peak needs a conserved,
collisionless, pressureless, pre-recombination clustering component with
Omega_c h^2 ~= 0.120.  The framework's 20% sterile-nu_R component can count.
The 80% R4 sector counts only if it can be derived as early w=0 dust, or if a
relativistic/AeST-like completion makes it mimic dust in the Boltzmann system.

Result
------
The EoS contradiction resolves, but the CMB completion does not:

* active microscopic R4 exhaust: w=+1/3, radiation-like, not CDM;
* active bound/halo R4: w_eff > 0 line-current/topological fluid, nonlinear and
  pressureful, not homogeneous recombination CDM;
* retired string-gas R4: w=-1/3, inactive;
* hypothetical early R4 dust: would close the budget exactly, but is not
  derivable from current scheduler/line-current mechanics.

exit 0 means the bookkeeping assertions pass and the no-go is the current
canon-level verdict.
"""

from __future__ import annotations

from dataclasses import dataclass


OMEGA_C_H2_NEEDED = 0.120
OMEGA_B_H2 = 0.0224
OMEGA_DM_H2 = 0.120
F_NU_R = 0.20
OMEGA_NUR_H2 = F_NU_R * OMEGA_DM_H2
OMEGA_R4_H2 = (1.0 - F_NU_R) * OMEGA_DM_H2


@dataclass(frozen=True)
class R4Regime:
    name: str
    active: bool
    w_label: str
    counts_as_cmb_cdm: bool
    reason: str


REGIMES = (
    R4Regime(
        name="microscopic unbound Landauer exhaust",
        active=True,
        w_label="+1/3",
        counts_as_cmb_cdm=False,
        reason="radiation-like emission/free-streaming pressure; redshifts as a^-4",
    ),
    R4Regime(
        name="bound halo line-current/topological fluid",
        active=True,
        w_label="w_eff > 0",
        counts_as_cmb_cdm=False,
        reason=(
            "scheduler-clocked Poisson line records with pressure/stiffness; "
            "a nonlinear halo response, not a homogeneous conserved dust fluid"
        ),
    ),
    R4Regime(
        name="retired constant-tension/string-gas branch",
        active=False,
        w_label="-1/3",
        counts_as_cmb_cdm=False,
        reason=(
            "retired with the negative-pressure/constant-tension Jeans route; "
            "not an active Item-132 mechanism"
        ),
    ),
    R4Regime(
        name="hypothetical early conserved R4 dust",
        active=False,
        w_label="0",
        counts_as_cmb_cdm=True,
        reason="would close the CMB budget, but no current R4/Kraus/scheduler theorem derives it",
    ),
)


DUST_CRITERIA = {
    "conserved comoving number": False,
    "pressureless/small sound speed": False,
    "independent pre-recombination clustering source": False,
    "homogeneous Boltzmann component before nonlinear halos": False,
}


def main() -> None:
    print("R4 EOS / CMB COMPLETION RESOLUTION")

    print("\n[1] Active-vs-retired R4 equation-of-state ledger")
    for regime in REGIMES:
        state = "ACTIVE" if regime.active else "inactive/retired"
        cmb = "counts" if regime.counts_as_cmb_cdm and regime.active else "does not count"
        if regime.name.startswith("hypothetical"):
            cmb = "would count if derived"
        print(f"  {state:16s}  w={regime.w_label:9s}  {regime.name}")
        print(f"    -> CMB: {cmb}; {regime.reason}")

    active_ws = {r.w_label for r in REGIMES if r.active}
    assert "+1/3" in active_ws
    assert "-1/3" not in active_ws
    assert "0" not in active_ws
    assert not any(r.counts_as_cmb_cdm for r in REGIMES if r.active)
    print("  [PASS] no active R4 branch is w=-1/3 or w=0 dust")

    print("\n[2] CMB third-peak anchor budget")
    anchor_current = OMEGA_NUR_H2
    shortfall = OMEGA_C_H2_NEEDED / anchor_current
    total_clustering_current = OMEGA_B_H2 + anchor_current
    total_clustering_lcdm = OMEGA_B_H2 + OMEGA_C_H2_NEEDED
    print(f"  required cold clustering component: Omega_c h^2 = {OMEGA_C_H2_NEEDED:.3f}")
    print(f"  nu_R anchor (20% of dark sector):    Omega h^2 = {OMEGA_NUR_H2:.3f}")
    print(f"  R4 non-dust sector (80%):            Omega h^2 = {OMEGA_R4_H2:.3f}  [not counted]")
    print(f"  current anchor shortfall:            {shortfall:.1f}x")
    print(
        f"  total clustering matter: framework {total_clustering_current:.4f} "
        f"vs LCDM {total_clustering_lcdm:.4f} ({total_clustering_current / total_clustering_lcdm:.2f}x)"
    )
    assert abs(OMEGA_NUR_H2 - 0.024) < 1e-12
    assert abs(OMEGA_R4_H2 - 0.096) < 1e-12
    assert 4.8 < shortfall < 5.2

    print("\n[3] Exact identity if R4 were early dust")
    anchor_if_r4_dust = OMEGA_NUR_H2 + OMEGA_R4_H2
    print(f"  nu_R + R4-as-dust = {anchor_if_r4_dust:.3f} h^2")
    print(f"  target            = {OMEGA_C_H2_NEEDED:.3f} h^2")
    assert abs(anchor_if_r4_dust - OMEGA_C_H2_NEEDED) < 1e-12
    print("  [PASS] the budget would close exactly if, and only if, R4 dust were derived")

    print("\n[4] Dust-derivability test under current R4 mechanics")
    for criterion, ok in DUST_CRITERIA.items():
        print(f"  [{'YES' if ok else 'NO '}] {criterion}")
    assert not any(DUST_CRITERIA.values())
    print(
        "  Current mechanics give non-conserved service records and late line-current/halo response. "
        "They do not give a conserved collisionless species in the recombination Boltzmann system."
    )

    print("\n[5] Verdict")
    print(
        "  The apparent R4 EoS contradiction is resolved by status separation: +1/3 is the active "
        "microscopic emission reading; w_eff>0 is the active bound/halo line-current reading; "
        "-1/3 is a retired string-gas/negative-pressure branch."
    )
    print(
        "  The CMB completion remains a budget-level no-go in current canon.  Only nu_R supplies "
        "cold pre-recombination clustering, giving a 5x shortfall.  The exact escape target is "
        "now narrow: derive an early conserved R4 dust/AeST-like mode, then run the Boltzmann "
        "system.  Without that new theorem, R4 cannot complete the CMB."
    )
    print("exit 0 -- R4 EoS contradiction resolved; CMB completion still open/no-go.")


if __name__ == "__main__":
    main()
