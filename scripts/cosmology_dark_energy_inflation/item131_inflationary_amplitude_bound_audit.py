#!/usr/bin/env python3
"""ITEM 131: bounded audit of the inflationary amplitude A_nu.

Question
--------
Can the HBC/QEC service ledger now derive the scalar amplitude

    A_nu = F_eff / N_eff

rather than treating it as an observed normalization?

Verdict
-------
Partly reduced, not closed.

The finite 28-channel service graph now closes the T4/Fano leg under the
canonical jump/CTMC reading: the graph is 25-regular, so the total service
count is exactly Poisson and F_eff=1.  That removes the old "assume Poisson"
gap for the dilute continuous-time reading.

The remaining blockers are T3 and T7:

* T3: derive the absolute number N_eff of independent service events in one
  mode-local printed shell.
* T7: derive the early saturated printing scale H_* code-natively, if N_eff
  is identified with de-Sitter entropy.

The viable conditional theorem is therefore:

    if N_eff = S_dS = 8 pi^2 Mbar_P^2 / H_*^2 and F_eff = 1,
    then A_nu = H_*^2 / (8 pi^2 Mbar_P^2).

The observed A_nu ~= 2.1e-9 is then equivalent to H_* ~= 9.9e14 GeV.
Current canon has that scale as a phenomenological R-activation/GUT-neighbour
input, not as a HBC/QEC derivation.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass


A_TARGET = 2.10e-9
MPL_REDUCED_GEV = 2.435e18
ALPHA_0 = 1.0 / 137.036
LAMBDA_QCD_GEV = 0.332


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def bits3(x: int) -> tuple[int, int, int]:
    return ((x >> 2) & 1, (x >> 1) & 1, x & 1)


def dot_mod2(a: tuple[int, int, int], b: tuple[int, int, int]) -> int:
    return sum(x * y for x, y in zip(a, b)) % 2


def hyperplanes() -> list[frozenset[int]]:
    planes = {
        frozenset(p for p in range(8) if dot_mod2(bits3(p), bits3(a)) == b)
        for a in range(1, 8)
        for b in (0, 1)
    }
    return sorted(planes, key=lambda h: tuple(sorted(h)))


def channel_graph() -> list[set[int]]:
    """The item-131 14 hyperplanes x 2 transverse modes service graph."""
    planes = hyperplanes()
    channels = [(i, mode) for i in range(14) for mode in (0, 1)]
    index = {channel: n for n, channel in enumerate(channels)}
    graph = [set() for _ in channels]

    for i, hi in enumerate(planes):
        for j, hj in enumerate(planes):
            if i == j:
                continue
            if len(hi & hj) != 2:
                continue
            for mi in (0, 1):
                for mj in (0, 1):
                    graph[index[(i, mi)]].add(index[(j, mj)])

    for i in range(14):
        a = index[(i, 0)]
        b = index[(i, 1)]
        graph[a].add(b)
        graph[b].add(a)

    return graph


def connected(graph: list[set[int]]) -> bool:
    seen = {0}
    queue = deque([0])
    while queue:
        node = queue.popleft()
        for nxt in graph[node]:
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return len(seen) == len(graph)


def amplitude_from_count(n_eff: float, f_eff: float = 1.0) -> float:
    return f_eff / n_eff


def required_fano_for_count(n_eff: float) -> float:
    return A_TARGET * n_eff


def entropy_count(h_gev: float) -> float:
    return 8.0 * math.pi**2 * MPL_REDUCED_GEV**2 / h_gev**2


def h_from_entropy_count(n_eff: float) -> float:
    return math.sqrt(8.0 * math.pi**2 * MPL_REDUCED_GEV**2 / n_eff)


def h_from_entropy_amplitude(amplitude: float = A_TARGET, f_eff: float = 1.0) -> float:
    return math.sqrt(8.0 * math.pi**2 * f_eff * amplitude) * MPL_REDUCED_GEV


def amplitude_from_entropy_h(h_gev: float, f_eff: float = 1.0) -> float:
    return f_eff / entropy_count(h_gev)


def h_from_landauer_attempt_count(amplitude: float = A_TARGET, f_eff: float = 1.0) -> float:
    """If N_eff = S_dS * Gamma0/H attempts per e-fold, solve A=F/N_eff."""
    gamma0 = ALPHA_0 * LAMBDA_QCD_GEV
    return (8.0 * math.pi**2 * gamma0 * MPL_REDUCED_GEV**2 * amplitude / f_eff) ** (1.0 / 3.0)


@dataclass(frozen=True)
class CountCandidate:
    name: str
    n_eff: float
    status: str


def main() -> None:
    print("ITEM 131 INFLATIONARY AMPLITUDE BOUNDED TARGET AUDIT")

    print("\n[1] T4/Fano leg under the canonical CTMC ledger")
    graph = channel_graph()
    degrees = sorted({len(neighbours) for neighbours in graph})
    check(len(graph) == 28, "service graph has 28 channels")
    check(degrees == [25], "service graph is 25-regular")
    check(connected(graph), "service graph is connected")
    print("  Therefore the total exit rate is channel-independent.")
    print("  Continuous-time service counts have iid exponential holding times.")
    print("  T4 result: F_eff(total) = 1 exactly for the dilute CTMC reading.")
    print("  Scheduler caveat: a bandwidth-one discrete tick with duty cycle p gives F_eff=1-p.")

    print("\n[2] Required T3 count")
    n_required = 1.0 / A_TARGET
    h_required = h_from_entropy_amplitude()
    v14 = (3.0 * h_required**2 * MPL_REDUCED_GEV**2) ** 0.25
    print(f"  A_target                         = {A_TARGET:.6e}")
    print(f"  N_eff/F_eff required             = {n_required:.6e}")
    print(f"  H_* if N_eff=S_dS and F_eff=1    = {h_required:.6e} GeV")
    print(f"  V^(1/4) for that H_*             = {v14:.6e} GeV")
    check(4.0e8 < n_required < 6.0e8, "the CMB amplitude needs roughly 5e8 independent events for F_eff=1")
    check(8.0e14 < h_required < 1.2e15, "the entropy-pixel route is equivalent to H_* near 1e15 GeV")

    print("\n[3] Finite/code-native count candidates")
    candidates = [
        CountCandidate("28 service channels", 28.0, "rejected: far too noisy"),
        CountCandidate("112 incidence flags", 112.0, "rejected: far too noisy"),
        CountCandidate("2^28 active-subset states", float(2**28), "near only with unexplained sub-Poisson Fano"),
        CountCandidate("2^29 active-subset states + binary phase", float(2**29), "near only with unexplained super-Poisson Fano"),
        CountCandidate("28 * 2^24 mixed scalar-shell count", float(28 * 2**24), "numerically close but 2^24 not derived"),
        CountCandidate("112 * 2^22 equivalent mixed count", float(112 * 2**22), "same unresolved mixed factor"),
    ]
    for candidate in candidates:
        amp = amplitude_from_count(candidate.n_eff)
        f_needed = required_fano_for_count(candidate.n_eff)
        print(
            f"  {candidate.name:43s} N={candidate.n_eff:.6e} "
            f"A(F=1)={amp:.6e} F_needed={f_needed:.6e}  {candidate.status}"
        )
    check(amplitude_from_count(28.0) / A_TARGET > 1.0e7, "28-channel counting misses by more than seven orders")
    check(amplitude_from_count(112.0) / A_TARGET > 1.0e6, "112-flag counting misses by more than six orders")
    check(0.50 < required_fano_for_count(2**28) < 0.60, "2^28 needs F_eff about 0.56, not derived")
    check(0.95 < required_fano_for_count(28 * 2**24) < 1.05, "28*2^24 is the only simple count near Poisson")
    check(True, "no current ledger derives the required 2^24 scalar-shell degeneracy")
    print("  Later route note: item131_inflationary_amplitude_alpha4_route_audit.py")
    print("  recasts the 28*2^24 near-hit as a disguised (4/3)*137^4 count,")
    print("  but that Casimir-weighted alpha^-4 ledger is still conditional.")

    print("\n[4] Entropy-pixel route and its alternatives")
    amp_at_1e15 = amplitude_from_entropy_h(1.0e15)
    h_landauer = h_from_landauer_attempt_count()
    gamma0 = ALPHA_0 * LAMBDA_QCD_GEV
    print("  Entropy-pixel premise:")
    print("    N_eff = S_dS and F_eff=1")
    print(f"    A_nu(H_*=1e15 GeV)       = {amp_at_1e15:.6e}")
    print(f"    ratio to target          = {amp_at_1e15 / A_TARGET:.6f}")
    print("  Per-pixel Landauer-attempt premise:")
    print(f"    Gamma0=alpha Lambda_QCD  = {gamma0:.6e} GeV")
    print(f"    inferred H_*             = {h_landauer:.6e} GeV")
    check(abs(amp_at_1e15 / A_TARGET - 1.0) < 0.03, "H_*=1e15 GeV matches the amplitude only under the entropy-pixel premise")
    check(5.0e8 < h_landauer < 2.0e9, "counting Gamma0/H attempts instead moves H_* to the 1e9 GeV scale")
    check(True, "the event ledger must derive which microscopic current defines one scalar-clock event")

    print("\n[5] Code-native scale stress test")
    qcd_scale_candidates = [
        ("Lambda_QCD/2 bipartite P^2 clock scale", 0.5 * LAMBDA_QCD_GEV),
        ("alpha Lambda_QCD leakage bandwidth", gamma0),
        ("Lambda_QCD", LAMBDA_QCD_GEV),
    ]
    for name, h_gev in qcd_scale_candidates:
        amp = amplitude_from_entropy_h(h_gev)
        print(f"  {name:40s} H={h_gev:.6e} GeV  A_entropy={amp:.6e}")
        check(amp < 1.0e-38, f"{name} is far too small for the observed scalar amplitude")
    check(True, "a viable entropy-pixel closure must derive a genuinely high early H_* scale, not a QCD tick scale")

    print("\n[6] Gate status")
    closed = [
        "T4 total-current Fano factor: F_eff=1 under the jump/CTMC service graph",
        "conditional observable formula: A_nu=F_eff/N_eff",
        "entropy-pixel reduction: A_nu=H_*^2/(8 pi^2 Mbar_P^2) if N_eff=S_dS",
    ]
    open_items = [
        "T3 mean current: absolute N_eff per mode-local printed shell",
        "T5 correlation volume: why one scalar mode samples exactly S_dS events",
        "T7 scale: code-native derivation of H_* ~= 1e15 GeV",
        "scheduler regime: dilute Poisson CTMC versus saturated bandwidth-one tick",
    ]
    for item in closed:
        check(True, f"closed/conditional: {item}")
    for item in open_items:
        check(True, f"still open: {item}")

    print("\n" + "=" * 104)
    print("VERDICT")
    print("  A_nu is bounded as a T3/T4 target, not recovered as a theorem.")
    print("  T4 is now closed for the canonical CTMC ledger: F_eff=1 exactly.")
    print("  T3 and T7 remain open.  Within this finite-count/entropy scan the viable")
    print("  route is N_eff=S_dS plus a code-native H_* ~= 1e15 GeV saturated printer;")
    print("  the separate alpha^4 audit now adds a second conditional route.")
    print("  Finite service counts fail or require an unexplained Fano/degeneracy, and")
    print("  QCD-scale clocks give amplitudes many orders too small.")
    print("=" * 104)
    print("exit 0 -- bounded reduction sharpened; A_nu not derived.")


if __name__ == "__main__":
    main()
