# Measurement, Lattice Consistency, and the Dissolution of Retrocausality

## How the Holographic Lattice Resolves the Measurement Problem and Apparent Time-Reversal

**D.G. Elliman — Neuro-Symbolic Ltd, February 2026**

Working note for *It from Bit, Revisited*

---

## 1. The Two Problems

Quantum mechanics has two deep interpretive puzzles that resist conventional resolution:

**The Measurement Problem.** A quantum system evolves smoothly (unitarily) until it is measured, at which point it appears to "collapse" discontinuously into a definite state. No one agrees on what constitutes a measurement, where collapse occurs, or whether collapse is physical or merely epistemic.

**Apparent Retrocausality.** In experiments such as Wheeler's delayed-choice, the entangled-state quantum eraser, and Aharonov-Bergmann-Lebowitz (ABL) weak-value measurements, the choice of what to measure *now* appears to determine what the system was doing *in the past*. The effect seems to precede its cause.

We argue that in the circlette framework, both problems dissolve simultaneously. They are artefacts of assuming that (a) the 3D spatial description is fundamental, and (b) the past state of a system is fully determined before measurement. On the holographic lattice, neither assumption holds.

---

## 2. The Lattice Picture of Measurement

### 2.1 One dynamics, not two

In the circlette framework, there is only one kind of dynamics: local lattice updates governed by the four constraints (R1–R4) and the CNOT rule at the bridge-isospin boundary. There is no separate "measurement postulate." What we call measurement is a specific physical situation — the coupling of a small number of circlettes (the system) to a macroscopic number of circlettes (the apparatus) — handled by the same update rules that govern everything else.

### 2.2 What measurement does at the bit level

A measurement apparatus is engineered to couple to specific bits of the target circlette:

| Measurement type | Bits probed | Sector |
|---|---|---|
| Spin (Stern-Gerlach) | χ (chirality) | Electroweak |
| Charge | Q = f(LQ, I₃) | Bridge + Electroweak |
| Colour (jet fragmentation) | C₀, C₁ | Colour |
| Flavour (oscillation) | G₀, G₁ | Generation |

The apparatus does not and cannot probe all 8 bits simultaneously. Each measurement type reads a *subset* of the ring, through a sector-specific coupling.

### 2.3 Constraint propagation and "collapse"

When the apparatus couples to specific bits, it imposes new boundary conditions on the lattice. The four local rules then propagate these constraints across the ring:

- **R3** (C₁ = 1 ⟹ LQ = 1): A colour measurement constrains the bridge bit.
- **R4** (LQ = 0 ∧ I₃ = 0 ⟹ χ = 0): A bridge constraint propagates to chirality.
- **R2** (χ = W): A chirality constraint locks weak participation.
- **R1** (G₀G₁ ≠ 11): Generation remains independently constrained.

This propagation is deterministic given the measurement outcome — the lattice relaxes to a state consistent with both the measured bits and all four rules. This is "collapse": not a mysterious projection, but local constraint satisfaction propagating across a connected bit-pattern.

### 2.4 Why collapse appears sudden

The apparatus contains ~10²³ circlettes, all coupled to each other through the lattice. When the target circlette's measured bits are pinned, the constraint propagation cascades through the apparatus patterns in ~O(1) lattice ticks (because each constraint spans at most 3 adjacent bits). The macroscopic apparatus cannot maintain superposition because the lattice lacks the computational bandwidth to track coherent amplitudes across ~10²³ entangled rings. Decoherence is not an approximation — it is a *computational resource limit* of the lattice.

---

## 3. The Key Insight: Lattice Consistency Is Global

### 3.1 The error-correcting nature of the lattice

The circlette is an error-correcting code: 45 valid states out of 256 possible 8-bit patterns. The lattice surface on which circlettes propagate must maintain code validity everywhere. When a measurement disturbs specific bits, the lattice must restore global consistency — not just at the spatial location of the measurement, but across the entire causal structure that connects the measured circlette to the rest of the lattice.

This is the crucial point: **the lattice's consistency restoration is not confined to a single time-slice.**

### 3.2 The holographic surface has no preferred time direction

On the 2D holographic boundary, "time" as experienced in the 3D bulk is encoded in the lattice structure — it is not a separate parameter. The lattice has spatial extent and update-step structure, but the constraint rules are *time-symmetric*: R1–R4 are purely structural (they constrain bit patterns, not their direction of evolution), and the CNOT rule is an involution (applying it twice returns to the original state: M² = I).

This means that when the lattice enforces consistency after a measurement, the constraint propagation runs in *all directions on the lattice* — including what, from the 3D perspective, we would call "the past."

### 3.3 What "propagation into the past" actually means

This is **not** retrocausation. The lattice is not sending signals backwards in time. What is happening is more subtle and more profound:

Before the measurement, the lattice configuration in the "past" region was *not fully determined*. It existed as a superposition of configurations — multiple valid bit-patterns consistent with all the information available up to that point. The lattice had not yet "committed" to a specific history because no interaction had forced it to.

When the measurement occurs, it provides new information — a new boundary condition. The lattice must now find a configuration that is consistent with:

1. The measured outcome (new boundary condition).
2. All four local rules (R1–R4) at every lattice site.
3. The CNOT dynamics at every update step.
4. All *other* boundary conditions (preparation, earlier measurements, apparatus settings).

The unique globally consistent configuration may differ from what we *expected* the past state to be, given only the preparation information. But the past state was never what we expected — it was undetermined until the full set of boundary conditions was established.

**The lattice does not change the past. It reveals that the past was never what we assumed.**

---

## 4. Resolution of Specific "Retrocausal" Experiments

### 4.1 Wheeler's Delayed-Choice Experiment

**Setup:** A photon enters an interferometer. After it has "already passed" the first beam splitter, the experimenter chooses whether to insert the second beam splitter. With the second splitter, interference is observed (wave behaviour). Without it, which-path information is available (particle behaviour). The photon appears to "decide" retroactively whether to be a wave or particle.

**Lattice resolution:** The photon circlette enters the interferometer and the lattice maintains superposition across both paths — the bit-pattern has amplitudes on multiple lattice sites. The "past" configuration (which path the circlette took) is *genuinely undetermined* at this point. The lattice has not committed.

When the experimenter inserts (or doesn't insert) the second beam splitter, this changes the boundary conditions at the output. The lattice finds the globally consistent configuration:

- **Second splitter present:** The consistent configuration has the circlette's propagation mode coherently spanning both paths → interference.
- **Second splitter absent:** The consistent configuration has the circlette's propagation mode localised to one path → which-path information.

The photon did not "go back and change its mind." The lattice simply had two different sets of boundary conditions, leading to two different globally consistent histories. The past was never determined until the full boundary was set.

### 4.2 Quantum Eraser

**Setup:** Entangled photon pairs are created. One photon (signal) goes through a double slit. The other (idler) is measured either in a way that preserves or erases which-path information. Erasing the idler's which-path information *restores interference* in the signal — even if the idler is measured after the signal has hit the screen.

**Lattice resolution:** The entangled pair consists of two circlettes whose bits are correlated on the lattice (their sectors share constraint linkages through the lattice surface). The signal circlette hits the screen, and its position is recorded — but the *full* lattice configuration is not yet determined, because the idler's bits have not been pinned.

When the idler is measured:

- **Which-path measurement:** The idler's bits are pinned in a way that constrains the signal's path information. The lattice finds a consistent history in which the signal traversed one slit. No interference.
- **Eraser measurement:** The idler's bits are pinned in a way that does *not* constrain path. The lattice finds a consistent history in which the signal's propagation mode spanned both slits. Interference is restored — but only visible in the coincidence counts (correlating signal and idler outcomes).

The signal's pattern on the screen was there all along — both the interference and non-interference sub-ensembles were present. The idler measurement simply determines which sub-ensemble you can *identify*. The lattice doesn't change the screen pattern retroactively; it determines which lattice history (interfering or non-interfering) is consistent with the full set of boundary conditions.

### 4.3 ABL Weak Values and Pre- and Post-Selection

**Setup:** A system is prepared (pre-selected) in state |ψᵢ⟩ and later post-selected in state |ψ_f⟩. Weak measurements between preparation and post-selection yield "weak values" that can lie outside the eigenvalue spectrum and appear to reveal bizarre intermediate states.

**Lattice resolution:** Pre-selection pins certain bits at the "early" boundary. Post-selection pins (possibly different) bits at the "late" boundary. The lattice must find a configuration consistent with *both* boundary conditions and all intermediate dynamics.

The weak value is the lattice's "best interpolation" between the two boundaries. It is the Fisher-information-weighted average of the bit configurations consistent with both constraints. Anomalous weak values (outside the eigenvalue range) arise when the two boundary conditions are nearly orthogonal — the lattice is being asked to interpolate between nearly incompatible constraints, and the resulting average is dominated by rare, extreme configurations.

This is not retrocausation. It is *bidirectional constraint satisfaction* on a lattice with time-symmetric rules.

---

## 5. The Bohmian Analogy

David Bohm's pilot wave theory resolves many of the same puzzles by introducing a global wave function ψ that "guides" particles along definite trajectories. The wave is non-local: it responds instantaneously to changes in boundary conditions (such as opening or closing a slit), and it effectively provides the holistic context that determines the particle's behaviour.

The circlette lattice plays an analogous role, but with important differences:

| Bohm's pilot wave | Circlette lattice |
|---|---|
| Wave function ψ defined on configuration space | Fisher information field F defined on the 2D lattice |
| Guides point particles along trajectories | Guides circlette bit-patterns through lattice updates |
| Non-local: changes everywhere when boundary changes | Non-local: constraint propagation across the lattice surface |
| Deterministic given initial conditions | Deterministic given *all* boundary conditions (initial + final) |
| ψ is postulated | F emerges from the lattice structure |
| Measurement changes ψ globally | Measurement imposes new boundary conditions; lattice re-satisfies globally |
| Empty branches persist forever | No branches — one lattice, one consistent history |

The key improvement over Bohm is that the lattice's "adaptation" has a clear physical mechanism — error-correcting code consistency enforcement — rather than being a mysterious action-at-a-distance of the wave function. The lattice *is* the wave function, written in its most fundamental language: bits on a surface with local rules.

### 5.1 How the lattice "adapts" like Bohm's wave

When Bohm's wave encounters a new boundary condition (a measurement), the entire wave function updates globally. This looks like instantaneous non-local action, which troubled Einstein.

On the lattice, the same effect occurs but with a clear mechanism:

1. **The measurement pins bits** at a specific lattice location.
2. **Constraint propagation** radiates outward from the measurement site at the lattice update speed (one cell per tick = one Planck length per Planck time = c).
3. **The propagation is causal** on the lattice surface — it respects the light cone.
4. **But** the holographic encoding means that lattice-surface causality maps non-trivially to 3D-bulk causality. What looks instantaneous or retrocausal from the 3D perspective is perfectly causal on the 2D surface.

This is the holographic resolution of non-locality: the lattice is the *more fundamental* description, and its causal structure is local and well-defined. The apparent non-locality of quantum mechanics (Bell violations, retrocausality, measurement collapse) arises because we are trying to describe a 2D-local process in a 3D projection that does not preserve all the lattice's causal relationships.

---

## 6. The Two-Boundary Interpretation

The circlette framework naturally leads to what might be called a **two-boundary interpretation** of quantum mechanics:

1. **Preparation** (initial boundary): pins certain bits of the circlette at the start of the experiment. This is the initial state |ψᵢ⟩.

2. **Measurement** (final boundary): pins certain bits at the end. This is the final state |ψ_f⟩ or the measurement outcome.

3. **The lattice interpolates** between the two boundaries, finding the minimum-action configuration consistent with both. This interpolation *is* the quantum evolution.

4. **Probabilities** arise because, before the final boundary is set, multiple interpolations are possible. The Born rule should emerge from the Fisher information measure on the space of consistent interpolations.

This is structurally similar to the two-state vector formalism of Aharonov, Bergmann, and Lebowitz (1964), and to the transactional interpretation of Cramer (1986). But it differs from both in having a concrete physical substrate (the lattice) and a specific mechanism (constraint propagation in an error-correcting code).

### 6.1 Why this isn't fatalism

A natural objection: if the future measurement partly determines the past state, doesn't this imply predetermination?

No. The future boundary is not "already there" waiting. It is established when and only when the measurement interaction occurs. Before the measurement, the future boundary is open — the lattice has not committed. The consistency enforcement is a constraint-satisfaction process that happens *when* the measurement coupling provides new information, not a pre-existing block universe that we merely discover.

The analogy is to a crossword puzzle: filling in a late clue constrains what the early answers must have been, but the puzzle was genuinely unsolved before you filled in the clue. The lattice is an ongoing crossword, with each measurement adding a new clue that retroactively constrains — but does not change — the partially filled grid.

---

## 7. Experimental Signatures

The lattice-consistency picture makes predictions that differ subtly from standard quantum mechanics:

1. **Constraint propagation has finite speed.** On the lattice, consistency enforcement propagates at c (one cell per tick). For experiments where the measurement's causal influence has not had time to reach the "past" region in question, the lattice cannot enforce full consistency. This would produce subtle deviations from standard QM predictions in *extreme* delayed-choice experiments with spacelike separation between the measurement event and the region whose "past" behaviour is being inferred. Current experiments are far from this regime, but it is in principle testable.

2. **Weak measurements as partial constraint propagation.** In the circlette picture, a weak measurement couples to only a few bits with low amplitude, causing *partial* constraint propagation that doesn't cascade to the full ring. This predicts a smooth transition between "no measurement" and "strong measurement" as coupling strength increases — consistent with existing weak measurement data, but with a specific bit-level model that could be tested against the detailed statistics.

3. **Sector-dependent decoherence rates.** Because different measurements probe different ring sectors, and the constraints linking sectors have different structures (some span 2 bits, some 3), the rate at which superpositions decohere should depend on *which property* is being measured. Electroweak sector decoherence (involving the CNOT dynamics) should differ from colour sector decoherence (involving only static constraints). This is testable in principle with high-precision interferometry.

4. **No retrocausal signalling.** The lattice consistency enforcement cannot transmit information backwards in time, because it only selects among *already-existing* amplitudes — it does not create new ones. This is consistent with all existing no-signalling theorems and distinguishes the circlette picture from genuinely retrocausal theories.

---

## 8. Summary

The measurement problem and apparent retrocausality are both artefacts of an incomplete description:

| Puzzle | Standard QM | Circlette resolution |
|---|---|---|
| What is measurement? | Undefined | Macroscopic coupling that pins specific bits |
| Why does collapse occur? | Postulated | Constraint propagation + computational bandwidth limit |
| Why is collapse sudden? | Mysterious | O(1)-tick propagation across locally-constrained ring |
| Why does the past seem to change? | Retrocausation? | Past was never determined; lattice finds globally consistent history |
| Why can't we signal backwards? | No-signalling theorem | Consistency enforcement selects among existing amplitudes |
| What is the wave function? | Ontologically ambiguous | Fisher information field on the lattice surface |
| Why does ψ adapt globally? | Non-local action | Lattice-local constraint propagation, holographically projected |

The circlette lattice provides a single, unified mechanism — error-correcting code consistency enforcement on a holographic surface — that accounts for measurement, collapse, decoherence, non-locality, and apparent retrocausality. No separate measurement postulate is needed. No many-worlds branching is required. No pilot wave is assumed. The lattice does what lattices do: maintain consistency, one cell at a time.

The connection to Bohm is instructive but not exact. Bohm's pilot wave *is* the lattice's Fisher information field, projected into 3D. But where Bohm postulated the wave and its guidance equation, the circlette framework *derives* both from the more primitive structure of bits on a surface with local rules. The lattice adapts as Bohm's wave adapts — holistically, responsively, maintaining global consistency — but it does so through a mechanism that is fully local on the fundamental (2D) description.

The apparent weirdness of quantum mechanics is the weirdness of mistaking a 2D-local process for a 3D-local one. On the lattice, everything is causal, consistent, and comprehensible. The price is accepting that what we call "spacetime" is not fundamental but emergent — a holographic projection of information patterns on a surface. Wheeler's intuition was correct: it really is "it from bit," all the way down.

---

*For the main paper, see Elliman (2026), "It from Bit, Revisited." For the Cheshire Cat analysis, see the companion working note on property stripping and the bridge bit.*
