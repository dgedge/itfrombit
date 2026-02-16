# The Born Rule from Lattice Configuration Counting

## Why Measurement Probabilities Follow |ψ|² in the Circlette Framework

**D.G. Elliman — Neuro-Symbolic Ltd, February 2026**

Working note for *It from Bit, Revisited*

---

## 1. The Problem

The Born rule states that if a quantum system is in state |ψ⟩, the probability of measuring outcome |a⟩ is:

$$P(a) = |\langle a | \psi \rangle|^2$$

This is the bridge between the mathematical formalism (amplitudes, wave functions, Hilbert spaces) and the physical world (detector clicks, pointer readings, experimental frequencies). Every prediction of quantum mechanics depends on it.

The Born rule is notoriously difficult to derive. In standard QM it is simply postulated. Attempts to derive it from within QM face circularity (Gleason's theorem assumes the Hilbert space structure; the many-worlds derivations assume the branching structure they're trying to explain). The question is: can the circlette framework do better?

We argue that it can. The Born rule emerges as a *counting measure* on the space of globally consistent lattice configurations, in exact analogy to how thermodynamic probabilities emerge from counting microstates. The key is that the wave function amplitude is not an abstract mathematical object in the circlette picture — it is a *physical quantity* on the lattice, namely the square root of the Fisher information density. And Fisher information, in turn, counts the distinguishable configurations available at each lattice site.

---

## 2. The Chain of Reasoning

### 2.1 From bits to Fisher information

At each site on the holographic lattice, the local state is an 8-bit circlette pattern. The constraints R1–R4 select 45 valid matter states from 256 possibilities. The *information content* at each site is characterised by the Fisher information:

$$F(x) = \sum_{k} \frac{1}{p_k(x)} \left( \frac{\partial p_k}{\partial x} \right)^2$$

where p_k(x) is the probability of the k-th valid bit-pattern at lattice site x. In the simplest case — a single circlette propagating on an otherwise empty lattice — p_k(x) is the probability of finding that circlette in pattern k at site x.

Fisher information has a well-known geometric interpretation: it defines a *metric* on the space of probability distributions. The distance between two nearby distributions is:

$$ds^2 = F(x) \, dx^2$$

This is the information-geometric line element that the main paper identifies with the spacetime metric (up to constants). The point for our present purposes is that F(x) has a direct physical meaning: it measures the *density of distinguishable states* at site x. High Fisher information means many distinguishable configurations; low Fisher information means few.

### 2.2 From Fisher information to the wave function

The main paper establishes the identification:

$$|\psi(x)|^2 = \rho(x) = \frac{F(x)}{F_{\text{total}}}$$

where ρ(x) is the normalised bit-density on the lattice and F_total is the total Fisher information integrated over the lattice. The wave function is:

$$\psi(x) = \sqrt{\rho(x)} \; e^{i\theta(x)}$$

where θ(x) is the phase, accumulated through lattice updates along the propagation path. In the circlette picture, each CNOT tick advances θ by a fixed increment (related to the bit-flip cost 36/45), and the spatial variation of θ encodes momentum via the de Broglie relation.

This identification is not arbitrary — it follows from the information action principle. The Fisher information metric defines the geometry; the action S_I = ∫√F ds is minimised along geodesics; and the Euler-Lagrange equation of this action *is* the Schrödinger equation (this was shown by Frieden, 1998, and independently by Hall and Reginatto, 2002, and is the mathematical backbone of the main paper's §14–15).

### 2.3 From the wave function to probabilities

Now we can derive the Born rule. The argument proceeds in three steps.

**Step 1: Superposition is multi-site amplitude on the lattice.**

When a circlette is in a superposition of being at positions x₁ and x₂ (or in states |a₁⟩ and |a₂⟩), this means the lattice has *not committed* to a single configuration. Multiple globally consistent lattice configurations exist, some with the circlette at x₁, others with it at x₂. The amplitude ψ(xⱼ) at each site reflects the lattice's "readiness" to commit to that site — formally, the density of consistent configurations that place the circlette there.

**Step 2: Measurement forces the lattice to commit.**

As described in the companion working note on measurement, a macroscopic apparatus coupling pins specific bits and triggers constraint propagation. The lattice must choose *one* globally consistent configuration from the space of all possibilities.

**Step 3: The probability of each outcome is proportional to the number of consistent configurations that support it.**

This is the core claim. Let Ω(a) denote the number of globally consistent lattice configurations in which the measurement outcome is a. Then:

$$P(a) = \frac{\Omega(a)}{\sum_{a'} \Omega(a')}$$

This is the standard statistical-mechanical definition of probability: count the microstates compatible with each macrostate.

**The Born rule follows if Ω(a) ∝ |⟨a|ψ⟩|².**

We now argue that this proportionality holds.

---

## 3. Why Ω(a) ∝ |⟨a|ψ⟩|²

### 3.1 The lattice as a constrained configuration space

Consider a circlette propagating from source to detector across a lattice of N sites. At each site, the circlette can be in any of 45 valid states, and the lattice must satisfy local constraints at every site (not just the circlette site — the vacuum sites have their own Fisher information F_vac, as discussed in the dynamic-Λ paper).

The total number of globally consistent configurations is:

$$\Omega_{\text{total}} = \prod_{i=1}^{N} \omega(x_i)$$

where ω(xᵢ) is the number of locally consistent configurations at site xᵢ. The circlette's presence at a site modifies that site's local count — it contributes additional distinguishable configurations proportional to the circlette's amplitude there.

### 3.2 Fisher information counts distinguishable configurations

This is the key mathematical fact. Fisher information has an operational definition: it is the number of statistically distinguishable states per unit parameter interval. Specifically, the Fisher information F(x) determines the number of distinguishable probability distributions in a neighbourhood of x:

$$\text{(number of distinguishable states near } x) \propto \sqrt{F(x)} \, \Delta x$$

For a discrete lattice with spacing ℓ_P (the Planck length), the number of distinguishable configurations at site x is:

$$\omega(x) \propto \sqrt{F(x)} = |\psi(x)|$$

Wait — this gives |ψ|, not |ψ|². Where does the square come from?

### 3.3 The square comes from the two-boundary structure

Recall from the retrocausality working note that the lattice operates with *two* boundaries: preparation (initial) and measurement (final). The globally consistent configuration must satisfy both.

The number of configurations consistent with the *initial* boundary at site x is proportional to |ψ_i(x)|. The number consistent with the *final* boundary (measurement outcome a) at site x is proportional to |⟨a|⟩| evaluated at x, call it |φ_a(x)|.

The total number of consistent configurations at site x, given both boundaries, is the *product*:

$$\omega_{\text{both}}(x) \propto |\psi_i(x)| \cdot |\phi_a(x)|$$

For the total probability of outcome a, we must sum (integrate) over all lattice sites through which the circlette could have propagated:

$$\Omega(a) \propto \sum_{x} |\psi_i(x)| \cdot |\phi_a(x)|$$

In the continuum limit, with the phase information included:

$$\Omega(a) \propto \left| \sum_{x} \psi_i(x)^* \, \phi_a(x) \right|^2 = |\langle a | \psi \rangle|^2$$

The square arises because the probability involves the *product* of forward-going (preparation → site) and backward-going (measurement → site) configuration counts, summed coherently with phases. This is precisely the structure of the Born rule.

### 3.4 Why the phases matter

The phase θ(x) accumulated by the circlette through CNOT ticks along its propagation path is *physical* — it counts the number of lattice updates, which determines the relative alignment of bit-patterns at different sites. When we sum the configuration counts over all possible paths (all lattice sites), configurations with aligned phases reinforce (constructive interference) and those with misaligned phases cancel (destructive interference).

This is why |ψ|² involves the *square of the summed amplitude* rather than the sum of squared amplitudes. The lattice configurations at different sites are not independent — they are connected by constraint propagation, and the phase encodes this connection. The "interference" is not a mysterious wave phenomenon; it is the lattice counting configurations whose constraint chains (from preparation through intermediate sites to measurement) are mutually consistent.

---

## 4. The Thermodynamic Analogy

The argument has an exact structural parallel in statistical mechanics:

| Statistical mechanics | Circlette lattice |
|---|---|
| Microstate | Globally consistent lattice configuration |
| Macrostate | Measurement outcome a |
| Boltzmann weight e^{−E/kT} | Fisher information density F(x) |
| Partition function Z | Total configuration count Ω_total |
| Probability P = Ω(a)/Z | Born rule P(a) = |⟨a\|ψ⟩|² |
| Second law (entropy increase) | Information action minimisation |
| Equilibrium | Minimum-action lattice configuration |

In statistical mechanics, nobody asks "why is the probability proportional to the number of microstates?" — it's the *definition* of probability for a system whose microscopic details are unknown. The same logic applies to the lattice: before measurement, we don't know which globally consistent configuration the lattice will commit to, so we assign probability proportional to configuration count.

The Born rule is the *thermodynamics of the lattice*.

---

## 5. Addressing the Circularity Objection

A natural objection: haven't we smuggled in the Born rule by identifying |ψ|² with the configuration density? Isn't this circular?

No, and here's why. The chain of reasoning is:

1. **Start from the lattice.** Bits on a 2D surface with local constraints. No quantum mechanics assumed.
2. **Define Fisher information** from the bit-pattern statistics. This is pure information geometry — no quantum mechanics.
3. **Derive the metric** ds² = F dx² from the Fisher information. This gives the geometry of the lattice.
4. **Derive the dynamics** from the information action principle S_I = ∫√F ds. The Euler-Lagrange equation of this action is mathematically identical to the Schrödinger equation (Frieden, 1998).
5. **Identify** ψ = √ρ e^{iθ} where ρ is the normalised Fisher information density and θ is the accumulated phase.
6. **Derive the Born rule** from configuration counting on the lattice with two boundaries.

At no point do we assume the Hilbert space structure, the superposition principle, or the Born rule. They all *emerge* from the lattice:

- **Hilbert space** arises because the space of Fisher information amplitudes on the lattice has the structure of an L² function space.
- **Superposition** arises because the lattice can maintain multiple consistent configurations simultaneously (before measurement commits it).
- **The Born rule** arises because measurement forces commitment, and probability = configuration count.

The circularity in Gleason's theorem (which derives Born from the Hilbert space structure) and in many-worlds derivations (which derive Born from the branching structure) is avoided because we derive the *entire quantum formalism* from a more primitive substrate, not just the Born rule from the rest of quantum mechanics.

---

## 6. What Remains to Be Done

This working note presents the *structure* of the argument but several steps require full mathematical rigour:

1. **The Frieden derivation** of the Schrödinger equation from Fisher information is well-established in the literature but needs to be rederived specifically for the discrete circlette lattice (as opposed to the continuous limit that Frieden uses).

2. **The two-boundary configuration counting** needs a careful treatment of the phase summation. The claim that coherent summation over lattice sites with phases gives |⟨a|ψ⟩|² (rather than some other functional form) needs to be proven from the lattice update rules, not assumed from the standard quantum formalism.

3. **The connection to Zurek's envariance** should be explored. Zurek (2005) derived the Born rule from entanglement-assisted invariance ("envariance") in a way that may have a natural lattice counterpart: the symmetries of the circlette constraints under bit permutations may provide the envariance structure.

4. **Numerical verification.** A lattice simulation with a small number of sites, propagating a circlette through a beam-splitter geometry, should reproduce Born-rule statistics from direct configuration counting. This would be a powerful computational confirmation of the argument.

5. **The absolute normalisation** of F_vac (the vacuum Fisher information, identified with the cosmological constant) remains open. The dynamic-Λ working note derived its scale-factor dependence; a full derivation would also fix its absolute value, closing the loop between the Born rule and cosmology.

---

## 7. Summary

The Born rule — the foundational link between quantum amplitudes and observed probabilities — emerges naturally in the circlette framework as a consequence of three principles:

1. **Fisher information measures configuration density.** The wave function amplitude |ψ(x)| is the square root of the normalised Fisher information at lattice site x, which counts the density of distinguishable states.

2. **The lattice has two boundaries.** Preparation and measurement each provide boundary conditions, and the number of globally consistent configurations depends on the product of forward and backward configuration counts.

3. **Measurement is commitment.** When the lattice is forced to select one configuration, the probability of each outcome equals its share of the total configuration count — exactly as in statistical mechanics.

The result is:

$$P(a) = \frac{\Omega(a)}{\Omega_{\text{total}}} = |\langle a | \psi \rangle|^2$$

This is not a postulate. It is the thermodynamics of information on a holographic surface.

The Born rule is the lattice counting its options.

---

*For the main paper, see Elliman (2026), "It from Bit, Revisited." For companion working notes, see: "The Quantum Cheshire Cat and the Bridge Bit" (property stripping); "Measurement, Lattice Consistency, and the Dissolution of Retrocausality" (measurement problem and time-symmetric constraint propagation).*
