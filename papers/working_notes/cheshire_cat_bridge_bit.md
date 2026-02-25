# The Quantum Cheshire Cat and the Bridge Bit

## How the Circlette Architecture Enables "Property Stripping"

**D.G. Elliman — Neuro-Symbolic Ltd, February 2026**

Working note for *It from Bit, Revisited*

---

## 1. The Phenomenon

In the Quantum Cheshire Cat (QCC) effect, a particle's physical properties — spin, charge, magnetic moment — are measured at a *different spatial location* from the particle itself. The "cat" is in Path A; its "grin" (spin) is in Path B.

The standard quantum-mechanical explanation treats this as a consequence of weak measurement and post-selection. The circlette framework offers a *structural* explanation: it tells us *why* properties are separable from particles in the first place.

## 2. The Ring Architecture

Recall the circlette layout:

```
G₀ — G₁ — C₀ — C₁ — LQ — I₃ — χ — W
 ↑___________________________________|
              (ring closes)
```

The ring partitions into three sectors:

| Sector | Bits | Encodes |
|--------|------|---------|
| Generation | G₀, G₁ | Flavour (1st, 2nd, 3rd) |
| Colour | C₀, C₁ | SU(3) colour charge |
| Electroweak | I₃, χ, W | Isospin, chirality, weak participation |
| **Bridge** | **LQ** | **Lepton/quark identity** |

Electric charge is *not* a single bit — it is a *derived quantity*: Q = f(LQ, I₃). This is the key to the Cheshire Cat.

## 3. The Bridge Bit as Property Router

The bridge bit LQ occupies a unique architectural position:

1. **It is the CNOT control line.** The unique update rule is I₃(t+1) = I₃(t) ⊕ LQ(t). The bridge bit *controls* the dynamics of the electroweak sector without itself being modified.

2. **It connects two independent sectors.** Colour (C₀, C₁) lies on one side; electroweak (I₃, χ, W) lies on the other. LQ is the *only* bit that participates in constraints spanning both sectors.

3. **It determines which sectors are "active".** For leptons (LQ = 0), the colour sector is forced to (0,0) — colourless. For quarks (LQ = 1), all three colour states are available. The bridge bit acts as a *routing switch*.

### The critical point

In the circlette model, a "particle" is not a marble with properties glued to it. A particle is a *pattern of bits on a ring*, and different bits encode different properties. The ring exists on a 2D holographic lattice; the 3D world we observe is the *error-corrected logical content* of that lattice.

This means:

- The **"location" of the particle** is the lattice position of the ring pattern as a whole — the coherent excitation that propagates across the surface.
- The **"charge" of the particle** is the value of the derived function Q = f(LQ, I₃) — specific bits within the ring.
- The **"spin" of the particle** is encoded in the circulation direction and the chirality bit χ.

These are *logically independent addresses* on the lattice. The pointer to "where the ring is" and the pointer to "what bits 4 and 5 say" are separate pieces of lattice bookkeeping.

## 4. The Cheshire Cat Mechanism

In an interferometer, the beam splitter creates a superposition of the ring pattern across two lattice paths. But superposition on the holographic lattice is not "the whole ring goes left OR the whole ring goes right." It is:

**The lattice creates a superposition of *sector-level* amplitudes across paths.**

Here is how the bridge bit enables this:

### Step 1: Sector decomposition at the splitter

When a circlette enters the interferometer, the beam splitter couples to the *propagation mode* of the ring — the coherent lattice excitation that carries the pattern forward. But because the ring has internal structure with sectors separated by the bridge bit, the splitter can create *different superposition coefficients for different sectors*.

Specifically, if the initial state is |ψ⟩ = |ring, Path A⟩, the splitter produces:

```
|ψ_split⟩ = α|propagation, A⟩ ⊗ |internal, A⟩ + β|propagation, B⟩ ⊗ |internal, B⟩
```

But "internal" decomposes across the bridge:

```
|internal⟩ = |generation⟩ ⊗ |colour⟩ ⊗ |LQ⟩ ⊗ |electroweak⟩
```

The bridge bit LQ acts as an *entanglement barrier*: sectors on opposite sides of the bridge can acquire independent path amplitudes because they are coupled only *through* LQ, not directly.

### Step 2: Post-selection separates the pointers

With appropriate post-selection (choosing only outcomes where the particle exits through a specific port), the weak measurement statistics reveal that:

- The **ring pattern** (identified by its propagation mode) has amplitude concentrated on Path A.
- The **spin/magnetic moment** (identified by the electroweak sector, specifically χ and I₃) has amplitude concentrated on Path B.

This is possible because the bridge bit *decouples* the sectors during lattice propagation. The generation+colour sectors carry the "identity" of the ring (which particle it is), while the electroweak sector carries the "spin character." The bridge bit routes these independently.

### Step 3: The CNOT structure makes this natural

The CNOT rule I₃ ⊕ LQ operates *only at the bridge boundary*. During propagation between ticks, the sectors evolve quasi-independently. The bridge bit synchronises them at each tick, but between ticks — which is when the interferometer's path separation acts — the sectors can be spatially separated on the lattice.

This is precisely analogous to how, in a quantum error-correcting code, *logical* qubits can be delocalised across *physical* qubits. The circlette IS an error-correcting code. Its properties are logical bits that need not be co-located on the physical lattice.

## 5. Why This Is Not "Weird" in the Circlette Framework

In the standard particle ontology, a particle is a point-like entity that *carries* its properties. Separating the properties from the carrier is paradoxical — hence the "Cheshire Cat" name.

In the circlette ontology:

| Standard picture | Circlette picture |
|---|---|
| Particle is a thing | Particle is a pattern |
| Properties are attributes of the thing | Properties are sub-patterns (sectors) |
| Location = where the thing is | Location = where the pattern propagates |
| Charge = attribute at that location | Charge = derived function of specific bits |
| Separation is paradoxical | Separation is *expected* for delocalised codewords |

The Quantum Cheshire Cat is not a paradox in the circlette model. It is a *prediction*. Any system where "identity" and "properties" are encoded in different sectors of an error-correcting code, connected by a bridge, should exhibit property-stripping under weak measurement with post-selection.

## 6. Specific Predictions

The circlette model makes testable predictions beyond the standard QCC:

1. **Sector-specific stripping order.** Properties encoded in the electroweak sector (spin, isospin) should be easier to strip than properties encoded in the colour sector (colour charge), because the bridge bit separates them. The generation sector (flavour) should be hardest to strip, being furthest from the bridge on the ring.

2. **Bridge bit as the hinge.** If a Cheshire Cat experiment could probe the lepton/quark identity bit independently, the bridge bit should always co-locate with whichever sector has the *larger* path amplitude — it acts as a router, not a passenger.

3. **Multiple simultaneous strippings.** The 2023 photon experiments showing multiple properties stripped simultaneously are natural in this framework: different sectors can take different paths. The maximum number of independently separable property groups should equal the number of sectors (three, plus the bridge).

4. **No stripping for leptons' colour.** Since leptons have LQ = 0 and colour = (0,0), there is no colour property to strip. Only quarks should exhibit colour-path separation, and only in conditions where confinement is relaxed (extreme energies).

## 7. Connection to the Information Action Principle

The information action S_I = S_external + S_internal tells us that the cost of propagation (external) and the cost of internal dynamics (internal) are *separately minimised*. This factorisation is exactly what allows the Cheshire Cat: the lattice can route the external (propagation/location) and internal (properties/sectors) components along different minimum-action paths when the interferometer geometry makes this energetically favourable.

The bridge bit's role as a read-only CNOT control line means it contributes to S_internal but not to S_external — it doesn't cost the lattice anything to propagate LQ, because LQ never changes. This asymmetry between the bridge and other bits is what makes property routing possible without violating conservation laws.

---

## Summary

The Quantum Cheshire Cat is a natural consequence of the circlette architecture:

- Particles are **patterns**, not point objects.
- Properties are **sector-level sub-patterns**, connected by the bridge bit.
- The holographic lattice can **route sectors independently** during propagation.
- The CNOT update rule **decouples sectors between ticks**, enabling spatial separation.
- The information action **factorises**, allowing external and internal components to take different paths.

The bridge bit LQ is the architectural feature that makes this possible. It is the "neck" of the Cheshire Cat — the narrow connection through which the grin (electroweak properties) is attached to the cat (particle identity). Weak measurement with post-selection simply reveals that the neck is a *logical* connection, not a physical weld.

---

*For the main paper, see Elliman (2026), "It from Bit, Revisited." Available at [GitHub link].*
