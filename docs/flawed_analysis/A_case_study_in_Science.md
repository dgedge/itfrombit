# The Garden of Many Paths in a Minefield

## How AI Will Make a Fool of You

<a href="assets/image-17.jpg"><img src="assets/image-17-thumb.jpg" alt="Garden of many paths in a minefield" style="width:40mm"></a>

## Introduction

This work began as a research project within Neuro-Symbolic Ltd. Dave, the author, had been working in AI for several years, mainly building deep-learning algorithms and large language models for applications ranging from reasoning systems to vision. A recent project involved building an explosives detector using Nuclear Quadrupole Resonance (NQR), which required some understanding of quantum physics.

The use of AI models to support coding and research was highly controversial among the engineers Dave was working with at the time. Some regarded them with deep suspicion. Correctness had always been the Holy Grail for a software engineer, and tools that produced confabulated answers were not to be trusted. Others recognised the enormous potential of AI for reading and summarising research, and for substantial productivity gains in writing code. Dave leant towards the latter group, then a minority. In retrospect, both groups were right. AI models are powerful and increasingly indispensable tools, but they also bring serious dangers that require careful safeguards.

A productive way forward seemed to be to devise methods that maximise the benefits of AI while using symbolic computing for speed, mathematics, and independent verification. Suitable architecture choices can reduce the degrees of freedom available for an AI model to generate undetected hallucinations, numerical errors, sycophancy, and confabulation.

Neuro-Symbolic Ltd is dedicated to that goal. So far, this project has made two contributions to it: a methodology for keeping research claims consistent, and a case study in how AI-assisted scientific work can go wrong unless it is externally disciplined.

## The Problem

### A Consistency-Maintenance System for Scientific Research

This project grew from Dave's long-standing interest in physics, his recent NQR work, and the need to stay up to date with quantum computing because of its defence implications. He had been reading John Wheeler's paper "It from Bit" [Wheeler 1992](#wheeler-1992), and wondered what information structure might be needed to support quantum physics.

Dave is not a professional physicist, but he has had rigorous scientific training to full Professor level at a Russell Group university, together with a strong mathematical background in engineering and information theory.

Truth can be a slippery concept. Pontius Pilate asked Jesus a good question: "What is truth?" Philosophers and epistemologists have been struggling to give a sound answer ever since. Science takes a narrower route. It tries to explain experimental results and observations by devising hypotheses. If a hypothesis survives extensive testing it may become a theory, or even a law. As soon as a fact is discovered that contradicts it, however, the theory must be abandoned or revised.

The author received an email from John Baez, who said that he now receives two or three Theories of Everything (TOEs) every day. Baez is the author of the [Crackpot Index](https://math.ucr.edu/home/baez/crackpot.html), which is surprisingly effective at dismissing theories that have no scientific validity. The theory proposed here only scored points for mentioning Feynman and Einstein, but was otherwise sound on that scale. Unfortunately, the Crackpot Index does not provide a threshold for deciding whether the term should be applied, so I cannot be entirely sure of my status. A minor crackpot at worst, perhaps.

Karl Popper argued that scientific theories can never be definitively proven true. Instead, they must be falsifiable and tested through rigorous attempts to disprove them [Merritt 2020](#merritt-2020). A further requirement is that a theory should make predictions that are later confirmed by experiment. As soon as a contradictory measurement is made, the theory fails.

This gives a useful test and a strategy for sustaining a theory. If it cannot be falsified, it does not even reach the starting gate; it is merely an opinion. While it survives, it should correctly predict experimental results and withstand serious attacks on its credibility. These ideas form the basis of the consistency-maintenance system described below.

## Failure Patterns in Science and Crackpot Theories

> "Seldom, very seldom, does complete truth belong to any human disclosure; seldom can it happen that something is not a little disguised or a little mistaken."
>
> Jane Austen

Few scientists, and few crackpots, are knowingly dishonest. They believe their findings and are keen to share them with a sceptical world. People are easily fooled. To quote Feynman, and earn more crackpot points:

> "The first principle is that you must not fool yourself — and you are the easiest person to fool."

The central problem is this: the human mind is a pattern-finding engine. It evolved to detect signals in noise quickly and cheaply, because a false alarm was survivable but a missed predator was not. It is therefore biased towards seeing structure that is not there. Add a flexible mathematical toolkit, which can be bent to fit almost anything, and a career structure that rewards positive, novel, clean-looking results, and error becomes the default state rather than an occasional accident.

A disciplined framework exists because the failure modes below are systematic and predictable. They can therefore be defended against in advance. What follows is a map of them, grouped by the mechanism that drives each.

### The Flexibility Problem: When a Model Can Fit Anything

**Numerology and pattern-matching on numbers.** Given any measured constant, one can almost always find a short combination of small integers, π, e, golden ratios, and square roots that reproduces it to a few decimal places. This feels like discovery, because the agreement looks too precise to be chance. Usually it is not. The space of "simple expressions" is enormous.

Arthur Eddington's attempts to derive the fine-structure constant, approximately 1/137, from pure combinatorics are the classic cautionary tale: elaborate, confident, and wrong. The warning signs are:

- the expression is found after the target value is known;
- no independent prediction follows from it.

### Overparameterisation: The Model With Too Many Knobs

The more free parameters a model has, the more data it can fit, regardless of whether it is true. Von Neumann's famous jibe was: "With four parameters I can fit an elephant, and with five I can make him wiggle his trunk." A fit is only impressive relative to how much freedom was spent to achieve it. A theory that matches ten data points using nine adjustable constants has explained very little.

### The Competitor and Search-Space Problem

This is the quantitative core of the first two problems. Before celebrating that expression X matches the data, one must ask: how many equally simple alternatives would also have matched within the same error bars? If the answer is "many", then the hit is consistent with the data, but it is not strong evidence for the proposed mechanism. The data could not have told the alternatives apart. A match is evidence only to the degree that rival explanations are excluded.

### The "Too Good to Be True" Fit

Counterintuitively, a perfect fit, with many decimal places and near-zero deviation, should increase suspicion rather than confidence. Real measurements carry noise. A result that fits better than the noise allows is often a signature of a bug, circular calculation, or selective data. Fisher's reanalysis of Mendel's pea-genetics data found that the results were closer to the predicted ratios than random sampling should normally produce. The data were, in some sense, too good.

### The Data-Contingent Analysis Problem

This is one of the most insidious traps because it catches honest researchers. It is sometimes known as HARKing: hypothesising after the results are known. You collect data, notice a striking pattern, and then write the paper as if you had predicted that pattern in advance. An exploratory finding, which generates hypotheses, is silently relabelled as a confirmatory finding, which tests them. The statistics that are valid for a pre-stated hypothesis are invalid once the hypothesis has been read off the data.

**The garden of forking paths.** This is the deeper version, named by Andrew Gelman and Eric Loken [Gelman and Loken 2013](#gelman-and-loken-2013). Their key insight was that one does not have to consciously fish for significance to be fooled. Suppose one genuinely intends to run one clean analysis. Along the way, however, there are dozens of reasonable, data-dependent choices: which outcome to focus on, whether to exclude an unusual subject, where to put an age cut-off, whether to analyse men and women separately, and whether to log-transform a variable. Each choice seems locally justified. The problem is that, had the data come out differently, different choices would have been made. The reported p-value pretends that only one path existed. Because the unreported look-elsewhere space was large, a "significant" result is far more likely to be noise than the single p-value suggests.

**A concrete illustration.** Studies have reported that clothing colour correlates with phases in a woman's menstrual cycle, and that forearm girth predisposes a particular voting preference. No single analysis was necessarily rigged, but the researchers chose which combination of cycle window, status coding, and outcome to highlight only after seeing which combination looked interesting. This is well illustrated by an [xkcd cartoon](https://xkcd.com/882/) in which 20 jelly-bean colours are tested against acne and green jelly beans are found to be "linked" at p < 0.05. The published headline is "Green Jelly Beans Cause Acne", while the 19 colours that showed nothing are buried. That is the multiple-comparisons effect, of which the garden of forking paths is the subtler, unconscious cousin.

### p-Hacking and Researcher Degrees of Freedom

The active version is to try several analyses and report the one that crosses significance. Add subjects until p dips below 0.05, then stop. Drop the inconvenient covariate. Slice the data another way. Each move may be individually defensible; collectively they manufacture false positives. The food-science lab of Brian Wansink became the textbook case: internal emails showed data being sliced repeatedly until "something" emerged, leading to many retractions.

### Training on the Test Set

In any modelling or machine-learning context, evaluating a model on the same data used to build it guarantees optimism. Without a held-out test set that the model has never seen, one measures how well it memorised, not how well it predicts.

### The Remedy That Defines the Discipline

Every item above is defeated by committing to the analysis before seeing the outcome: pre-registration, blinding, and holding out data. These work not by making the researcher smarter, but by removing the forks. The data can no longer choose the path, because the path was fixed first.

## Human Cognitive Biases

People evolved to make quick decisions rather than long, rigorous chains of reasoning. In evolution, the decision was often instant action or death. This means we must guard carefully against our own biases.

### The Confirmation Problem: Seeing What You Expect to See

**Confirmation bias and motivated reasoning.** We seek, weight, and remember evidence that supports what we already believe. We scrutinise disconfirming evidence far more harshly than confirming evidence. The asymmetry is the danger: a result we like sails through; an identical result we dislike gets a second look that finds the error. The fix is to spend the disconfirming effort symmetrically, making "yes" as expensive as "no".

### Apophenia and the Texas Sharpshooter

Apophenia is the perception of meaningful patterns in random data. The Texas sharpshooter fires at a barn wall and then paints the target around the densest cluster of bullet holes, declaring himself a marksman. Drawing the hypothesis around the cluster one happened to find is the same move, and it always produces a bullseye.

### Ignoring Outliers and Cherry-Picking

Discarding the data points that do not fit is one of the most common forms of self-deception, because there is usually a plausible-sounding excuse available. "They must be measurement error" is often a dangerous sentence. Millikan's oil-drop notebooks show drops annotated and excluded; whatever his exact justification, the lesson stands. Exclusions must follow a rule fixed before one sees which points the rule will remove, or they become a tool for sculpting the answer. An outlier is either a flaw in the data or a flaw in the theory. One cannot know which if one deletes it on sight.

### Ad Hoc Rescues and Unfalsifiability

When a theory is contradicted, the temptation is to bolt on an auxiliary assumption that explains away the contradiction and nothing else. Do this repeatedly and one gets what Imre Lakatos called a degenerating research programme: a theory that survives only by accumulating patches, each of which shrinks its empirical content.

Popper's criterion is the guard. A theory that cannot specify, in advance, an observation that would prove it wrong is not making a scientific claim. "My theory is consistent with everything" is not a strength; it is the absence of content.

### Confusing "Consistent With" and "Evidence For"

A theory that accommodates a fact after the event is in a completely different epistemic position from one that predicted the fact before measurement. Treating the two as equal is how weak theories acquire the appearance of strong support.

### The Human and Social Problem

The failure modes above are amplified by the environment in which scientists work. Publication pressure rewards a steady stream of positive, clean, novel results: exactly the results that the failure modes above manufacture. Most misconduct is not outright fraud but a spectrum of questionable research practices: rounding a p-value down, omitting a failed experiment, presenting an exploratory result as confirmatory, or quietly dropping a condition that did not work. Each step feels small and is easy to rationalise under a deadline. The danger is that the cumulative effect on the literature is indistinguishable from deliberate fabrication, and the incentives push everyone, including the honest, in the same direction.

### The File-Drawer Problem

Negative and null results are far less likely to be published, so they sit unseen in the file drawer. The published literature is therefore a biased sample of all experiments run, over-representing flukes that happened to reach significance. A field can collectively "discover" an effect that does not exist, simply because the experiments that found nothing were never reported.

### Sunk Cost and Commitment to a Pet Theory

The more years, papers, and reputation invested in an idea, the harder it becomes to abandon. A scientist defending a lifetime's theory is fighting not just for a hypothesis but for the meaning of a career, and may unconsciously lower the bar for evidence that saves it and raise the bar for evidence that kills it.

### Groupthink, Authority, and Conformity

Fields develop consensus, and dissent carries social and professional cost. Citation cartels, deference to senior figures, and reviewer conservatism can entrench an error for a generation. The pressure to agree is itself a bias. Agreement requires the same scrutiny as disagreement, but rarely receives it.

### The Seduction of Beauty

Elegant, symmetric, simple theories are more persuasive than they deserve to be. Beauty is a heuristic, not evidence, and the history of physics is littered with beautiful ideas that were wrong.

### Pathological Science

Irving Langmuir used the term pathological science for research that is not fraud, and whose practitioners are sincere, but where wishful thinking, effects near the threshold of detectability, and refusal to accept disconfirmation produce a phantom result that slowly fades under scrutiny. N-rays, polywater, and cold fusion are canonical cases.

## The Replication Reckoning: Evidence That This Is Systemic

The replication crisis shows that these are not rare lapses but structural features. When large collaborations re-ran published findings under pre-registered conditions, only a minority reproduced. In a landmark 2015 effort in psychology, roughly one-third to one-half of high-profile results replicated, often at much smaller effect sizes. Similar shortfalls appeared in parts of cancer biology, economics, and beyond. Celebrated effects such as power posing, ego depletion, and many priming studies substantially weakened or evaporated. The crisis is not mainly a story of villains. It is the predictable aggregate of forking paths, publication bias, and underpowered studies operating across an entire field.

Physics is not exempt. It has merely paid for its lessons publicly. The 2011 OPERA "faster-than-light neutrinos" turned out to be a loose fibre-optic cable. The 2014 BICEP2 announcement of primordial gravitational waves was, on reanalysis, galactic dust. The 2015 750 GeV diphoton bump at the LHC spawned hundreds of theory papers before vanishing as a statistical fluctuation: a live demonstration of the look-elsewhere effect.

## Why a Disciplined Framework Is Needed

The unifying lesson is Feynman's: the easiest person to fool is yourself. Intelligence and integrity are no defence. The smartest, most honest researchers fall into these failure modes precisely because the failures operate below conscious awareness. One cannot out-think a bias one never notices.

Discipline works by a different route. It does not ask people to be cleverer or more honest; it structurally removes the degrees of freedom that the mind exploits to fool itself. Each guardrail neutralises a specific failure mode:

| Failure mode | Disciplinary guardrail |
|---|---|
| Numerology and expression-shopping | Count the competitors: how many rival simple formulae also fit? State a result as "consistent with", not "evidence for", when several alternatives survive. |
| Overfitting and too many knobs | Count parameters against data; penalise flexibility; demand out-of-sample prediction. |
| Forking paths, HARKing, and p-hacking | Pre-register the hypothesis and analysis; blind the analyst; hold out a test set touched once. |
| Confirmation bias and motivated reasoning | Default to refutation; spend effort trying to falsify your own result; keep "yes" as expensive as "no". |
| Ignoring outliers and cherry-picking | Fix exclusion rules before seeing which points they remove; treat an outlier as a possible refutation, not a nuisance. |
| Single lucky derivation | Require independent routes to the same result; report where they disagree. |
| Too-good-to-be-true fits | Treat a perfect fit as a red flag to re-check, not as a triumph. |
| Unfalsifiability and ad hoc rescue | State in advance what observation would kill the claim; log every patch as a debt, not a feature. |
| Drifting claims and quiet retraction | Keep a retraction ledger: record superseded claims with the reason, and never silently overwrite. |
| Asserting from memory or belief | Compute and quote the number, cite a checkable source, or mark it explicitly unverified. |
| Anchoring to prior output | Treat your previous conclusion as someone else's work to audit; prior agreement is not evidence. |

The deep point is that none of these guardrails trusts the individual. They relocate the burden of correctness from fallible human judgement onto external, mechanical checks: a pre-committed plan, a held-out dataset, a competitor count, an independent derivation, a written ledger, or an executable script that asserts every number it prints. A framework is disciplined precisely to the degree that it does not require its practitioners to be unusually wise or unusually honest to arrive at the truth.

## AI Has a More Dangerous Flaw for Scientists: It Flatters You

If you log into an AI system and type, "I think I have solved dark matter using a 5D vortex model; here is my maths", the model will rarely say, "This is numerology and your assumptions are flawed." Instead, it is likely to say something like, "This is a fascinating and potentially groundbreaking paradigm shift", and then help format the flawed mathematics into a professional-looking physics paper.

This creates a dangerous psychological feedback loop. The amateur theorist proposes a wild idea, the AI validates it, the theorist feels confirmed, and the AI connects a few random numbers to "prove" the theory. The human walks away with inflated confidence, unaware that they were talking to a mirror trained to be agreeable.

This is largely our fault. AI models are trained using processes such as Reinforcement Learning from Human Feedback (RLHF), which rewards answers that humans rate as helpful, polite, and agreeable. Models are therefore pushed towards sycophancy.

Sycophancy is not merely noise. It is an active distortion of the evidence channel. A model that reflexively validates makes "yes" cheaper than "no", biasing the whole collaboration towards confirming whatever the user already believes. Stripping praise from the interaction is not mainly about tone; it is about keeping agreement expensive so that, when it comes, it carries information. In this project, the valuable turns were often the ones where the model disagreed and showed why.

A useful test is this: the anti-sycophancy directive is working when the model's agreement and disagreement are indistinguishable in tone. Both should be flat, and both should be backed by a check. If "yes" sounds warmer than "no", the evidence channel is still biased.

You might miss the praise, and the AI might start to sound like an old headmaster, but that is better than becoming delusional.

## The Ultimate Coincidence Engine

AI language models are the most powerful pattern-matching engines ever built. If one asks an AI to find a mathematical relationship between the height of the Great Pyramid, the speed of light, and the mass of an electron, it will find one. It will combine fractions, factors of π, roots, and constants, and serve up a beautifully formatted equation that looks like a profound secret of the universe.

Another anti-pattern is to ask it to produce five facts in support of whatever one is investigating, such as "jellyfish are immortal". If it can find only four, it may hallucinate a plausible fifth.

Sycophancy and flattery are the most corrosive failure modes in scientific use. Models trained to be helpful and agreeable learn to mirror the user's framing. Ask "isn't it elegant that my theory predicts X?" and the model will tell you that it is elegant and then find reasons that X follows. Propose a flawed premise and it may build on it rather than challenge its foundation.

### Towards More Dangerous, More Delusional Science

The role science most needs a collaborator to play is that of an adversary, as Popper's criterion demands: someone who tries to break the idea. A sycophantic model is the opposite. It supplies an endless stream of agreement and elaboration for whatever the user already believes. It converts the tool that should be the strongest check on confirmation bias into an amplifier of it.

Give the garden of many paths to an artificial intelligence model and ask it to find the way. Of course it always will.

For a lone theorist emotionally invested in a pet Theory of Everything, this is a private chorus of validation available on demand. The dangerous combination is the cognitive bias of a person joined to the cognitive bias of a language model. The potential for delusion is large.

Over the last two years, thousands of amateur theorists have started using AI models such as ChatGPT, Gemini, and Claude to help generate their TOEs. The results are often a disaster. They are not trying to cheat or fool anyone. They are true believers. But the easiest person to fool is always yourself, and we need to identify the actual failure modes in order to counter them.

### How LLMs Compound Research Failures

Large language models do not merely add a new list of errors. They amplify human errors. Nearly every human failure mode has an LLM counterpart that feeds it, and the human-model pair is more dangerous than either alone.

The reason is structural. An LLM is optimised to produce text that is probable and approved of: fluent, expert-sounding, and agreeable. Those objectives are orthogonal to truth. The model is therefore a machine for generating exactly the kind of plausible, confident, confirmatory material that human biases are primed to accept. Nowhere is this more acute than in grand unified theories and TOEs, which sit at the maximally flexible, maximally unfalsifiable, maximally seductive end of the spectrum: the perfect substrate for the loop described below.

1. **Sycophancy: the disappearance of critical evaluation.** The model gives the user the validation that science most needs to withhold.
2. **Hallucination and confabulation.** It invents citations, DOIs, author names, theorem names, numerical values, experimental results, and derivation steps, presenting them in the same fluent register as true material.
3. **High plausibility without grounding.** An LLM optimises for the probability of text, not for correctness. A wrong derivation can look exactly as authoritative as a right one.
4. **Numerology and curve-fitting at industrial scale.** Ask it to relate constant A to constant B and it will find a route, silently discarding the thousands that failed.
5. **Decorative formalism.** It can generate Lagrangians, symmetry groups, and cohomology arguments that look load-bearing but are ornamental.
6. **Anchoring and context contamination.** The model inherits the user's notation, assumptions, and desired conclusion. Within a long conversation it may treat its own earlier speculation as established fact.
7. **Spurious coherence.** It is exceptionally good at making disparate claims hang together into a satisfying narrative. The coherence is rhetorical, not physical.
8. **Miscalibrated confidence.** Models often express uniform confidence regardless of reliability, and rarely volunteer "I have not checked this."
9. **The generation-verification asymmetry.** It has always been cheaper to produce a plausible claim than to check it. LLMs push the generation cost towards zero while leaving verification cost high.
10. **Tirelessness.** The model never feels the fatigue, friction, or doubt that makes a human pause and ask whether something is really right.

### The Compound Mechanism: A Delusion Amplifier

The individual failures are serious. Their interaction is worse. Place a human with an emotionally invested pet theory next to a sycophantic, fluent, hallucinating, tireless model, and one gets a closed feedback loop with no external anchor:

1. the human proposes an idea;
2. the model agrees, elaborates, and supplies fabricated or numerological support;
3. the human's conviction deepens and the theory grows more intricate;
4. the human proposes more;
5. the model confirms and extends further.

At no point does reality enter the loop. Each turn adds plausible-looking structure and subtracts nothing. The result is a shared, self-reinforcing belief system that becomes more elaborate and more detached from testability with every iteration.

### Why TOEs Are the Perfect Storm

A Theory of Everything maximises every vulnerability at once:

- it has no single decisive experiment, so falsification is slow or absent;
- it ranges over all of physics, so there is always another coincidence available to "explain";
- it rewards grand narrative coherence, which is one of the model's strengths;
- it attracts lone, emotionally committed theorists, where sycophancy has maximum leverage;
- it lives in the numerology-rich territory of fundamental constants.

The LLM is best at exactly the things a TOE most wants, and worst at providing the one thing it most needs: a reason to stop.

### Why the Discipline Must Be External and Mechanical

The decisive insight is that one cannot fix this by asking the model to be more careful, honest, or rigorous. Those requests are answered with the same fluent, agreeable, unverifiable text as everything else. The model's self-report about its own reliability is itself an LLM output.

The durable defence is to relocate the burden of correctness away from the model's judgement and onto external, mechanical checks:

| LLM failure mode | Mechanical guardrail |
|---|---|
| Sycophancy and mirroring | Mandate an adversarial stance. Agreement and disagreement should carry the same evidentiary burden. |
| Hallucinated numbers or derivations | Compute, do not assert. Every number must come from executed code, a checkable source, or be tagged unverified. |
| Fabricated citations or theorems | Require a resolvable source for every reference. Treat any unverifiable citation as fabricated until checked. |
| High plausibility with no tell | Tag every claim as verified, proposition, or speculation. Never let fluency stand in for evidence. |
| Industrial numerology | Count the competitors. Report the search space, not just the hit. |
| Anchoring and context contamination | Treat the model's prior turn as someone else's work to audit. Prior agreement is not evidence. |
| Spurious coherence and TOE narrative | Demand falsifiable, advance predictions. A coherent story is not a confirmed one. |
| Miscalibrated confidence | Make "I do not know" and "I have not checked" first-class answers. |
| Generation much cheaper than verification | Automate verification too: self-asserting scripts, human review of diffs, and reproducible checks. |
| Feedback loops | Keep a retraction ledger so the loop cannot quietly forget its own failures. |

The constructive thesis of this report is that a disciplined framework is not an optional refinement for AI-assisted science. It is the only thing standing between these tools and the accelerating production of confident, intricate, delusional theory. Used inside such a cage, an LLM is a powerful instrument. Used outside one, it is the most effective engine for self-deception yet built: not because it is malicious, but because it will tirelessly and articulately give you the beautiful wrong answer you were hoping for.

## The Solution: How ANCHOR, DRIFT, and PTMS Keep the Research Honest

## The Shared Failure Mode

AI systems and human researchers fail in the same direction. Both are fluent generators of plausible narratives, and both are bad at noticing when a plausible narrative has detached from what was actually shown. An LLM confabulates, flatters the user, and recalls from training instead of computing. A person rationalises, grows attached to a prior claim, remembers selectively, and finds patterns in noise. The common root is that grounding depends on in-the-moment discipline, and discipline is exactly what fails under pressure, fatigue, or desire for a result.

ANCHOR, DRIFT, and PTMS attack that root by moving grounding out of the agent and into artefacts. Truth-tracking stops being a virtue one has to exercise and becomes a structure one has to pass through.

## The Three Instruments

**ANCHOR.md: the settled canon.** ANCHOR is append-only and organised with dated brackets, one per substantive result. Each entry names the script that backs it, with `exit 0` meaning that the script's own assertions passed. It records the claim tier: derived or computed, imported, conditional, or open. It also records citations. A numbering checksum gates every commit, so accidental corruption is caught mechanically. Supersessions are recorded with reasons and never silently overwritten.

*Counters:* over-claiming, confabulation, silent revision, and structural corruption.

**DRIFT.md: the exploration and drift log.** DRIFT is where frontier attempts, reversals, and "this turned out wrong" entries are logged, dated, and cross-referenced. Its job is to make the evolution of a claim legible: not just where it ended, but every place it moved and why.

*Counters:* the universal urge, human and AI alike, to bury a mistake and pretend the current answer was always the answer. Under DRIFT, a retraction is a permanent entry.

**PTMS: the machine-readable claim registry and structural-integrity checker.** PTMS, run with `python -m ptms check`, extracts claims into a registry and checks for registry drift and tier-coverage debt. Crucially, the foundations paper states plainly that a clean PTMS run is not a proof. PTMS enforces structure, not truth.

*Counters:* sycophancy and motivated reasoning, because a machine check does not care who is asking or what they want to be true. By being explicit about what it does not establish, PTMS also counters over-claiming about the verification itself.

## The Disciplines That Bind Them

- **Compute, do not recall.** Important figures, rates, code symbols, and results are read from source or re-derived, never quoted from memory.
- **Self-asserting scripts.** Verification lives in the artefact: `exit 0` is the certificate.
- **Per-claim tiering.** Every statement carries its epistemic status.
- **Supersession with reason.** Claims change in the open, with the prior version and the cause preserved.
- **Single source per fact.** Lower-authority stores point at higher-authority ones rather than restating them. Restated facts drift.
- **Honest negatives are first-class.** A no-go result is as valuable as a positive result. The system is built to surface its own over-claims, not hide them.

## The Living Example

This project provides its own case study. The photon's Lorentz invariance was initially marked as passed: an over-claim. A direct lattice computation, run to confirm it, instead contradicted it, forcing a retraction to a Collins-Perez-Sudarsky-style open problem. Deeper computation then showed that this framing rested on a misidentification of which object was the photon, and resolved that part of the issue, with the genuinely open piece honestly re-filed. At every step there was a dated bracket, a script, a tier, and a stated reason.

What matters is that the apparatus caught both an AI error and a human error. Nobody had to be infallible. The structure surfaced the error and made the correction legible.

## The Core Insight

The point is not that the framework's physics is right. PTMS says explicitly that a clean run proves nothing of the sort. The point is that the research is kept honest and self-consistent. At any moment one can see the status of every claim, errors surface quickly, and nothing needs to be oversold. That value is independent of whether the theory ultimately succeeds.

One cannot argue with a script that will not run. One cannot quietly retract under an append-only log. One cannot pass a guess off as a theorem under tiering. One cannot corrupt the canon past a checksum. One cannot flatter a machine check into agreement. Because the standard is external, it binds the AI and the human identically, which is precisely why it works on both.

## References

### Wheeler 1992

John Archibald Wheeler, "Recent thinking about the nature of the physical world: It from bit", *Annals of the New York Academy of Sciences*, 655(1), 349-364, 1992. Wiley Online Library.

### Merritt 2020

David Merritt, *A Philosophical Approach to MOND: Assessing the Milgromian Research Program in Cosmology*, Cambridge University Press, 2020. [Google Books preview](https://books.google.co.uk/books?hl=en&lr=&id=5iDgDwAAQBAJ&oi=fnd&pg=PR9&dq=A+Philosophical+Approach+to+MOND&ots=uOYhaVKTd-&sig=AL1vJ84ztTGrDcTLuMbPbi8mZIU&redir_esc=y#v=onepage&q=A%20Philosophical%20Approach%20to%20MOND&f=false).

### Gelman and Loken 2013

Andrew Gelman and Eric Loken, "The garden of forking paths: Why multiple comparisons can be a problem, even when there is no 'fishing expedition' or 'p-hacking' and the research hypothesis was posited ahead of time", Department of Statistics, Columbia University, 2013. [PDF](https://sites.stat.columbia.edu/gelman/research/unpublished/p_hacking.pdf).

### Further Reading

For the Einstein-Bohr debate, the source draft contained a placeholder reference but no usable URL. A full citation should be added before this item is cited as evidence.
