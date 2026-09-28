# specs/039 — Interconversion cycles as conserved pools (design; not yet pre-registered)

**Status:** design, 2026-09-27. Awaiting Adam's decision on the fixed
assumption (φ₀) and on where detection lives, before pre-registration.

## Why

specs/034 §11–12:
- our loops are mostly curated interconversion cycles (A ⇄ B);
- the solver multiplies around them (gain 1, the knife-edge, the RAF collapse);
- MP-BioPath's curators cut such loops by hand and added "reverse enzyme ⊣
  modified form" inhibitions.

The catalog has **259 curated interconversion pairs in 55 pathways**.
- 93 have a catalyst in both directions, 150 in one, 16 in neither.
- Phospho 89, ubiquitin 59, GTP/GDP 32, binding/dissociation 31, methyl 28.

## Method: a blind double derivation (agreed with Adam)

The same problem statement (`problem.md`) went to two models independently:
- this session (Opus 5.5): `derivation-opus.md`, written before the other was
  read;
- a Fable 5.1 subagent: summarised below.

## The two derivations agree on the rule

A steady-state rate balance on a conserved pool, A + B = total:

    r = (u + ε) / (v + ε)       D = φ₀·r + (1 − φ₀)
    B_fold = s · r / D           A_fold = s / D
    F_fold = R_fold = J = s / ((1 − φ₀)/u + φ₀/v)      (the cycle flux)

- u is the forward reaction's drive, evaluated with the existing AND/divide
  semantics with A held at baseline.
- v is the same for the reverse reaction, with B held at baseline.
- s is the supply fold of the cycled protein: its pin, else the OR-mean of its
  external producers, else 1.
- φ₀ is the baseline fraction of the pool in the modified form.

**Multi-state cycles:** J = s / Σⱼ (φⱼ / uⱼ), and Aᵢ = J / uᵢ.

**Properties:**
- exactly baseline when unperturbed;
- no gain-1 loop and no all-zero root;
- bounded: B ≤ s/φ₀;
- the direction of every enzyme perturbation is correct for any φ₀ in (0, 1).

## They differ on φ₀, and the difference matters

| φ₀ | 80x through a 4-tier cascade (RAS → RAF → MEK → ERK) | source |
|---|---|---|
| 0.5 (Opus, uninformative) | 1.98, 1.33, 1.14, **1.07: ERK reads NORMAL** | the compression flagged in both derivations |
| **0.1** (Fable; modified forms a minority at rest: RAS ~5–10% GTP-loaded, basal phospho-occupancy ≤10–20%) | stays UP for about 20 tiers | knockouts are exact 0 at any depth, for any φ₀ |

**Proposed: φ₀ = 0.1, fixed a priori** on the biological justification.
- Orientation: the product of the catalysed, modifying direction is B.
- φ₀ = 0.5 would be a declared sensitivity arm, never used to choose.

Fable notes that the φ₀ → 0 limit (B = s·u/v, A = s) is exactly MP-BioPath's
hand cut plus its added inhibitions, with no parameter at all. That holds only
where the reverse reaction is curated (24% of the added inhibitions, specs/034
§12b).

## Added by the Fable derivation

- **Reaction nodes carry the flux J**, so co-outputs (ADP, Pi) read flux.
- **Drop the E_f ⊣ A depletion edge inside a cycle**; the rule already lowers A.
- **Autocatalysis** (B catalyses F): bistable. Take the nonzero root that is
  continuous from baseline.
- **Detection** by shared member leaves (forms of one reference entity).
  Degradation must not be mistaken for a reverse step.
- **Failure modes:**
  - φ₀ near 1 (constitutively modified sites);
  - unmodified-form readouts read NORMAL under an enzyme perturbation at
    φ₀ = 0.1 (1.11);
  - a set-member knockout gives u = (N−1)/N, so B reads DOWN only for N ≤ 6;
  - overlapping cycles need the Markov-chain form;
  - a pinned member fixes the other form by conservation.

## Fixtures (φ₀ = 0.1, s = 1)

| perturbation | A | B | F = R |
|---|---|---|---|
| none | 1 | 1 | 1 |
| E_f KO | 1.111 | 0 | 0 |
| E_f 80x | 0.1124 | 8.989 | 8.989 |
| E_b KO | 0 | 10 | 0 |
| E_b 80x | 1.1096 | 0.01387 | 1.1096 |
| gene 80x (KO) | 80 (0) | 80 (0) | 80 (0) |
| both enzymes 80x | 1 | 1 | 80 |
| two tiers, kinase₁ 80x | A₂ = 0.556 | B₂ = 4.997 | 4.997 |

## Open decisions (Adam)

1. φ₀ = 0.1 as the fixed assumption.
2. Detection in the generator, which ships `cycles.csv` (forward, reverse, A, B
   per virtual-reaction pair), rather than in the solver.

## Decisions (Adam, 2026-09-27)

1. **φ₀ = 0.1** for every cycle, fixed a priori. φ₀ = 0.5 runs only as a declared
   sensitivity arm and is never used to choose.
2. **Detection in the generator.** It detects each protein's pool once (its
   forms, by shared reference entity, and the curated reactions converting
   between them) and ships it with the network. DeltaSignal solves the pool on
   every run.
3. **General form: π = πP** (Adam's framing).
   - Each pool is a Markov chain: the states are the protein's forms, and each
     transition rate is its baseline rate × its drive.
   - The stationary distribution gives the split between forms.
   - The two-state closed form above is the special case.
   - For a branched pool, the baseline rates need one extra stated assumption:
     an equal split of flux between exits, or detailed balance. It is to be
     chosen and stated in the pre-registration.

The explanation for readers is in `docs/MODEL.md` §4.

## Pre-registration (2026-09-27, before any code or arm)

**Pool shapes** in Neo4j across the 92 pathways, as connected groups of forms
linked by curated forward/reverse pairs: **213 pools**, of which 178 have 2
forms, 35 are chains or branches of 3–9 forms, and none is a simple ring.

**Detection (generator, `pools.csv`):**
- A pool is found at **node (uuid) level**. For a curated pair (F: A → B,
  R: B → A in one pathway), the node-level loop A_u → F → B_u → R → A_u must
  exist: the same A node feeds F and receives R's output.
- Stid-level pairs whose occurrences the uuids separated are NOT pooled, so
  artefact merges are not re-created.
- Forms joined by such loops form one pool.
- The file lists, per pool, each form node (and whether it is the base form)
  and each transition (from form, to form, reaction node).

**Base form** (the least modified):
- the fewest modified residues summed over the form's leaves (Neo4j);
- then the fewest components;
- remaining ties go to the smaller stId, and are counted and reported.

**Baseline** (the one assumption, **φ₀ = 0.1**, Adam): detailed balance, with
each step away from the base form holding (φ₀ / (1 − φ₀)) = 1/9 of its
predecessor's share.
- 2 forms: 0.9 / 0.1.
- A 3-form chain: 0.890 / 0.099 / 0.011.
- The baseline rates follow from this; branched pools need no extra assumption.

**Solve (`DS_CYCLE_MODE=off` default | `balance`, with `DS_CYCLE_PHI=0.1`):**
- Each transition's drive is its reaction evaluated with the existing semantics
  and its source form held at baseline, divided by baseline.
- Rates are baseline rate × drive.
- π is solved from πQ = 0, Σπ = 1.
- The supply s is the OR-mean of the forms' producers outside the pool, else 1;
  a pinned form sets the pool.
- Form value = baseline × s × π / π₀. Each transition's reaction node carries
  its flux, relative to baseline.
- This is applied inside the component iteration every sweep, and is
  order-invariant.
- The catalyst ⊣ source-form depletion edge inside a pool is not applied (the
  balance already contains it).
- Autocatalytic pools take the root continuous from baseline.
- A pool the solve cannot handle falls back to iteration, and is counted.

**Arms:**
- **Build:** one build from the generator with `pools.csv`; its networks are
  otherwise identical to the canonical build.
- `ctrl039`: `DS_CYCLE_MODE=off`.
- `bal01`: `balance`, φ₀ = 0.1. **Primary.**
- `bal05`: `balance`, φ₀ = 0.5. A declared sensitivity arm, never used to choose.

**Adopt `bal01` only if all hold** (against `ctrl039`, same build):
- curator held-out net > +15, p < 0.05;
- experimental net ≥ +15, p < 0.05, excluding RAF (RAF reported separately);
- no pathway, tuning or held-out, loses more than 10;
- gains span at least 2 pathways and 5 perturbations.

**Predictions:**
- convergence improves (fewer non-converged solves);
- RAF's collapse stops (RAF reported);
- the knockout and overexpression of reversal enzymes (phosphatases, GAPs)
  gain cases.

## Amendment 1 (2026-09-27, after the adversarial review, before any arm)

The first build (`20260927-2257_c210f21_pools039`, 102 pools) was reviewed
before any arm ran. The review found that detection did not implement what this
pre-registration describes, plus defects in the solve. No arm was run on that
build, and it is not used.

**Detection, as found.** Only 24 of 102 pools were pure interconversions of one
protein.
- 61 were enzyme binding cycles (E + S → E:S → E + P).
- 17 had forms sharing no reference entity.
- 11 used proteasomal degradation as the "reverse" step.

The query asked only for "A → B and B → A". It never checked "forms of one
protein", which this pre-registration's own text requires.

**Detection, amended (generator `1358d03`):**
- Each transition reaction's only non-small-molecule input is the source form,
  and its only non-small-molecule output is the other form. Ubiquitin counts as
  a co-substrate. This excludes binding, release, degradation, and reactions
  taking two forms.
- The two forms share a reference entity of an EWAS leaf.
- A reaction node that would be a transition more than once is dropped, with
  its loops. Its node would be written by two fluxes, which gave a
  label-dependent value and no convergence (R-HSA-5654736, residual 0.79).
  They are counted.

**Base form, amended.**
- The first build decided a donor rule first, then components, then residues.
  The donor rule was not pre-registered, and the order was not the one
  pre-registered. 30 of 78 exchange pools were oriented backwards; all 24 pure
  pools were oriented correctly.
- Now: residues first (as pre-registered). Then the donor rule, added because
  RAS:GDP and RAS:GTP have no residue difference: the product of the direction
  that consumes a group donor (ATP, GTP, SAM, acetyl-CoA, NAD+, ubiquitin) is
  the modified form. Then components.
- Undecided pools fall back to the smaller stId, and are counted.

**Solve, amended:**
- **Supply** comes only from producers outside the pool's cyclic component.
  A producer inside it can be fed by the pool itself (11 pools).
- **Regulariser:** rates are k·(u + 1e-9), not k·u + 1e-12. With every drive at
  zero, the split now tends to π₀. Before, the constant chose it: a GEF+GAP
  double knockout read CDC42:GTP at 5.0x.
- **Depletion and inhibitor edges into a form from outside the pool** are
  applied as a fold on that form. Before, pool management silently dropped them
  (13 edges). The depleted share is removed, not redistributed.
  - The pool's own catalyst ⊣ source-form edges stay unapplied, as
    pre-registered.
- **Transition nodes** carry their flux as drive × the fold of the source form.
  This equals s·u·π_i/π₀_i when no modifier applies.
- **`DS_SCC_METHOD` `pool*` / `minimize`** never run the sweep that updates
  pools, so under them every pool is reported unmanaged. Before, they were
  miscounted as solved.

**Stated limitations (not fixed):**
- Influence scores do not see the pool rule.
- The autocatalysis rule is not implemented. No detected pool has a form
  catalysing its own pool's transition; this is checked on the build and
  reported.
- Several pools in one component are updated in a fixed order each sweep. The
  result is order-invariant only at a unique fixed point.

The arms, gates and predictions above are unchanged. They run on a rebuild
(`pools039b`) from generator `1358d03`. Its pool counts are recorded below
before any arm.

**Rebuild `20260927-2322_1358d03_pools039b` (recorded before any arm):**
- **25 pools**, all two-form (50 forms, 204 transition-reaction copies), in 12
  pathways. RHO GTPase cycle (R-HSA-9012999) has 12 of them; RAF/MAP (R-HSA-5673001) has 1.
- 0 orientation ties; 0 multi-use reactions dropped.
- At the stId level the amended query finds **55 form pairs in 19 pathways**.
  The node-level check keeps the 25 whose loop exists in the network.
- **Why so few:** the loose query found 203 distinct form pairs. 164 fail the
  one-protein-interconversion test, almost all because they are the enzyme
  step E ⇄ E:S (USP9X ⇄ Ub-SMAD4:USP9X, PPM1A ⇄ p-SMAD:PPM1A, …).
  - Reactome usually curates a modification cycle in several steps: bind,
    modify, release.
  - The substrate cycle SMAD4 → … → Ub-SMAD4 → … → SMAD4 runs through
    enzyme:substrate complexes, and a one-reaction A → B → A pool cannot see it.
  - These **multi-step modification cycles** are the larger class. This rule
    does not cover them, and treating them needs a design of its own: the
    enzyme is a second conserved species. It is not attempted here.
- **Consequence for the gates:** with 25 pools in 12 pathways, the
  held-out > +15 and experimental ≥ +15 gates may be out of reach. The gates are
  not lowered. A result under them is reported as it is.
