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
