# specs/042 problem: a removal has no route to what it removes

For two blind derivations, as for specs/039 and 040.

## Context

DeltaSignal propagates fold-changes through a Reactome-derived logic network
(see docs/MODEL.md):
- a reaction is the AND of its inputs, divided by its inhibitors;
- an entity is the OR-mean of its producers.

**Inputs only activate reactions.** A reaction never lowers what it consumes.
Lowering happens only through explicit negative edges:
- curated negative regulators;
- generator-derived `depletion` edges (catalyst → substrate, for phosphatases
  with Pi as an output and ubiquitin ligases with Ub as an input; specs/011
  bounds them).

Curated interconversion cycles (A ⇄ B) are now conserved pools (specs/039).

## The defect (specs/041 step 3, motif M1; specs/034 §12d)

When a reaction **removes** an entity and does not give it back, nothing in
the network lets the remover lower it:
- **WNT:** the free β-catenin root has no depletion in-edge. Every other CTNNB1
  copy is a *product* of destruction-complex reactions assembled from APC and
  AMER1. So APC KO reads 0 and APC 80x reads 100, **both signs inverted**.
  That is about 8 experimental and about 70 curator errors.
- **HDR:** branch competition. HR and SSA, SDSA and resolution compete for one
  intermediate; raising one branch should lower the other.
- **p27 (CDKN1B) and RB1:** a binder or repressor whose only edge is a
  positive input to its sequestering reaction.
- **MP-BioPath's curators added these by hand** as "brakes": about 470 cases on
  their networks (binders and traps, DUBs, GAPs, repressors, endocytosis;
  specs/034 §12d). Every case a brake fixes reads NORMAL without it.

## Prior results that must not be repeated

- **"General consumption" (blanket).** "More of one input consumes the
  co-input", as a divide-form depletion on every reaction: NEGATIVE
  (2026-07-16).
- **Naive substrate depletion,** catalyst → input depletion on every catalysed
  reaction: −14pp (2026-05-29).
- **specs/019:** three sink-bridge interventions, all null or negative.
- **specs/035:** holding conserved loop inputs at baseline lost experimental
  −122.
- **Depletion edges are load-bearing.** Removing them is held-out −90.

## Questions

1. State a rule, from mass balance (turnover at steady state), for when a
   reaction *removes* an entity, and how the remover then acts on it.
   - Candidate: if X is supplied at rate s and removed by non-returning exits
     with drives u_i, then X = s / Σ w_i·u_i, relative to baseline. A remover
     acts on X by division.
   - Which exits count as "non-returning"? Degradation, a dead-end complex,
     modification into a form nothing consumes. What about consumption whose
     product goes on to signal?
   - How are the w_i set at baseline without a fitted parameter?
2. How does it differ from the failed blanket rules, precisely, and why should
   it not repeat their losses?
3. How does it compose with specs/039 pools (interconversion), specs/011's
   depletion bound, and specs/022/040's self-contained inhibitors?
4. Can the generator detect it from Reactome, or must the solver? State the
   detection, and how to census it before any arm.
5. Work the cases numerically: WNT APC KO / 80x on nuclear β-catenin; HDR
   (raising HR on SSA); p27 KO on CDK2 activity; RB1 KO on E2F targets.
   Expected from biology: APC KO → β-catenin UP, APC 80x → DOWN; RB1 KO →
   E2F targets UP.
6. **What can go wrong.** Where would the rule wrongly lower something? What
   does it predict outside the four motivating pathways? Name the census and
   the checks.

Constraints: no parameter fitted to the evaluation; baseline exact; switchable
and byte-identical when off. Judged on the pathways it was not derived from,
on both axes, and against the specs/041 expected-biology review.
