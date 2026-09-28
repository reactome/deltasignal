# Problem: a propagation rule for interconversion cycles in a signed logic network

## The model as it is

DeltaSignal propagates perturbations through a signed logic network built from
Reactome.

- **Nodes** are physical entities and reactions. Each node carries an activity,
  expressed here as a **fold-change relative to baseline** (baseline = 1; 0 =
  absent; capped at 100).
- **Perturbations are pins.** An overexpressed gene's ROOT node is fixed at 80;
  a knockout is fixed at 0. Everything else is solved.
- **Reaction node** value =
  `(product of its positive AND inputs' folds) × (per inhibitor: 1 / inhibitor fold)`,
  clamped to [0, 100].
  - Positive AND inputs are its substrates, catalysts and positive regulators.
  - This is "hill_sat" AND with "divide" inhibition.
  - Several producers of one entity are combined by their mean (OR).
- **Entity node** value = the mean of the reactions producing it (a single
  producer passes straight through).
- **Solving:** the network is condensed into strongly connected components,
  processed in topological order. Inside a cyclic component the solver iterates
  x ← F(x) with damping 0.5 until the residual is below 1e-6 or 500 sweeps. The
  result must not depend on node labels or sweep order.
- **Readouts** are compared with curator and experimental expectations of
  UP / NORMAL / DOWN, using fold cutoffs of about 0.85 and 1.15.

## The structure in question

Reactome curates **interconversion cycles**, pairs of reactions within one
pathway:

    F: A (+ co-inputs) --[catalyst E_f, regulators]--> B (+ co-outputs)
    R: B (+ co-inputs) --[catalyst E_b, regulators]--> A (+ co-outputs)

Examples:
- phosphorylation / dephosphorylation (kinase / phosphatase);
- GEF-driven GDP→GTP exchange / GAP-driven hydrolysis;
- ubiquitination / deubiquitination;
- methylation / demethylation;
- binding / dissociation.

**In the catalog:** 259 such pairs in 55 pathways.
- 93 have a catalyst on both directions, 150 on one, 16 on neither.
- Multi-state cycles (A → B → C → A) also occur.

In the network, A → F → B → R → A is a directed cycle. Because every edge
multiplies folds, a cycle has gain exactly 1 at baseline (a knife-edge):
- a small push rails the cycle to 100, or drains it to the absorbing all-zero
  state;
- which one happens can depend on sweep order;
- one pathway (RAF/MAP kinase) collapses to 0 under every upstream
  perturbation.

**Hand-curated equivalent.** The earlier tool, MP-BioPath, had loops cut by
hand. Its curators also added an inhibition "reverse enzyme ⊣ modified form" in
about 219 places (phosphatase ⊣ phospho-protein, GAP ⊣ RAS:GTP, …). Only about
24% of those correspond to a curated reverse reaction catalysed by that enzyme.

## What the ground truth expects (curator reasoning and cell experiments)

1. **Forward enzyme** E_f knocked out → modified form B DOWN; overexpressed → B UP.
2. **Backward enzyme** E_b knocked out → B UP; overexpressed → B DOWN.
3. **The cycled protein's gene** (A and B share it) overexpressed → **both A and B
   UP**; knocked out → both DOWN.
4. An **upstream activator of E_f** propagates through B to downstream readouts
   with its sign.
5. **Single-member effects count.** When E_f is a set of kinases, knocking out
   one member is expected to change the readout; a rule treating members as
   interchangeable measured badly.

## What is asked

Derive a rule for computing A and B (and the reaction nodes F and R) inside
such a cycle. It must be:

- **(a)** consistent with expectations 1–5;
- **(b)** exactly baseline when nothing is perturbed;
- **(c)** free of the gain-1 knife-edge and the all-zero trap;
- **(d)** local, computable from the cycle's own inputs (E_f, E_b, co-inputs,
  regulators, the supply of the cycled protein) within the fold-change
  framework;
- **(e)** defined for a one-direction-catalysed pair, for a multi-state cycle,
  and for a cycle nested inside a larger cyclic component;
- **(f)** free of parameters fitted to the evaluation data. Any constant must be
  a stated assumption with a biological justification, such as the baseline
  fraction in the modified form.

**Please give:**
1. the rule in closed form, with its derivation (for example, from a
   steady-state rate balance);
2. how it composes with the existing AND/OR/divide semantics at the cycle's
   boundary;
3. the assumptions, and the sensitivity of the predicted UP/NORMAL/DOWN calls to
   each assumption;
4. the failure modes: cases where the rule would predict the wrong sign;
5. a small set of test fixtures with expected fold values.

Be concrete, and challenge the framing if it is wrong.
