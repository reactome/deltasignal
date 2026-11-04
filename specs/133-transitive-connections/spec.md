# specs/133: the pathway alone where it can, all of Reactome's curated connections where it cannot

**Status:** registered 2026-11-03, before any code or count for it. Follows specs/132 R3a. Amendments 1–2
(2026-11-04) fixed two construction defects found by tracing. Measured (research R5): A2+F and T2+F pass the gate on
balanced accuracy (+0.015 [+0.004, +0.027]) while korkut raw accuracy falls; every changed case traced, none wrong in
construction; depth past one hop adds only PI Metabolism's rails. Not adopted; nothing changes a default.

## The decision (Adam, 2026-11-03)

"Use the pathway's answer when it can, and if it can't, use the combined pathways from all of Reactome. Then look at
each one to see why it helped or broke things, to figure out if it is being done right or not."

Pathways are still connected only through Reactome's curated connections (specs/132's principle): a cross-pathway
`precedingEvent` whose preceding reaction outputs the entity the following reaction takes as input, catalyst or
regulator. Never by shared molecules.

## What specs/132 showed had to change

specs/132 R3a traced every case its rule changed. Three things were wrong in how the combination was made, not in
the connections:

1. **Cell-to-cell signals.** IL4 (and IL21) handed from *Affinity selection of immunoglobulins* are cytokines a T
   follicular helper cell secretes onto a B cell. A one-cell model cannot take one cell's output as another's input.
   Reactome's cell-type annotation labels no event (specs/132 R5); the compartment does: of 48 bridged entities, these
   are the only two in the extracellular region.
2. **Bridges moved roots in the solver.** The bench already reads roots without context bridges, but the solver's
   structural rules (specs/084 component twins, specs/059 recycled twins, specs/100 RNA twins, carrier loops) read the
   bridge as an in-edge, so a starting copy stopped being one (ATM in *Cell Cycle Checkpoints*, 4 experimental cases;
   21 korkut AKT readouts with no bridge on their route).
3. **One hop only.** specs/132 added a pathway's direct sources and never their sources. "All of Reactome" means
   following connections back as far as they go.

## The rule (one candidate, T)

For a scored pathway P:

- **Connections followed transitively, as a tree.** Depth 0 is P. A pathway Q joins at depth d+1 when a curated
  handoff runs from a reaction in Q into a root copy, in a pathway X at depth d, that X keeps (for X = P every root;
  otherwise an ancestor, inside X, of X's own bridge sources). Each pathway joins once, at its least depth; bridges run
  only from depth d+1 into depth d, so the context has no cycle and never feeds back into P. No depth limit.
  Handoffs are specs/132's L1 (any `precedingEvent` across pathways, not limited to encapsulated pathways). A source
  pathway is the lowest pathway with a diagram holding the preceding reaction (as specs/132); every such pathway
  reached is built with the frozen generator (tag `freeze-2026-10-26`).
- **Never bridged:** cofactors, drugs, pathogens (as specs/132), and **any entity in the extracellular region**
  (GO:0005576; a set is excluded when any member is). This also drops an autocrine ligand; accepted, since no
  annotation separates the two.
- **Solver, `DS_CONTEXT_BRIDGE_ORDER`**: `first` (default; every rule sees bridges, as specs/056 and 132 ran) | `last`
  (context bridges are set aside before every structural rule, as recycle edges are (specs/081), and added back after
  them as `pos,or` in-edges). A bundle without bridges solves identically under both. Reports
  `context_bridge_order` and `context_bridges`.
- **Fallback F, which is the decision itself:** P's own answer wherever P alone routes from the perturbation's pins to
  the readout; the combined answer elsewhere.

## Arms (scripts/run_arm.sh, all four axes)

- **c0**: each pathway alone, from the new base build (must reproduce the frozen model case for case).
- **A**: one hop (specs/132 L1) with the extracellular exclusion, under `DS_CONTEXT_BRIDGE_ORDER=last`. Attribution
  only: what the two fixes alone do.
- **T**: the rule above, under `last`. The candidate; read with F.

## Order of work

1. Census (no solving): the closure: pathways reached, depth distribution, bridges per depth, entities excluded as
   extracellular, and how many pathways must be newly built. The ceiling (control-wrong cases with no route alone but
   a route through T) on experimental, korkut, decryptM, curator.
2. Build, run the three arms, score with `combine_132.py`.
3. **Trace every case T+F changes** with `trace_132.py`, grouped by the chain of connections it crossed, and classify
   each group: the combination done right (a curated same-cell chain carrying the signal), wrong in construction
   (fix it, re-register), or right in construction but the receiving answer wrong (the error belongs to a single
   pathway, as PIP3's AKT readouts in specs/132).

## Gate

As specs/132: pooled 3-class balanced accuracy over experimental (RAF set aside) and korkut, T+F minus c0, bootstrap by
perturbation (10,000, seed 132), interval above 0 and neither dataset's point estimate below 0. Accuracy and macro-F1
reported beside it. decryptM and curator held-out reported, not gated. korkut also reported once per (drug set,
readout) changed, so one answer copied into several pathways is visible.

**Honest limit.** specs/132's traces of experimental and korkut cases chose fixes 1 and 2, so this gate is no longer
blind on those datasets. It can reject; a pass is not confirmation. Confirmation would be one scoring on data not
used here (Lun 2019, or the CPPA test portions), and opening either is Adam's decision.

## Predictions

- **Q1.** ATM's 4 experimental cases in *Cell Cycle Checkpoints* return to the control's (correct) answer under A.
- **Q2.** The IL-4/13 curator losses (66 under specs/132 L1+F) are gone under A and T.
- **Q3.** Transitivity adds pathways mostly at depth 2–3 and reaches few new wrong cases: the ceiling rises by less
  than half over specs/132 L1 on experimental.
- **Q4.** More than half of T+F's korkut changes are still AKT readouts that repeat *PIP3 activates AKT signaling*'s
  own answer; that error is not addressed here.
- **Q5.** The gate fails (the experimental gain is a handful of PTEN cases; korkut loses on raw accuracy).

Nothing here changes a default; `DS_CONTEXT_BRIDGE_ORDER` defaults to today's behaviour. Lun 2019 stays sealed.

## Amendment 1 (2026-11-04, before any score of A or T was read)

**Found while checking arm A** (research R3): under `last`, the bridged starting copy of p-AKT1/2/3 in *Transcriptional
Regulation by TP53* is a starting copy again, and specs/059's recycled-twin rule, which joins copies of one entity by
identity, gave PIP3's imported copy of p-AKT an in-edge from TP53's copy. With the bridge restored that is a loop
PIP3 → TP53 → PIP3: the context fed back, which the rule forbids, and the two pathways were joined by a shared molecule,
which specs/132's principle forbids. It moved 5 PTEN cases (4 to wrong). Under `first` the same kind of edge can form
from any unbridged starting copy, so specs/056 and specs/132 arms may hold some too (counted in R3).

**Change.** A context bundle records each node's pathway (`context_origin.csv`, uuid → pathway, written by
`context_build.py`). Under `DS_CONTEXT_BRIDGE_ORDER=last` the solver drops every edge a structural rule added between
nodes of two different pathways, so each pathway's rules act on its own copies, as when it is built alone; pathways
meet only at curated bridges. Reported as `context_cross_edges` (edges dropped). `first` is unchanged. Checked by
running the rules on each pathway's part separately and comparing with the filtered result, bundle by bundle.

A and T are re-run with this; the first runs are kept and reported (R3) as the record of the defect.

## Amendment 2 (2026-11-04, after the corrected runs A1/T1 were scored, before any rule choice)

**Found tracing A1/T1** (research R4): every PTEN overexpression case in *Transcriptional Regulation by TP53* still read
wrong, now DOWN in both directions. The bridge sources, *PIP3 activates AKT signaling*'s p-AKT1/2/3, sat at the 100×
ceiling under both PTEN knockout and overexpression, while PIP3 solved alone reads PTEN overexpression DOWN (0.17×) and
the knockout UP. The PIP3 inside the bundle was not PIP3: the builder (specs/056's design, kept by specs/132 and this
spec) removes from a source pathway every reaction of the pathways it feeds, and TP53 transcription holds PTEN's
transcription and translation, which are what make PIP3 alone read PTEN right (specs/095, 100). Pruning a source to the
ancestors of its bridge sources is also not exact once the solver's twin rules add in-edges from non-ancestors.

That breaks the decision itself: "use the pathway's answer" must mean the source pathway's own answer too.

**Change.** `context_build.py --whole`: a source pathway enters the bundle whole, not reduced and not pruned, so it solves
exactly as alone except for the bridges it receives itself (deeper tree levels under T). Copies of a shared reaction in
two pathways are distinct nodes and meet only at curated bridges, as amendment 1 requires. The tree is unchanged (which
roots may receive deeper handoffs is still decided by ancestry, now over the whole pathway). Checked: for the one-hop
bundle, every source pathway's nodes take the same values as that pathway solved alone, under the same pins.

Arms `A2`, `T2` (bundles `.A2`, `.T2`); A1/T1 are kept as the record. Gate, predictions and F unchanged.
