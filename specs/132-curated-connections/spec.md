# specs/132: connecting pathways through Reactome's curated connections

**Status:** registered 2026-11-03, before any code or count.

## The principle (Adam, 2026-11-03)

**Pathways are connected only through the connections Reactome curates. Connecting them any other way would never
work.** Two pathways that contain the same molecule are not thereby connected; only a curator's statement says that
one pathway hands something to another, and through which reaction. This holds for tests and for regular use.

Every earlier combination joined pathways by entity identity: May 2026 (stable ids merged), July 2026 (all 92 pathways
merged on identical nodes), specs/056 (a copy produced in any pathway bridged to a root copy of the same entity). None
used a curated connection. Their losses (shared hubs flooded; free EGF read at the 10× de-repression ceiling in
Signaling by EGFR and carried into PIP3; free KDM1A drained to 0 and carried in) are what identity matching produces.
They are not evidence about combining pathways as Reactome connects them, which has never been measured.

## The curated connections (Neo4j, Release 97, read 2026-11-03)

- **Pathway level, `hasEncapsulatedEvent`**: the process nodes a curator draws in a pathway diagram pointing to another
  pathway (6,614 links). Examples: Mitotic G1 phase and G1/S transition → 6 pathways (PIP3 activates AKT signaling,
  SMAD2/3 transcription, Signaling by PTK6, three senescence pathways); PIP3 activates AKT signaling → 19 (EGFR, ERBB2,
  MET, MTOR, TP53 regulation, ...).
- **Reaction level, cross-pathway `precedingEvent`**: a reaction whose curated preceding reaction sits in another
  pathway (about 4,600 human pairs, counted as pairs whose direct parent pathways differ). This names the reaction that
  hands off, so it names what is handed off.

## Rules (each a separate catalog build from the frozen build, `bench/context_build.py` extended; default off)

For a pathway P being solved:

- **L1 (reaction level).** For each pair (r in P, s in pathway Q, s precedes r): a bridge `pos,or,context_bridge` from
  Q's copy of each entity that is an **output of s** to P's copy of the same entity that is an **input, catalyst or
  regulator of r** and a root in P. Q contributes only the ancestors of those bridge sources (as specs/056's pruning).
  No other entity is bridged.
- **L2 (pathway level).** Q must also be one of P's `hasEncapsulatedEvent` pathways (L1 restricted to pathways the
  curator drew in P's diagram).
- **Fallback (F).** Applied with L1 and with L2: P's own answer stands for any (perturbation, readout) that P alone
  routes; the connected build is read only where P has no route from the perturbed gene's pins to the readout.

Fixed now, not swept: cofactors, drugs and pathogens are never bridged (as specs/056); bridges run only into P (feed-in,
no cycle through P); neighbours of a neighbour are not added; P's readouts are scored on P's own nodes. Neighbour
pathways not in the frozen 92-pathway build are built with the frozen generator (tag `freeze-2026-10-26`); the solver
is unchanged.

## Order of work, and what is read when

1. **Ceiling count (no solving).** On the FIT datasets from cells (experimental, korkut, decryptM, HPN-DREAM): of the
   frozen model's wrong cases, how many have no route in P but a route through L1, through L2, and through identity
   matching (specs/056's bridges, for comparison). If L1 reaches fewer than 2% of wrong cases on every dataset, stop:
   curated connections cannot be what the errors lack.
2. **Census of the bridges L1 makes**: how many, which entities, and how many are free ligands or released forms
   (the specs/056 loss class), before any arm runs.
3. **Arms** with `scripts/run_arm.sh` on each dataset: control (frozen), L1, L2, L1+F, L2+F. The identity-matched
   specs/056 `ctxin` is re-run on the frozen build as a comparison only; it is never a candidate.

## Gates

- **Candidate gate (FIT, measured):** balanced accuracy (3-class) over experimental (RAF set aside) and korkut pooled,
  arm minus control, bootstrap by perturbation (10,000, seed 132), interval above 0; and neither dataset's point
  estimate below 0. decryptM and HPN-DREAM reported, not gated.
- **Curator guard (SELECT):** curator answers are scoped to one pathway by construction (a curator predicted from
  inside P), so the curator set cannot judge a connection. It is used only as a guard: under F, a case P routes alone
  must not change (checked, it is the rule's definition); without F, curator held-out losses are reported and every
  broken case in a sample of 20 is classed as connection-explained (the curated connection carries the change) or not.
- At most one candidate (the best by the candidate gate's point estimate; ties to the narrower rule, L2 over L1, F
  over no F) goes forward to specs/125 stage 5: scored once on Lun 2019 and the CPPA test portions beside the frozen
  model, under specs/127's protocol.

## Regular use

The pathway alone stays the default view, since that is what curators reviewed. A rule that passes becomes an option
("include connected pathways"), and every readout whose answer came through a bridge names the curated connection it
used.

## Predictions

- **P1.** L1 reaches more wrong cases than specs/056's identity bridges per bridge added (it adds far fewer bridges).
- **P2.** L1's bridges carry few free ligands or released forms, so the specs/056 loss class mostly does not recur.
- **P3.** Gains concentrate in the cell-cycle pathways (Mitotic G1, S phase, checkpoints), as every earlier
  combination's did.
- **P4.** F matters more on korkut than on experimental (korkut readouts sit in many pathways at once).

Nothing here changes a default or reads a TEST set; Lun 2019 stays sealed.
