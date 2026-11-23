# Feature Specification: Perturbation at the gene's own copies, wherever they sit (option B)

**Feature Branch**: `specs/171-stock-copy-perturbation`

**Created**: 2026-11-23

**Status**: REGISTERED 2026-11-23 (design and gate). Nothing is coded or run yet. Research resolves the open items below
before the first arm, and any change to the gate is recorded as an amendment before measuring.

**Input**: Adam, 2026-11-23 (confirmed in conversation): "Then register full option B as its own
spec." His rule:
- Roots are not real. They mark where a curator cut the pathway.
- Roots and terminals are computed on whatever network we have, curated plus literature.
- Perturbations and drugs act on nodes that are not roots.

## The change

**Today (specs/023, 024, 058, the protocol of record):** a gene perturbation pins the nodes that are, or contain, the
gene and have no counted in-edge. If there are none, it pins the specs/024 own-cycle form. A gene whose own copies all
have a producer is not perturbed at all, and its cases are invalid.

**Under B:** a perturbation clamps the gene's **stock copies** whatever feeds them, and the clamp cuts their in-edges
for that solve.
- **Stock copies** are the copies the abundance readout already reads (specs/134 amendment 3): the gene's starting and
  translated copies, and sets or pools of them. They exclude any copy that the pathway makes from the protein (bound,
  modified or released forms).
- **Fallback:** a gene with no stock copy (for example, one that occurs only inside complexes) keeps today's rule.

So "the gene" means one set everywhere: what a perturbation clamps, and what a total-protein readout reads.

**Roots** then mean only where the network's inputs enter, which matters for stimuli and the baseline, and they are computed on the network as served.

Four structural rules in the solver also identify a gene's own copy as "no in-edge":
- recycled twins (specs/059)
- closed carrier twins (specs/049)
- component twins (specs/084)
- RNA twins (specs/100)

Under B they ask the same identity question, so their "starting copy" becomes the stock copy.

## Arms (on the served catalog, curated networks only; no literature)

| arm | bench | solver |
|---|---|---|
| R (control) | today's protocol | defaults |
| BP | stock-copy pins | defaults |
| BF | stock-copy pins | the four rules read stock copies (new flag, default off) |

All arms are run with `scripts/run_arm.sh` at one commit, on curator, experimental and korkut.

## What is measured, and the gate

B changes which cases are scored, so every comparison is paired on the cases scored in both arms. The cases only B
scores (genes newly perturbable) are reported separately, with their accuracy and change-call precision.

**Gate (registered):**
- on the primary pool (korkut with own markers excluded, plus MP-BioPath experimental, paired cases), net(BF vs R) ≥ 0;
- the paired sign test of fixed against broken does not show harm (two-sided p ≥ 0.05, or net > 0);
- a sensitivity guard: experimental change-call recall in BF is not more than 0.02 below R (the false-alarm rule of
  2026-11-16).

Curator held-out is reported, never gated. It is a SELECT set chosen under today's protocol, so a protocol change is
judged against data that was never used to choose rules (2026-11-20 rule, applied by analogy).

Every changed case on the primary pool is listed. At least the five largest groups (perturbation × pathway) are traced
to the pin or rule that moved them.

**Prediction:** I expect few changes on the curated networks, because their stock copies are mostly roots already.
Changes come where a gene's free form is also translated in the pathway (TP53 targets such as MDM2), so a knockout
also zeroes the translated copy.

**What PASS would mean:** B becomes a candidate protocol of record, for Adam to decide. Adoption would re-baseline every
headline, keeping the old numbers labelled as the root protocol (as specs/023 did). Literature bundles would then be
re-run under B (specs/168 and 170 arms).

## Open items for research (resolved before any arm, recorded in research.md)

1. The exact stock-copy set, computed in both the bench and the engine (`abundance_copies` in `src/api/readouts.jl`,
   the bench's `ABUNDANCE_RULES`). In particular, whether an mRNA copy is counted as the protein. If it were, a
   translated copy would read as "made from" the protein.
2. Clamp values for copies that are not roots. A knockdown pins at `DS_PERTURB_UI_DOWN` and an over-expression at
   min(80, 0.8 × ceiling), as today. Is clamping a translated copy and its starting copy to the same value right, or
   should the dose split across copies? Decided from the record before the first arm.
3. Engine parity. `/api/pins` (specs/140) must serve the same pins as the bench under B. The parity gate is re-run on
   the frozen catalog.
4. The cycle-balance pools: a clamped copy inside a pool. Confirm from the code that a pin there is handled as a pinned
   root in a pool is now.

## Requirements

- **FR-001:** the rule, tests that run real code on fixtures (a translated copy is clamped; a bound form is not; the
  fallback fires only with no stock copy), and the scorer are committed before any arm is run.
- **FR-002:** the BF solver flag defaults off. Under `R` the solve is bit-identical to today's, checked by replay.
- **FR-003:** results are recorded in `research.md` whatever they say. Adoption is Adam's decision, never a default
  flipped by this spec.

## Non-goals

- No literature bundles in this spec's gate. They are re-run under B only after the curated measurement.
- No stimulus or baseline input. That is a separate spec once the context format exists.
- No change to served pins or defaults.
