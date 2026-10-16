# specs/064: duplicate variant reaction copies (generator bug fix)

**Status:** measured 2026-10-15; the fix is LNG PR #114 (unmerged). **No gate
was registered before measuring**: this is a bug fix found by the week
review (docs/reviews/2026-10-15-week-review.md, finding 2), measured with
identical solver code on both catalogs, and the decision is Adam's.

## The defect

The generator's variant emitter enumerated an output complex's BOUND slot
when the complex also had a slot shared with the inputs, then overwrote it:
|members| identical copies of each real reaction copy, each its own node.
Canonical build `20261012-1401_01206ca`: 2,906 exact-duplicate reaction nodes
(same inputs, outputs, modifiers) across 36 pathways. Under OR-mean a
duplicated reaction weighs 2–6× in its outputs; 4 reactions were pooled only
because duplicates pushed them over `LNG_MAX_VARIANTS`.

## Measurement

Variant build `20261015-2251_5e49bc9_dedup` (generator `5e49bc9` = LNG main
`01206ca` + the fix): 12 exact duplicates remain (from another cause, not
traced). Solver `76d7306` (specs/061 amendment 5) on both catalogs:
`20261012-1401_01206ca/results/76d7306/rtwin58g6` and `20261015-2251_5e49bc9_dedup/results/76d7306/dedup_g6`.

| | net | fixed / broke | p |
|---|---|---|---|
| curator held-out | −4 | 0 / 4 | 0.125 |
| curator tuning | 0 | 0 / 0 | |
| experimental | 0 | 0 / 0 | |
| convergence | curator 1,654 (1,653), experimental 222 (221) | | |

86 cases are not paired (their pinned uuid count changed: RUNX1 76, NoRC 2,
plus 8); each is right in both or wrong in both, so they move nothing.

**The 4 losses** are one perturbation: SH3GL1 (endophilin) knockout in
Signaling by EGFR, readouts 180523, 179838, 109783 and 180500 (expected UP).
The canonical network duplicated GAB1 phosphorylation by EGFR (R-HSA-177930),
its dephosphorylation by SHP2 (R-HSA-177924) and CBL binding/ubiquitinating
Sprouty (R-HSA-183089), 14 extra copies each, doubling their weight in the
EGF:EGFR complexes they share. Canonical: UP at every cap (5.2×–97×). Without
the duplicates the loop collapses instead (DOWN at every cap, value
unsettled, 0.0–0.22×). The canonical right answer rested on the duplicate
weighting.

## Recommendation

Merge LNG #114 and make the deduplicated build canonical: the duplicates are
a generator bug, the measured cost is one perturbation in one loop (not
significant), nothing else moves, and convergence improves slightly.
