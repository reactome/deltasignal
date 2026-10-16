# specs/061: an entity-own depleter may de-repress only while its source copy holds

**Status:** registered 2026-10-14, before any code or arm. R1 (`rtwin58g`) and R2 (`ownsrc`) fail: the source test was circular when the target is its own source. **Amendment 1 (`rtwin58g2`, R3) passes the registered gate, and amendment 2 (`rtwin58g3`, R4: two review fixes) reproduces it exactly:** held-out +32 (33 / 1), tuning +112, experimental +7 (7 / 0), worst pathway −1. **Adopted as the default 2026-10-15** (Decision, end of file); post hoc, korkut2015 check inconclusive (specs/062), Release 98 is the replication. Post hoc: the
third refinement in the chain specs/058 → 059 → 060, designed after seeing
which cases each arm moved. Adam approved continuing (2026-10-14).

## Why

specs/060 R2 found that SMAD2:SMURF2 -| SMAD2 and MDM2:TP53 -| TP53 are the
same motif with opposite answers:
- SMAD2 knockout: the complex falls because the TARGET's own entity fell, so
  SMAD2 must not rise (`rtwin58s`: TGF-beta −12, SMAD2/3 knockouts read UP);
- MDM2 knockout: the complex falls because the OTHER partner fell, so TP53
  should rise (`rtwin58e`: tuning −33, experimental −17 once the entity test
  forbade it).

No structural test separates them. Which partner fell is known only during
the solve.

## The rule (solver, `DS_DEPLETION_OWN_PRODUCT=suppress_only_source`, not a default)

As `suppress_only` (specs/016/060) for node-copy own products: such a
depleter may suppress its target, never de-repress it. In addition, for each
depletion edge D -| T whose depleter's entity CONTAINS T's entity (the
specs/060 amendment-1 pairs, by the bundle's containment table) and that the
node-copy test does not already flag:

- the **source copies** of the pair are the nodes with T's entity (same
  stable id) found walking backward from D along non-depletion edges,
  stopping at each such copy (the copies D is built from; T itself may be
  one);
- during the solve, the edge's factor is capped at 1 (suppress only) while
  ANY source copy reads below its baseline (by more than a relative 1e-6),
  and is the normal divide-form factor otherwise;
- a pair with no source copy is not capped (as `full`).

Structural part computed once per network; the test reads the current
iterate, like every other input. The solve reports `own_source_depletions`
(pairs gated).

## Census (before any code)

`bench/analysis/061/own_source_census.jl` on the canonical build
(`bench/analysis/061/census_01206ca.tsv`): **315 pairs in 31 pathways**; 290
have one source copy, 11 several, 14 none; in 202 the target is its own
source. Largest: TGF-beta (R-HSA-170834) 48, R-HSA-195721 43, R-HSA-73894
37, EGFR 34. The two motifs:
- SMAD2:SMURF2 (98920560) -| each SMAD2 [cytosol] copy: source c854b46a, a
  root copy; a SMAD2 knockout pins it to 0, so the rule caps;
- MDM2:TP53 pool (80e7a4cf) -| TP53 (8c4aa191): source is TP53 itself, which
  rises under an MDM2 knockout, so the rule does not cap.

## Arms

| arm | build | solver / bench | against |
|---|---|---|---|
| control | canonical `20261012-1401_01206ca` | defaults (`results/64d23cb`) | — |
| `rtwin58g` (primary) | canonical | `--server DS_RECYCLED_TWIN=on --server DS_DEPLETION_OWN_PRODUCT=suppress_only_source --bench DS_PIN_ROOT_IGNORE=depletion` | control; also `rtwin58s`, `rtwin58e` |
| `ownsrc` (record) | canonical | `--server DS_DEPLETION_OWN_PRODUCT=suppress_only_source` | control; also `ownsup` |

`--pair-on-valid`, RAF set aside, both readings reported.

## Gate (as specs/055–060)

Against the control: curator held-out net > 0 at p < 0.01; experimental net
loss (RAF set aside) not significant; no pathway losing more than 10 on
either axis; every experimental loss outside RAF traced. Also reported:
convergence counts, since the rule switches on a threshold.

**Because this is post hoc**, passing the gate is not enough to recommend a
default: the held-out and experimental sets have both now steered three
rules. If it passes, I recommend an independent check before adoption
(Adam's call): the korkut2015 TEST set (specs/051, values unread), scored
once, or a fresh Reactome release.

## Predictions (before any code)

- `rtwin58g` against `rtwin58s`: TGF-beta +12 (the SMAD2/3 knockouts), TP53
  unchanged (MDM2 knockout keeps de-repression), so against the control
  about held-out +33, tuning about +112, experimental about +7.
- Against `rtwin58e`: tuning and experimental recover the MDM2 knockout
  (+48 curator, +12 experimental).
- Unknown: the other ~250 pairs. Where the target is its own source (202),
  the rule only forbids a target from being de-repressed by its own fall,
  which I expect to move little. RNF168 overexpression (DSB) and the CBL
  class lost in `rtwin58e` may or may not return.
- `ownsrc` alone: close to `ownsup` (held-out 0, tuning +46, experimental
  +4); SMAD2/3 does not move without the twin edge and pins.
- Risk: the threshold makes the factor discontinuous at baseline;
  convergence may drop in pathways with these pairs inside loops.

## Amendment 1 (2026-10-14, after R1, before any amended code): the target is not its own source

R1 traced the MDM2 loss to a circular test: where the target is its own
source copy (202 of 315 pairs), its below-baseline value both results from and
decides the cap, so a dip locks in (TP53 8c4aa191 under the MDM2 knockout).

**Rule (`DS_DEPLETION_OWN_PRODUCT=suppress_only_source`, changed in place;
the R1 form is recorded at `cb53f3f`):** the source copies exclude the target
itself. The walk does not stop at the target: it passes through it (so a
copy reached only through the target still counts). Pairs left with no source
are not capped. By the census, 191 pairs lose their only source and 124
remain gated, including the three SMAD2 recycled copies (source c854b46a).

**Arm `rtwin58g2`:** as `rtwin58g`, at the amended commit. Same control, gate
and pairing; also reported against `rtwin58g` and `rtwin58s`.

**Predictions:** MDM2 KO returns to canonical (+42 curator, +9 experimental
against `rtwin58g`); TGF-beta stays 0 against the control; AKT1/AKT2 move
toward `rtwin58s` (+38 each). If so, against the control: held-out about +24
to +33, tuning about +100, experimental about +7, which would pass the gate
unless RNF168 (DSB −16, a coin) still loses more than 10 in one pathway — it
did in R1, so **the gate most likely fails on DSB**. If only that fails, I
report it as the specs/057 coin and do not change the gate.

## Amendment 2 (2026-10-15, after an adversarial code review, before any amended code)

The review of the code behind `rtwin58g2` found two places where it does not
do what this spec and specs/059 say. Both are fixed in place; the measured
`rtwin58g2` result stays recorded at `9a0b839`.

1. **A carrier free form is not a root copy** (specs/059). `apply_recycled_twins`
   runs after `apply_carrier_loops`. Under `closed_twin`, a closed free form
   with no root twin has lost every in-edge, so it passed the root test and
   could feed a recycled copy with a node held at baseline. Fix: a form listed
   in `carrier_loops.csv` is never a root copy, as it is never a recycled one.
2. **Source copies are what the depleter is BUILT FROM**: the walk back from
   the depleter follows only material edges (input, output, and the derived
   supply edges), not catalyst or regulator edges. Before, a copy of the target
   that only inhibited or catalysed the depleter's reaction counted as a source.

Also: `own_source_depletions` now counts only the pairs actually gated (pairs
the node-copy test already caps are left out), and the census is re-run on the
network as the solver sees it, after the carrier and twin edges.

**Arm `rtwin58g3`:** as `rtwin58g2`, at the amended commit. Same control, gate
and pairing; also reported against `rtwin58g2`.

**Prediction:** few predictions move against `rtwin58g2` (the defects need a
carrier free form with no root twin next to a recycled copy, or a copy of the
target that is only a modifier upstream of a depleter). If more than a handful
move, each is traced before anything else is concluded. The gate decides as
before; korkut2015 (now FIT) is reported for information only.

## Decision (2026-10-15): adopted as the default

Adam: "The decision is to include something that passes more tests. so I
would say yes." Adopted on the gate result (R3, reproduced by R4 after the
review fixes), knowing it is post hoc and that the korkut2015 check (specs/062)
was inconclusive. Three defaults change together:

- solver `DS_DEPLETION_OWN_PRODUCT=suppress_only_source` (was `full`);
- solver `DS_RECYCLED_TWIN=on` (was `off`, specs/059);
- bench `DS_PIN_ROOT_IGNORE=depletion` (was unset, specs/058); `none` restores
  the old pins.

The canonical numbers are re-benchmarked at the adopting commit and must
reproduce `rtwin58g3` case by case. Reactome Release 98 is the replication;
if the rule loses there, it is reverted.

## Amendment 3 (2026-10-15, after R5, before any code): decide the gate once, not on every sweep

R5: the gate is a switch read on the live iterate; inside a cyclic component
it toggles between sweeps, the solve never converges (TP53, TGF-beta), and
the scored value is the iteration cap's choice.

**Rule (`DS_DEPLETION_OWN_PRODUCT=suppress_only_source`, changed in place; the
amendment-2 form is recorded at `9e93c11`):** two passes. First the network
is solved with the entity-own pairs not gated (as `suppress_only`: node-copy
own products capped, entity-own pairs free). A pair is GATED for the second
pass if any of its source copies reads below baseline × (1 − 1e-6) in that
first solution. The second pass solves with the gated pairs capped at 1 and
the others free, fixed for the whole iteration. Only solves with at least one
entity-own pair run twice. The solve reports `own_source_depletions` as the
pairs gated in the second pass.

**Why this form:** "which partner fell" is a property of the perturbation,
not of a sweep; the first solve answers it without the cap feeding back into
its own condition. No threshold is read during the iteration, so no switch
can chatter. SMAD2 KO: the source copy is pinned to 0 in both passes, so the
pair is gated, as before. MDM2 KO: the source (excluding the target) does not
fall, so it is free, as before.

**Arm `rtwin58g4`:** the adopted defaults at the amended commit, against the
control `64d23cb` (pre-adoption) and against `9e93c11` (adopted). Same gate,
plus a convergence condition registered now: **curator and experimental
convergence must be at least the control's (1,652 / 215) minus 5**, and every
discordant case is re-checked at caps 500–4000 with `capcheck.jl`; only stable
cases are counted as fixed or broken in the decision (both counts reported).

**Predictions:** held-out keeps the 33 stable fixes (SMAD2/3, MET, CBL, PTCH
classes do not depend on the toggling pairs); SMAD7 KO converges and returns
to UP (the control's answer), so held-out about +33. The TP53 AKT/MDM4
knockouts converge to one answer; I do not know which, so tuning and
experimental could move either way from the control by up to the 111 / 6
unstable cases. Convergence returns to about the control's. If the stable
held-out gain survives and convergence holds, the default stays (amended);
if not, the adoption is reverted to the previous defaults.

## Amendment 4 (2026-10-15, after arm `rtwin58g4`, before any amended code): no cap is applied during the iteration

`rtwin58g4` (amendment 3) moved no class against the adopted defaults and
convergence stayed at 1,645. Per-case convergence (new in the case dumps)
showed TGF-beta SMAD7 KO now converges, but TP53 still oscillates. Traced
(`bench/analysis/061/tp53iso.jl`, AKT1 KO, readout 6801636): `full`
converges (residual ~1e-6) and so does `suppress_only_entity`; `suppress_only`
alone oscillates (residual 0.03), with or without twins, under either pin
set. So the cause is the node-copy cap itself, `min(factor, 1)` on a depleter
built directly from its target: a kink at baseline that the damped
Gauss-Seidel sweep circles. specs/060's TP53 tuning gains under
`suppress_only` were probably the same artifact.

**Rule (`suppress_only_source`, changed in place; amendment 3's form is at
`b46783c`; the other `suppress_only*` modes are left as recorded):**
1. Pass 1: solve with no own-product cap at all (as `full`).
2. A depletion edge D -| T is REMOVED for pass 2 if it would de-repress
   (D reads below baseline × (1 − 1e-6) in pass 1) and either
   (a) D is a node-copy own product of T (T reaches D within two activator
   hops, the existing structural test), or
   (b) D is an entity-own depleter of T (containment) and one of its source
   copies (amendment 1–2 walk) reads below baseline in pass 1.
3. Pass 2: solve with those edges removed and no cap (as `full`).
No threshold or `min` is evaluated during either iteration. Removing an edge
is what the cap does at a fixed point for a de-repressing edge (factor 1).
The solve reports `own_source_depletions` as the edges removed.

**Arm `rtwin58g5`:** the defaults at the amended commit, against the control
`64d23cb` and the adopted `9e93c11`; the amendment-3 convergence condition
(curator ≥ 1,647 and experimental ≥ 210 converged) and the cap-stability
re-check (`capcheck.jl` on every discordant case) apply, and only stable
cases count in the decision.

**Predictions:** convergence back to about the control's (TP53 and TGF-beta
settle). The 33 stable held-out fixes stay (they do not depend on TP53).
TP53's AKT/MDM4 knockouts converge to one answer; most likely their
control-like NOT, so tuning falls from +112 to near the 2 stable fixes and
experimental to about +1. If the held-out gain holds with convergence, the
default stays (amended); otherwise it is reverted.

## Amendment 5 (2026-10-15, after a probe of amendment 4, before the arm): drop the node-copy cap

Before running `rtwin58g5`, the TP53 probe (`tp53iso.jl`, AKT1 KO) under
amendment 4 still oscillates (residual 0.04 at caps 500 and 2000): removing
the de-repressing node-copy edges leaves TP53's loop without the restoring
force, so pass 2 does not settle either. Under `full` it converges. So the
node-copy part of the rule (inherited from `suppress_only`, specs/016/060)
prevents convergence in TP53 in every form tried: a live `min`, and removal.
specs/060 R1 already showed it moves **no** held-out case (its tuning +46 were
all TP53, and R5 shows the TP53 movements are cap artifacts).

**Rule (`suppress_only_source`, changed in place; amendment 4's form was
probed in the working tree, not committed or run as an arm):** as amendment 4, but node-copy own-product
pairs are NOT capped or removed (they act as under `full`). Only entity-own
pairs whose source copy fell, and whose depleter fell, in the uncapped first
solve are removed for the second. `rtwin58g5` is not run; arm `rtwin58g6`
replaces it, with the same comparisons, convergence condition and cap
re-check.

**Predictions:** convergence about the control's (TP53 behaves as under
`full`); held-out keeps the 33 stable fixes (SMAD2/3 through the entity-own
SMAD2:SMURF2 pairs, MET/CBL/PTCH through pins and twins); tuning about +2,
experimental about +1. If held-out stays positive at p < 0.01 with
convergence, the default is this amended form; otherwise the adoption is
reverted.
