# specs/061: an entity-own depleter may de-repress only while its source copy holds

**Status:** registered 2026-10-14, before any code or arm. R1 (`rtwin58g`) and R2 (`ownsrc`) fail: the source test was circular when the target is its own source. **Amendment 1 (`rtwin58g2`, R3) passes the registered gate:** held-out +32 (33 / 1), tuning +112, experimental +7 (7 / 0), worst pathway −1. Not adopted yet: post hoc, so an independent check is recommended first (Adam's call). Post hoc: the
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
