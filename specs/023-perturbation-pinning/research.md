# 023 — Pin where the gene enters the network, and let the rest propagate

## What the benchmark does today (measured 2026-09-25, build `20260925-1039_d4f4f64`)

A perturbed gene resolves to the entities that reference it (Neo4j
`referenceEntity`, `neo4j_gene_to_stids`). **Since 2026-07-14 (`0bd4565`)** each
of those stids also matches every node whose `member_leaves` include it, so the
benchmark pins **every complex that contains the gene, wherever it sits**:

- **Scale:** 856 perturbations pin 9,378 nodes. Only **1,023 (11%) are roots**
  (in-degree 0).
- **Mid-pathway pins:** 8,355, made up of 6,371 complexes, 1,333 dissociation
  sinks and 651 entity forms produced by a reaction (e.g. a glycosylated
  secreted form).
- **Coverage:** 805 of 856 perturbations pin at least one mid-pathway node. 77
  have no root occurrence at all.
- **Example:** WNT5A in Signaling by WNT pins 1 root entity and 36 downstream
  nodes. These include generic set complexes (`WLS:WNT`, `WIF1:WNT`), and
  pinning them at 0 knocks out *every* WNT ligand. That is the specs/021
  reversal.

A hard pin overrides the node's own producing reaction. A downstream complex
is then not *computed* from the perturbed gene; it is *set*, including for
partners the knockdown should not affect.

## Intended protocol (Adam, 2026-09-25)

*"the idea is to pin the root input that is either a particular entity or a
complex containing that entity … breaking apart the complexes into root
entities where the root would look like A entity and B entity -> AB complex ->
reaction. In this way both the entity that is not part of the complex and the
entity inside the complex can be perturbed."*

## The rule (`DS_PIN_SCOPE=entry`; default `all` = today)

Let G be the nodes the gene resolves to today. Pin only the **entry
occurrences**: the members of G that cannot be reached from any other member
of G. Roots in G are always entries.
- A gene first produced mid-pathway (e.g. by transcription) gets that first
  occurrence pinned.
- When every member of G is reachable from another (G lies inside one cycle),
  all of G is pinned, as today.
- Everything else propagates.

## Pre-registration (committed before the arm runs)

- **Arms:** control = production scoring (`ae84de9`, byte-identical to the
  `ffb3aa1` control in specs/022); arm = `DS_PIN_SCOPE=entry`, code defaults
  otherwise.
- **Decides:** curator held-out net and macro-F1, McNemar p, pathways moved,
  distinct genes. **Reported:** tuning, experimental axis, per-pathway net,
  pinned-node counts, and the WNT5A dose ladder.
- **This is a protocol correction, not a model change.** The primary question
  is whether the benchmark measures what it claims to, so it is adopted if the
  held-out result is not worse than the noise floor (−15). Adam's statement of
  intent is the justification, not the accuracy. A gain is reported, not
  required. A loss beyond the noise floor is traced before any decision.
- **Prediction:** large churn, because 94% of perturbations change their pin
  set. Direction unknown. The false-change rate should fall, because
  downstream complexes stop being set directly.

## Result, arm 1: `DS_PIN_SCOPE=entry` alone — a large loss, traced

Build `20260925-1039_d4f4f64`, bench `ec81afd`, solver unchanged; control =
production `ae84de9`. Same 23,268 valid cases in both arms.

| | control | entry pinning | net | fixed / broke | p |
|---|---|---|---|---|---|
| curator held-out | 0.8687 / mF1 0.8299 | 0.8281 / mF1 0.7617 | | | |
| curator all | 0.8477 / 0.8145 | 0.7994 / 0.7401 | **−1,164** | 179 / 1,343 | 1e-220 |
| experimental | | | **−137** | 3 / 140 | 9e-38 |

- **Breadth:** 71 pathways got worse and 2 got better.
- **Direction:** 1,149 of the 1,343 breaks (86%) are overexpression cases,
  and 1,148 are UP → NORMAL.

**Traced (KMT2C OE, Chromatin modifying enzymes; readout 100 → 1.0).**
- KMT2C, ASH2L, RBBP5, WDR5 and DPY30 feed the MLL3 complex by `assembly`
  edges. This is exactly the decomposed-root structure the protocol intends.
- The old protocol pinned the **MLL3 complex itself** at 80x.
- Under entry pinning only KMT2C is pinned. `DS_ASSEMBLY_LIMITING=1` makes a
  complex its scarcest subunit: min(80, 1, 1, 1, 1) = 1. The overexpression
  cannot pass through any complex.
- **So the broad pin was compensating for a modelling choice.** It set every
  complex directly, which is how 89% of pins came to be mid-pathway, and how
  overexpression cases scored.

**Status:** the pinning is now what the protocol intends, but the model cannot
carry a single-subunit overexpression through an assembly. Adopting arm 1 alone
would make the benchmark honest and the model worse at the thing curators score
most.

## Pre-registration, arm 2 (committed before it runs)

`DS_PIN_SCOPE=entry` + `DS_ASSEMBLY_LIMITING=0`, so assembly inputs combine by
the AND mode (`hill_sat`, multiplication capped at 100) instead of min.
- A knockdown still pulls a complex down (0 × anything = 0).
- An overexpressed subunit now raises it.
- The specs/009 measurement of limiting off (−196) was taken under broad
  pinning, which masked exactly this effect, so it does not carry over.

**Reading:** compare against the production control on the same columns as
arm 1. If held-out is within the noise floor (±15) of production or better,
the honest protocol plus multiplication is adoptable as a pair, pending Adam's
call on the biology. That call is whether overexpressing one subunit should
raise a complex when its partners are at baseline. If held-out is below −15,
record it and trace.

## Result, arm 2: entry pinning + `DS_ASSEMBLY_LIMITING=0`

Same build, bench `ec81afd`. The container's own environment was verified as
`DS_ASSEMBLY_LIMITING=0` before the run.

| | production | arm 1 (entry) | **arm 2 (entry + multiply)** |
|---|---|---|---|
| curator held-out acc / mF1 | 0.8687 / 0.8299 | 0.8281 / 0.7617 | **0.8475 / 0.7952** |
| curator all acc / mF1 | 0.8477 / 0.8145 | 0.7994 / 0.7401 | 0.8281 / 0.7846 |
| held-out net vs production | — | | **−402** (88 / 490), 36 of 38 moved pathways worse |
| tuning net | — | | −69 |
| experimental net | — | −137 | **−51** (1 / 52) |

Letting an overexpressed subunit raise its complex recovers about 60% of arm
1's loss, but the honest protocol still scores well below production on both
axes. The remaining breaks are 339 OE and 151 KD on held-out. They are not yet
traced.

**What this means for every number reported since 2026-07-14.** Production
accuracy (held-out 86.87%) is measured under a protocol that sets every complex
containing the perturbed gene directly. Part of that accuracy is the pin, not
the propagation. Until the protocol question is decided, **quote the broad-pin
numbers with that caveat.** This affects docs/RESULTS.md and the manuscript
claims.

**Decision needed (Adam):**
1. **Which protocol is the benchmark of record?**
   - (a) Broad pins, as today: comparable to every result since July, but it
     sets complexes by hand.
   - (b) Entry pins: what the protocol intends, but about 2pp worse today, and
     it exposes the modelling gaps.
2. **Under (b), should an overexpressed subunit raise its complex** when its
   partners are at baseline? Limiting-reactant biology says no; the curators'
   expectations say yes.

Nothing is adopted until then. The flag stays, default `all`.

## Decision (Adam, 2026-09-25): root pinning is the protocol of record

*"These tests from the mpbiopath publication perturb only root inputs that are
or contain an entity. So that is what we should go with."*

`DS_PIN_SCOPE=root` is now the benchmark default. It pins only nodes with no
incoming edge that are, or contain, the gene. This is narrower than `entry`: a
gene with no root form is not perturbed, and its cases become invalid rather
than scored (77 perturbations had no root occurrence). `all` is kept only to
reproduce results from 2026-07-14 to 2026-09-25.

## Pre-registration: the new baseline and the assembly question under it

The protocol is decided; it is not being A/B'd. What is measured:

- **Arm R (new baseline):** root pinning, solver code defaults. Its numbers
  replace the production baseline in docs/RESULTS.md, reported against the old
  broad-pin baseline for continuity, **on the shared valid set**, with the count
  of cases that became invalid.
- **Arm R0:** root pinning + `DS_ASSEMBLY_LIMITING=0`. This asks whether a
  single overexpressed subunit should raise its complex, a modelling question
  the broad pins had hidden. It is decided against arm R on curator held-out
  (net, McNemar, concentration), with the experimental axis reported.
  - Adopt limiting-off as the solver default only if held-out net > +15 with
    p < 0.05, not concentrated, and experimental no worse than −15.
  - Otherwise keep the default and record.
- Both run through `scripts/run_arm.sh`.

## Results under root pinning (build `20260925-1039_d4f4f64`, all through `scripts/run_arm.sh`)

All comparisons are paired on the cases valid in both arms
(`bench/analysis/arm_compare.py`). Root pinning makes 1,144 curator cases
invalid (a gene with no root form), so accuracies here are on the valid set and
are not comparable to earlier holdout_report figures, which scored invalid
cases as NORMAL.

| arm (commit / name) | curator held-out acc / mF1 | vs previous row, held-out net | experimental acc / mF1 |
|---|---|---|---|
| broad pins, production (`ae84de9`) | 0.8795 / 0.8452 on its own 18,238; **0.8816 / 0.8428 on the shared 17,420** | — | 0.7224 / 0.6501 |
| root pins, code defaults (`1eca749/root_baseline`) | 0.8352 / 0.7625 (n 17,420) | −809 (134 / 943) | 0.5528 / 0.5163 |
| + `DS_ASSEMBLY_LIMITING=0` (`1eca749/root_nolimit`) | 0.8582 / 0.8031 | **+401** (565 / 164), p 2e-52, 50 of 58 pathways up | 0.6622 / 0.5809 (+89) |
| + `DS_SELF_INHIBITOR_WEIGHT=0.1` (`1a10ed2/root_nolimit_selfinh_w0.1`) | **0.8634 / 0.8135** | **+90** (127 / 37), p 1e-12, 10 of 11 up | **0.6732 / 0.5896** (+9) |

Checks:
- The merged commit's root baseline (`bf6cfc1`) reproduces the pre-merge one
  byte-for-byte on the scored columns.
- w = 0.1 with limiting on (`bf6cfc1/root_selfinh_w0.1`) is +59 held-out, so
  the weight helps under either assembly rule.

**Concentration of the w = 0.1 gain.**
- The top held-out pathway is RUNX1 (+33).
- Without it, held-out is **+57** (91 / 34, p 3.5e-7) over 24 perturbations.
- The next largest are Pre-NOTCH (+23), IL-4/13 (+11), ROCKs (+8) and AP-2
  (+8).
- It passes the specs/022 rule that failed under broad pins, where it was one
  gene.

**Adam predicted this**, as *"this change we made earlier should have fixed a
bunch of cases"*. It did not show under broad pins because the pins set the
very complexes whose inhibitors the rule corrects.

**Adopted as solver defaults under the root protocol:**
- `DS_ASSEMBLY_LIMITING=0`, by the arm 2 rule above;
- `DS_SELF_INHIBITOR_WEIGHT=0.1`, by the specs/022 rule, now met.

Both axes improve at every step after the protocol change. **On identical
cases, the new defaults are below the old broad-pin numbers**: held-out
0.8634 vs 0.8816 on the same 17,420, net **−318** (104 / 422), and
experimental −40 (9 / 49). Part of the old accuracy came from the pins. Two
further caveats (review of PR #72):
- **Survivorship.** The 1,144 cases lost to invalidity were 79% correct under
  the broad pins (908 / 1,144), so valid-only accuracy flatters. The
  every-case row in docs/RESULTS.md is the honest one.
- **Concentration.** 33 of the 37 held-out breaks from w = 0.1 are in RUNX2,
  which alone is 31 fixed / 33 broke (66 changed scored predictions). Excluding
  both RUNX1 and RUNX2 leaves 60 fixed / 1 broke. The McNemar p treats about
  25 perturbations' correlated readouts as independent, so it is optimistic.

## Where the adopted configuration still fails (`bench/analysis/failure_structure.py`)

Arm `1a10ed2/root_nolimit_selfinh_w0.1`, which is equivalent to the new
defaults. Held-out: 17,420 scored, 2,380 wrong (13.7%): 1,403 missed, 770 false
change, 207 wrong direction. Each case is classified on the region between its
pinned roots and its readout.

**1. No route at all: 1,137 missed changes (48% of held-out errors).**
- Only 122 of them (11%) were right under the old broad pins, so this class
  mostly predates the protocol change.
- It is concentrated: Interferon α/β 256 (the severed branch of specs/016 and
  019), MET 80, PDGF 80, RUNX1 78.
- It spans 267 perturbations and 250 readouts.
- This is connectivity, not propagation; a propagator change cannot reach it.

**2. With a route (5,712 cases, 21.8% wrong).** Loop and structure effects
are given within pathway, because pooled loop effects have been
between-pathway confounds before (the retracted "loops are the lever").

The estimator is the committed Mantel-Haenszel risk difference in
`failure_structure.py` (`mh_risk_difference`), run on
`d4f2bda/new_defaults`. An earlier version of this table used ad hoc code
whose values depended on the estimator; the review of PR #72 found the
negative-loop sign flipped between estimators.

| on the route | pooled error rate | MH within-pathway difference | pathways |
|---|---|---|---|
| any loop | 28.4% vs 14.8% | +7.5% | 39 |
| loop with only activating edges | 28.1% vs 14.8% | **+8.0%** | 30 |
| loop containing any inhibitory edge | 28.7% vs 14.8% | +2.2% | 16 |
| giant loop (≥ 100 nodes) | 28.6% vs 14.8% | +9.1% | 9 |
| welded loop (only derived edges close it) | 24.2% vs 14.8% | +10.4% | 6 |
| self-contained inhibitor (solver's definition) | 25.1% vs 21.0% | **+12.2%** | 18 |
| assembly step | 25.1% vs 18.4% | +7.7% | 62 |

- **All-positive loops carry most of the loop excess within pathway.** Loops
  containing an inhibition carry a smaller one (+2.2%). "Contains an
  inhibitory edge" is not the same as negative feedback: the label is now
  `has_inhibition`.
- False change dominates the failures with a route.
- **Self-contained inhibitors remain the largest structural excess** after
  w = 0.1.

**Traced: CREBBP knockdown, DDX58/IFIH1 interferon induction** (IFNB1
expected DOWN, predicted 2.75x):
- *CREBBP, EP300 binds p-IRF3 dimer* has two variants. The CREBBP variant
  correctly reads 0.
- The EP300 variant's negative regulator is **CREBBP:NS1**, influenza NS1
  sequestering CREBBP. It contains the knocked-down gene but not this variant's
  input (EP300), so it is not self-contained per variant. It falls to 0 and
  de-represses the EP300 variant to the 10x ceiling.
- OR-mean of the two variants gives 5; with the IRF7 branch, the readout
  reads 2.75.
- **The mechanism: a sequestration complex falls because the sequestered
  protein fell, and that is read as the sequestering agent disappearing.** Here
  it happens across sibling variants of one Reactome reaction, which the
  per-variant rule cannot see.
- It is the same class as WNT3A:sFRP in the WNT5A trace.
- A second issue in this case: **NS1 is a viral protein** held at "normal"
  level in an uninfected human pathway.

**Candidate next rule (not implemented):** for a complex inhibitor, only the
change in its *non-shared* components (the sequestering agent, NS1 or WIF1)
should count as inhibition. A drop caused by the sequestered partner should
not de-repress. That is the specs/022 split generalised from "the reaction's
own input" to "any component that is an input of the same Reactome reaction".
Two things are needed first:
- a definition of the shared fold when that component is not this variant's
  input;
- Adam's call on whether to go by reaction or by variant.


## The self-inhibitor rule after two reviews of PR #72

The rule was revised twice after it was measured. Each revision is
re-measured through `scripts/run_arm.sh` on the same build:

| version | what changed | commit | vs the version before |
|---|---|---|---|
| 1 | inhibitor fold divided by the shared input's fold | `1a10ed2` / `44e733b` | (the +90 above) |
| 2 | + requires the input to reach the inhibitor; + weaken-only clamp; + 1e-3 fold floor | `d4f2bda` | **0 scored changes**; 213 raw values moved |
| 3 | split by the **log-fold overlap** of inhibitor and input | `5f7b9c7` | **0 scored changes**; 173 raw values moved (max 1.5x) |

- **Why version 3 exists.** Version 1 amplified where the inhibitor did not
  track the input (L = 5x read 21x). Version 2 fixed that, but where an
  inhibitor tracked only partly, its clamp deleted the inhibitor outright
  (722 of 1,214 probed catalog cases, 588 in PIP3). That is the edge deletion
  specs/012 measured as harmful, not Adam's rule.
- **What version 3 does.** It attributes to the input only the overlap of the
  two changes, keeps that at weight w, and keeps the rest whole. A fully
  tracking inhibitor still reads x^(1-w).
- **What the measurement means.** All three versions give identical class
  predictions on this catalog. The scored gain (+90 held-out over limiting-off)
  comes entirely from fully tracking inhibitors, where the versions agree. The
  revisions change behaviour off the benchmark's cases: partial tracking,
  non-tracking inhibitors, OR inputs, and knockdowns near 0. That behaviour is
  pinned by `test/test_self_inhibition.jl` (80 assertions).
- **Canonical numbers are unchanged.**
