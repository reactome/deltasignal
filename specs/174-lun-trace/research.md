# specs/174 research (2026-11-26)

Registration `fab220a5`. Code: `bench/analysis/174/trace_174.py`, `k562_174.py` (tests `test_174.py`). Per-case
tables, traces and counts are outside the repo (`~/deltasignal-catalogs/analysis/174/`). Lun half B was not read.

## Lun 2019 half A: classes (registered)

3,127 measured pairs, the specs/173 call and predictions unchanged.

| class | | cases | on background readouts |
|---|---|---|---|
| R1 right change | route, model calls the measured direction | 143 | 86 |
| R2 wrong direction | route, opposite direction | 40 | 28 |
| R3 called, silent | route, model calls a change, nothing moves | 365 | 116 |
| R4 missed on route | route, model NOT, readout moves | 61 | 42 |
| R5 right no-change | route, both NOT | 183 | 82 |
| N1 changed, no route | | 500 | 331 |
| N2 no change, no route | | 1,835 | 661 |

Readout types (registered rule): 10 **background** readouts (AKT, AMPKα, ERK1/2, STAT3, STAT5, p38, p90RSK, FAK, SRC,
p70S6K: each changes for 25–43% of proteins with no route to it), 4 **route-specific** (MAPKAPK2, MKK3/6, MKK4,
PLCγ2), 14 quiet.

**Q1. 60% of the right change calls (86 of 143) are on background readouts**, above the registered half: most of the
model's right calls are on readouts that move under many unrelated proteins. Only 29 right calls are on route-specific
readouts.

**Q2. Failures are spread out.** The five readouts with the most wrong change calls (β-catenin 44, p-AKT 44, p-GSK3β
37, p-MAPKAPK2 36, p-CREB/ATF1 35) hold 48% of the 405 R2 + R3 cases.

## Traces (registered groups)

Shortest signed path from the perturbed protein to the readout in each network that reaches it (readout nodes matched
from `leaf_refs.csv`; p70S6K's and CREB's sites did not match by position, so those groups are partly untraced).
Mechanism labels are hypotheses from the trace.

| class: readout | cases | what the trace shows | label |
|---|---|---|---|
| R3: β-catenin | 40 | The model reads total β-catenin as every complex containing it (Lun predictions were made on 2026-11-03, one day before specs/134 made total readouts read abundance). Routes end in AKT2:β-catenin or destruction-complex forms, about 2×. | **readout** |
| R3: p-GSK3β | 35 | All 31 traced go through PIP3 → PDPK1 → AKT2 → GSK3β. 25 read **45×** whatever the perturbation, including paths with an odd number of negative steps; 6 read about 0.5–0.75×. That is the PIP3 switch (specs/136) railing. | **other: PIP3 switch** |
| R3: p-AKT | 28 | 20 of 22 traced go through PIP3; the switch rails at 46–58×. Paths are long (median 22 steps) and cross modules: BRAF → ERK → FGFR complexes → PI3K; ALK → TGFBR3 → FGF. | **other: PIP3 switch** |
| R2: p-AKT | 16 | Receptor tyrosine kinases (EGFR, NTRK1, FLT3, FLT4, ALK), SRC, BTK, PI3K subunits (PIK3CG, PIK3R3) and PIP5K1A: the model reads 47–58× UP, measured DOWN. The curated route exists and is short (6–16 steps); nothing curated lets sustained over-expression lower basal p-AKT (feedback, regulatory-subunit excess). | **missing negative**, amplified by the switch |
| R3: p-CREB/ATF1 | 35 | 11 traced (the CREB1 site did not match by position; ATF1's did): the curated ERK/p38 → MSK1 (RPS6KA5) → ATF1 route, 4–36 steps, about 4–5×. Measured: no change, even for MSK1 itself (5.7×). | **inactive in cell** |
| R3: p-MAPKAPK2 | 35 | Of 31 traced: 13 run through T-cell or B-cell receptor modules (PRKCQ → CARD11 → … → p38 → MAPKAPK2, 42–60 steps) and 3 through VEGF receptor complexes; 12 by other long routes (FGFR/GRB2 complexes, 24–37 steps); 3 are short (MAP2K6, TAK1 complex, MAPK14: 4–10 steps) and read 14–100×, but MAPKAPK2 does not move. | **inactive in cell** (immune routes in HEK293T); short ones **other** |
| R3: p-ERK1/2 | 21 | Of 20 traced, 8 go through PIP3 → PDPK1:PKCζ → ERK; most others through RAF/KSR scaffold complexes (the RAF switch, specs/054), 2–40 steps. | **other: switch** |
| R4: STAT5, p70S6K, AMPKα | 31 | Routes are very long (median 47 steps) and the model reads about 1×: the signal fades to nothing along the path. 24 of the 61 R4 cases come from the hub proteins below. | **other: background** |
| N1: AKT, STAT5, ERK, SRC, FAK, p90RSK, STAT3 | 256 | No route. Mostly UP (STAT5 38 of 42, p90RSK 31 of 31), except p-AKT, which falls (32 of 48). | **other: background** |

**Passes on route-specific readouts (R1: 16 cases in the registered groups, 12 traced):** short, canonical kinase cascades, all in curated immune or MAPK
modules. SYK → PLCγ2 (4 steps, 99×); BTK, ITK, TEC, LCK → PLCγ2 signalling complexes (6–8 steps); TAK1 (MAP3K7),
IRAK2, ASK1 (MAP3K5) → MKK6 → p38 → MAPKAPK2 (8–10 steps). NTRK1 → … → MAPKAPK2 (24 steps) and CSF1R/SRC by longer
routes are also right.

**Route length (post hoc, all 792 pairs with a route; 591 traced):**

| steps (entity and reaction nodes) | pairs | changed | model change calls | right direction |
|---|---|---|---|---|
| ≤ 6 | 106 | 33% | 98 | 26% |
| 7–12 | 80 | 39% | 73 | 33% |
| 13–24 | 142 | 23% | 113 | 15% |
| > 24 | 263 | 27% | 135 | 24% |

Without background readouts: ≤ 12 steps 26% change, > 12 steps 17%, no route 13%. Long routes are close to no route;
short routes carry some signal; but even on short routes the model's direction is right on only 26–33% of its calls.

**Hub proteins (post hoc; the definition uses outcomes, so this is descriptive only).** 22 of 134 proteins change at
least 75% of the background readouts they are measured on: ALK, BTK, CSF1R, EGFR, EPHA2, EPHA4, FES, FGFR1–3, FLT3,
ITK, JAK3, MAP3K8, NTRK1, PPP1CC, RET, RPS6KA6, SRC, STK3, TAOK3, TEC (mostly tyrosine kinases). They change 62% of
their readouts with a route and 61% without. Route does not matter for them. They account for 176 of the 500 N1 cases
and 67 of the 143 right change calls. Without them, pairs with a route change 22% and pairs without 16%.

## K562 split A: curated regulator pairs (registered check)

Served Release 98 build `20261113-1311_baf608b_r98fix`; symbols mapped by the reviewed human proteome FASTA. 80 scored
(knocked-down regulator, target) pairs from 34 genes (specs/169's Neo4j version found 104; the two definitions differ
in what counts as an expression event). **2 of 80 changed ≥ 1.5×** (STAT5B, E2F7).

Of the 78 that did not change:
- 19 come from knockdowns that change none of their scored targets (SMAD4, MOV10, NANOG, RBPJ, ESR1, TCF7L2 and others):
  the regulator is not active, or not needed, in K562;
- the other 59 come from knockdowns that do act. Median knockdown effect 6.1% of scored targets, against 0% for a
  typical gene, and knockdown efficiency is good (median 22% left). These regulators act in K562, but not on their
  curated target. Most are broad co-repressors and co-activators or RNA-silencing machinery (RBBP4, HDAC1, RBL1,
  CREBBP, AGO2/3, TNRC6A), where redundancy and context decide the target.

## What this means

1. **Much of the apparent skill is hub proteins and background readouts.** Over-expressed tyrosine kinases switch on
   most of the panel regardless of route, and 60% of the model's right change calls are on readouts that move under
   many unrelated proteins.
2. **Real route signal exists, but it is narrow:** short, canonical kinase cascades (SYK/BTK/TEC → PLCγ2; TAK1/ASK1 →
   MKK6 → p38 → MAPKAPK2). Long routes (over about 12 steps) carry almost nothing.
3. **The wrong calls have four sources, none of them fixed by adding edges:**
   - the PIP3 switch rails AKT, GSK3β and part of ERK at 45–58× under almost anything upstream (specs/136, known);
   - curated routes that are silent in these cells (MSK1 → CREB/ATF1; immune receptor modules in HEK293T);
   - no curated way for sustained over-expression to lower AKT;
   - a readout defect: total β-catenin is read as all its complexes. The specs/134 abundance rule was adopted the day
     after the Lun predictions were made.
4. **K562:** the curated regulators mostly act in the cell but not on the curated target. That is context, not a
   missing route.

None of this is adopted. Candidate registrations, each to be judged once on half B:
- re-predict with the abundance rule (a readout fix that was already adopted, no fitting);
- a route-length or route-type limit on change calls;
- a hub-protein term.
