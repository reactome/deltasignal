# specs/133 research

## R1: the closure (census, no solving; 2026-11-03)

`bench/analysis/closure_133.py` (Neo4j only) follows specs/132's L1 handoffs back from the 92 scored pathways until a
round adds nothing. Output (private): `~/deltasignal-catalogs/analysis/133/` (`links.tsv`, `s_q.tsv`, `closure.tsv`).

- **275 pathways** are reached: the 92, then 83, 50, 22, 15, 7, 2, 1, 1 and 2 by round (nine rounds). 216 of them
  receive at least one handoff; 2,433 handoff rows over 407 entities (input 1,952, catalyst 352, positive regulator
  115, negative regulator 14). These rounds are the superset the builder draws from; each scored pathway keeps only
  its own tree.
- **Extracellular, never bridged: 26 rows, 18 entities**, every one a signal between cells or a solute outside them:
  serotonin, noradrenaline, L-glutamate, choline, LTC4, bile salts, complement C3a and C5a, angiotensin fragments and
  alamandine, GLP-1, acyl ghrelin, IL4, IL21; and H+, Na+, H2O2 outside the cell. The rule removes the IL4 class
  without naming it.
- 33 of the 407 entities are sets (79 members; `set_members_132.py` on this table).
- 100 of the 275 pathways were not in specs/132's base build; the base is rebuilt whole with the frozen generator
  (`8b6527e`, tag `freeze-2026-10-26`) as variant `ctx133`.

## R2: the bundles (census, no solving; 2026-11-03)

Base build `20261103-2155_8b6527e_ctx133` (275 pathways, all built; not `current`). Bundles from it with
`bench/context_build.py` (commit 910c5a7), each over the 92 scored pathways:

| bundle | rule | pathways with a bridge | bridges |
|---|---|---|---|
| `.c0` | no links (each pathway alone; the control) | 0 | 0 |
| `.A` | specs/132 L1, extracellular rows removed | 25 | 72 |
| `.T` | `--transitive` (the tree) | 25 | 163 |

- **A is specs/132 L1 minus exactly the two cell-to-cell bridges** (L1: 27 pathways, 74 bridges): IL4 into *Interleukin-4
  and Interleukin-13 signaling* and IL21 into *Interleukin-2 family signaling*. Nothing else in the one-hop set is
  extracellular.
- **T joins the same 25 pathways, not more.** A pathway only gets a tree when one of its own roots receives a handoff,
  and the closure's depth adds sources to those, never new receivers. Of 2,109 handoff rows the builder examined, 270
  are bridged (by depth: 115, 76, 54, 13, 12); 1,540 land on a copy that is not a root in the receiving pathway (the
  entity is there but made inside it), 232 come from a pathway already joined at a lesser depth, 59 carry a held entity
  (cofactor, drug, pathogen), 6 are extracellular (4 at depth 5 in *Costimulation by the CD28 family*'s tree, and the
  IL4/IL21 rows).
- **Depth is shallow for most.** 16 of the 25 pathways stop at depth 1 (T equals A there). The deep trees: *Transcriptional
  Regulation by TP53* (28 pathways joined, depth 5), *PIP3 activates AKT signaling* (11, 3), *Intrinsic Pathway for
  Apoptosis* (10, 4), *Signaling by ERBB2*, *S Phase*, *RUNX2* (9 each, 4), *Mitotic G1-G1/S* (8, 4), *CD28 family
  costimulation* (7, 4).

So "all of Reactome" in practice is these nine trees; the other 83 scored pathways answer as they do alone. The ceiling
(cases with no route alone and a route through T) is read from the arm results (R3).

## R3: the first runs of A and T, and the defect they exposed (2026-11-04)

Control `s133_c0` (build `…ctx133.c0`, code 910c5a7) reproduces specs/132's control case for case on every axis
(experimental 796, curator held-out 18,481, korkut 11,062, decryptM 1,802 cases: 0 predictions and 0 routes differ), so
the rebuilt 275-pathway base changed nothing for the 92 scored pathways.

**Reporting gap.** The API passed solver diagnostics through a fixed list that lacked `context_bridge_order`; the first
runs' logs read "order unknown". The containers did run `DS_CONTEXT_BRIDGE_ORDER=last` (`server.env`), and the case
changes below confirm it. Fixed in 4800ab3.

**Q1 holds.** Against specs/132 L1 (same bridges into Cell Cycle Checkpoints, order `first`), arm A returns ATM's four
experimental cases to the correct answer (knockout DOWN, overexpression UP at both readouts). **Q2 holds**: the 91
curator and 34 korkut cases of *Interleukin-4 and Interleukin-13 signaling* that L1 moved are back.

**The defect.** A also moved five PTEN overexpression cases in *Transcriptional Regulation by TP53* away from L1, four to
wrong. Traced by solving the pathway under both orders (scratch Julia; rule edges compared): the only rule difference is
three recycled-twin edges (specs/059), from TP53's own bridged copies of p-AKT1/2/3 [cytosol] (R-HSA-198356, -202062,
-3009367) to the p-AKT copies imported from *PIP3 activates AKT signaling* that feed the bridges. With the bridges set
aside the imported copies look recycled, and the rule, which matches copies by entity alone, joined them to TP53's.
Restored, that is a loop PIP3 → TP53 → PIP3: a shared molecule joining two pathways, which specs/132's principle forbids,
and feedback into the context, which this spec's tree forbids. Registered as amendment 1 (spec.md) before the
first runs' scores were read; implemented in 4800ab3: a bundle records each node's pathway (`context_origin.csv`) and,
under `last`, the solver drops every rule edge between two pathways.

`bench/analysis/origin_check_133.jl` checks that the drop equals running the structural rules on each pathway's part
alone. It does, on all 25 bridged pathways of both bundles (0 differ). The edges the rules had drawn across pathways were
not few: Cell Cycle Checkpoints 140 (both bundles), TP53 transcription 46 under A and 116 under T, 0–9 elsewhere. Under
`first` the same identity joins can form from any unbridged starting copy, so specs/056 and specs/132 arms likely held
some too; not measured, since neither is a candidate.

First-run scores, kept as the record of the defect, not as results (rule+F against c0; gate as registered):

| arm | experimental fixed / broken | curator held-out | korkut fixed / broken | decryptM | pooled gain [95%] |
|---|---|---|---|---|---|
| A (defective) | 5 / 0 | 12 / 4 | 69 / 310 | 0 / 22 | +0.0115 [+0.0014, +0.0225] |
| T (defective) | 5 / 0 | 12 / 6 | 63 / 312 | 0 / 22 | +0.0097 [−0.0001, +0.0208] |

The corrected arms are `s133_A1` and `s133_T1` (bundles `.A1`, `.T1`, identical to `.A`, `.T` but for
`context_origin.csv`; code 3541e6a). The control is not re-run: c0 has no bridges and no origin table, so the amendment
cannot change it.

## R4: A1 and T1 (amendment 1), and the source pathways that were not themselves (2026-11-04)

Arms `s133_A1`, `s133_T1` (code 3541e6a), against c0. Logs: "order last", 3,908–7,490 cross-pathway rule edges dropped
over all solves. Rule + F (the decision):

| arm | experimental fixed / broken / wrong→wrong | curator held-out | korkut | decryptM | pooled gain [95%] |
|---|---|---|---|---|---|
| A1+F | 5 / 0 / 10 | 12 / 4 | 71 / 250 | 0 / 22 | +0.0141 [+0.0037, +0.0254] |
| T1+F | 5 / 0 / 10 | 12 / 6 | 71 / 250 | 0 / 22 | +0.0141 [+0.0037, +0.0254] |

Accuracy / macro-F1 (rule+F minus c0): experimental +0.0063 / +0.0064, curator +0.0003 / +0.0005, korkut −0.0162 /
+0.0008, decryptM −0.0122 / −0.0065. The gate passes on balanced accuracy while korkut's raw accuracy falls: by unique
(drug set, readout) korkut is 30 fixed / 115 broken. Q4 holds: 306 of the 355 korkut changes are AKT readouts, and 324
of the 337 that have a PIP3 case repeat *PIP3 activates AKT signaling*'s own answer, which is wrong there (PIP3 alone
calls 655 korkut readouts DOWN, 515 measured unchanged). T adds nothing over A on the gated sets.

**These are not results either.** Every experimental case T1+F changes is PTEN in *Transcriptional Regulation by TP53*,
and against specs/132 L1 the overexpression cases moved to wrong (BAX, NDRG1, BTG2 DOWN; BIRC5 UP; measured the
reverse) while knockout read DOWN too. Traced (scratch Julia: every rule edge compared between orders, then the route):

1. Between `first` and `last`+amendment 1 the only differences are cross-pathway identity edges present under `first`:
   RNA twins joining TP53's and PIP3's PTEN mRNA both ways, recycled twins of GDP, GTP, H+, O2, cardiolipin and PE
   between pathways, carrier and component twins from TP53 into *HDR* copies. specs/132's "PTEN 9 / 0" leaned on
   those.
2. With them gone, the bridge sources, PIP3's p-AKT1/2/3 [cytosol], sit at the 100× ceiling under PTEN knockout AND
   overexpression. PIP3 solved alone reads PTEN overexpression DOWN (0.17×) and knockout UP (83.5×), correctly. The PIP3
   in the bundle was not PIP3: the builder (specs/056's design) removes from a source pathway every reaction of the
   pathway it feeds. TP53 transcription holds PTEN's transcription and translation, so the imported PIP3 lost exactly
   the PTEN machinery specs/095 and 100 fixed. (`--no-prune` changes nothing here: `main_curated` always pruned.)

Amendment 2 (`--whole`) enters each source pathway whole. `bench/analysis/source_check_133.jl` pins three starting
copies of each source pathway in turn and compares every one of its nodes with the pathway solved alone: in `.A2` the
worst |log ratio| over 49 (bundle, source) pairs is 2.2e-4 (solver tolerance). The same check on `.A1` (reduced and
pruned): 13 of the 49 pairs differ by more than 1% on nodes both hold, the worst |log ratio| 16 (*Signaling by
Insulin receptor* into *Signaling by Type 1 IGF1R*'s bundle), 4.6 in four, 1.1 for PIP3 in TP53's. The TP53 trace under `.A2`: PTEN
knockout puts p-AKT at 50–100×, TP53 tetramers fall to 0.28–0.53×, BAX 0.44× (measured DOWN); overexpression puts p-AKT
near 0, TP53 tetramers 1.15–1.32×, BAX 1.21× (measured UP). The signal now runs the curated route, PTEN ⊣ PIP3 → AKT
(PIP3 pathway) → bridge → AKT → MDM2 ⊣ TP53 → BAX (TP53 pathway). Arms `s133_A2`, `s133_T2` (code 5869e98).

## R5: A2 and T2 (both amendments), every changed case traced (2026-11-04)

Arms `s133_A2` (code 5869e98) and T2 (`s133_T2` korkut and decryptM at 5869e98; curator and experimental re-run as
`s133_T2b` at e67114c with `--max-edges 150000`, because T2's TP53 transcription bundle has 101,254 edges and the default
100,000 skipped it; `src/` and `bench/` are identical between the two commits; merged for scoring in
`analysis/133/s133_T2m`). Logs: order last; 62,000–74,000 cross-pathway rule edges dropped over all solves (whole source
pathways share many entities, so amendment 1 carries real weight). Against c0:

| | experimental fixed / broken / w→w | curator held-out | korkut | decryptM | pooled gain [95%] | gate |
|---|---|---|---|---|---|---|
| A2 rule | 10 / 14 / 9 | 13 / 18 | 97 / 395 | 2 / 33 | +0.0052 [−0.0053, +0.0161] | fails |
| **A2+F** | **6 / 0 / 5** | **12 / 4** | **71 / 250** | **0 / 22** | **+0.0147 [+0.0039, +0.0265]** | passes |
| T2 rule | 9 / 14 / 19 | 13 / 32 | 97 / 395 | 2 / 33 | +0.0046 [−0.0059, +0.0156] | fails |
| **T2+F** | **5 / 0 / 15** | **12 / 6** | **71 / 250** | **0 / 22** | **+0.0141 [+0.0037, +0.0254]** | passes |

Accuracy / macro-F1, rule+F minus c0: A2 experimental +0.0075 / +0.0069, curator +0.0004 / +0.0006, korkut −0.0162 /
+0.0008, decryptM −0.0122 / −0.0065; T2 experimental +0.0063 / +0.0072, curator +0.0003 / +0.0005, korkut and decryptM as
A2. Korkut once per (drug set, readout): 162 changed, 30 fixed / 115 broken / 17 wrong→wrong (both arms). Without F the
combination breaks 14 experimental cases the pathway alone routed: F is doing real work.

### The traced groups (`trace_132.py`; private `analysis/133/trace_A2F.tsv`, `trace_T2F.tsv`)

1. **PTEN → PIP3 → p-AKT → TP53 (experimental, TP53 transcription): done right.** A2+F: PTEN knockout BAX, PUMA, BTG2
   DOWN and BIRC5, CDC25C UP; overexpression BAX UP; all as measured (6 fixed). The route runs PIP3's own answer
   through the curated handoff (AKT phosphorylates MDM2). The 5 wrong→wrong are NOXA (PMAIP1), CASP6, TP53I3, PLK2 and
   the CDKN1A mRNA, which the experiment measured UP under PTEN knockout while BAX, PUMA and BTG2 fell in the same
   experiment; no one curated route gives both, and the model gives the curated one. A ground-truth question, not a rule.
2. **TP53 → BAX/BID (curator, Intrinsic Pathway for Apoptosis): done right.** 12 fixed. The 2 "broken" are Activated BAX,
   which the curators call NORMAL under TP53 knockout and overexpression while calling its downstream (cytochrome c
   release, apoptosome) DOWN and UP: their own answers disagree.
3. **AKT into RUNX2 (curator): right in construction, judged against a single pathway.** SRC (and, in T2, MAPK1) reach
   MMP13 through PI3K/AKT in PIP3 and the curated AKT → RUNX2 phosphorylation; the curators, judging RUNX2 alone,
   said NORMAL. 2 broken (A2), 4 (T2).
4. **p-AKT handoff into TP53, Mitotic G1, RUNX2, ERBB2, S Phase (korkut): right in construction, the answer is PIP3's
   and PIP3 is wrong.** 306 of 355 changes are AKT readouts; 324 of the 337 with a PIP3 case repeat PIP3's own answer
   exactly. PIP3 alone calls 655 korkut readouts DOWN, 515 of them measured unchanged. The error is in one pathway.
5. **VEGF's PIP3 into Rho GTPases (decryptM, 19 of 22): right in construction, the cell is not signalling.** PI3K
   inhibitors and dasatinib cut VEGFR2:PI3K's PIP3 that feeds the WAVE complex, and every downstream site reads DOWN;
   measured unchanged. decryptM's cells are not stimulated with VEGF; the model treats every pathway as active at rest.
6. **Depth 2 imports PI Metabolism's rails (T2 only): right in construction, the answer is PI Metabolism's and it is
   wrong.** *PI Metabolism* hands PI4P and PI5P to PIP3's PIP2 synthesis. Solved alone, it puts PI4P at 62× under PTEN
   knockout and 23× under overexpression (PI5P 99× and 2.4×): its lipid interconversion loops rail either way
   (specs/107 rule 1, not adopted). Through PIP3 that flips PTEN overexpression back to wrong (BAX DOWN) and adds PML's
   cases inverted (PML → PTEN regulation → PI Metabolism → PIP3). Every case T2 changes beyond A2 is this.

No group is wrong in construction any more. The two defects that were (identity twins across pathways; source pathways
that were not themselves) are amendments 1 and 2, each with a check that passes.

### Predictions

- **Q1** (ATM back under A): holds (R3).
- **Q2** (IL-4/13 losses gone): holds (R3).
- **Q3** (experimental ceiling rises by less than half over L1): **refuted** for T, 15 → 24 (+60%); holds for A (15).
  None of T's nine extra reachable cases is fixed (group 6).
- **Q4** (most korkut changes are AKT readouts repeating PIP3's answer): holds (group 4).
- **Q5** (the gate fails): **refuted.** A2+F and T2+F pass on pooled balanced accuracy. The pass rests on korkut, where
  raw accuracy falls (−0.016) and by unique (drug set, readout) 115 changes are broken to 30 fixed; balanced accuracy
  rewards the few correct DOWN calls more than it charges the many wrong ones. Experimental is a handful of PTEN cases.
  And the gate is not blind: these datasets' traces chose the fixes.

### What this says

The combination is now done right: each pathway answers as itself, pathways meet only at Reactome's curated handoffs,
and every changed case follows a curated chain. Done right, it mostly carries single-pathway errors further: PIP3's AKT
overcalls, PI Metabolism's lipid rails, pathways the cell is not running. Following the chain past one hop (T) adds only
group 6. The few gains it makes are real biology (PTEN through AKT to TP53's targets; TP53 to BAX/BID).

Not adopted; nothing changes a default. Adoption is Adam's decision. A confirmation would be one scoring on data these
traces have not touched (Lun 2019 or the CPPA test portions), which is also Adam's decision. The levers this points at
lie inside single pathways: PIP3's AKT readouts under drugs (korkut), and PI Metabolism's interconversion loops.

## Code review 2026-11-05 (builder and solver fixes; nothing re-run)

Defects found after R5, fixed in code and covered by tests. The A/A1/A2 and T/T1/T2 bundles were built and solved
before them, so they may have **missed** bridges (1, 2) or solved a carrier-loop case differently (5); none of the
arms is re-run here, and their numbers stand as recorded, not as the corrected rule's.

1. The builders counted a copy as a root only with no in-edge of any type; the bench pins depletion- and recycle-fed
   starting copies as roots. Both curated planners now use the bench's `root_indegree`.
2. A set-valued catalyst or regulator (a set pool node) was never bridged; a set handoff now enters the pool's members,
   member to same member. A pair with copies but no same member is `no_same_member_pair`, not `bridged`.
3. `--transitive` ignored `--l2`; the one-hop planner ignored the `extracellular` flag. Both now honour both, and
   BUILD.json records `context_l2`, `context_extracellular` and `context_root_ignore`.
4. `closure_133.py` left a tie between source pathways at equal depth to Neo4j's row order; it now takes the smaller
   stId and prints ties.
5. Solver under `last`: another pathway's starting copy suppressed a closed carrier's specs/084 component edges and the
   cross edge drawn instead was dropped; twin decisions now read the form's own pathway only. Rule edge counts are
   reported after the cross-pathway drop; a bundle with bridges and no origin table reports `last: no origin table`
   and the bench refuses it; influence scores resolve the pair rules without bridges, as the solve does.
