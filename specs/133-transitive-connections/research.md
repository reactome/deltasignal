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
