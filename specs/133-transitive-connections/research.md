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
