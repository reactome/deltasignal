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
