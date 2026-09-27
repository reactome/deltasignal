# specs/034 — How catalysts and regulators are wired: inventory

**Status:** analysis, a design input. Nothing is proposed for adoption here.

**Why (Adam, 2026-09-26):** "we might need to put a lot more thought to how
these regulators are connected ... there are a lot of potential
configurations." Four findings in a row traced to regulator wiring:
- set members wired as separate required inputs (specs/033);
- an arbitrary copy chosen for joins (specs/030, STAT1);
- drug complexes as inhibitors (specs/032);
- recycled catalysts closing loops (RAF).

Before designing, this inventory counts which configurations occur and which
cases route through them.

**Tool:** `bench/analysis/regulator_inventory.py --catalog BUILD --cases ARM`.
- Build `20260926-1221_590301c` (canonical); cases from arm
  `results/039ce74/poolbase`.
- The role comes from the Regulation subclass in Neo4j.

## 1. Annotations (Neo4j, 92 pathways)

| role | protein | complex | complex with set | set | set of complexes | small molecule | total |
|---|---|---|---|---|---|---|---|
| catalyst | 263 | 406 | 427 | 156 | 108 | 0 | 1,361 |
| positive gene-expression regulation | 42 | 141 | 120 | 4 | 4 | 1 | 312 |
| negative regulation | 46 | 59 | 108 | 13 | 9 | 2 | 237 |
| positive regulation | 28 | 50 | 38 | 15 | 4 | 12 | 147 |
| negative gene-expression regulation | 4 | 53 | 64 | 1 | 3 | 0 | 125 |
| requirement | 3 | 1 | 10 | 1 | 1 | 9 | 25 |

**The generator has three roles** (`_CATALYST_CYPHER`, `_POS_REG_CYPHER`,
`_NEG_REG_CYPHER`), and Neo4j's labels fold the rest into them:
`Requirement` and `PositiveGeneExpressionRegulation` both carry
`PositiveRegulation`, and `NegativeGeneExpressionRegulation` carries
`NegativeRegulation`. So a Requirement ("needed, not sufficient") is wired
exactly like a positive regulator: an AND input whose rise multiplies the
reaction.

## 2. Edges in the network, by role and wired copy

`copy` is the SOURCE node of the edge:
- **root:** no incoming edge;
- **upstream:** produced, and not reached from the reaction;
- **downstream:** produced by something the reaction itself reaches, so the edge
  closes a loop.

| role | edges | root | upstream | downstream |
|---|---|---|---|---|
| catalyst | 26,457 | 20% | 48% | **32%** (8,481) |
| negative | 1,954 | 16% | 82% | 2% |
| positive | 1,948 | 16% | 83% | 1% |
| positive gene expression | 409 | 8% | 88% | 5% |
| negative gene expression | 158 | 4% | 95% | 1% |
| requirement | 118 | 29% | 42% | 29% |

**A third of all catalyst edges are fed by a copy the reaction regenerates
downstream.** That is catalytic recycling (the loop taxonomy of 2026-06-12),
and it is rare for regulators. By shape:

| catalyst shape | edges from a downstream copy |
|---|---|
| set | 6,575 |
| complex with set | 1,435 |
| set of complexes | 369 |
| complex | 89 |
| protein | 13 |

The cells with the most edges are catalyst / complex with set / upstream
(7,210), catalyst / set / downstream (6,575) and catalyst / set / root (3,124).

## 3. Cases routed through each configuration

A case counts toward a cell if its pinned nodes reach the edge's source and the
edge's reaction reaches its readout. **A case can count in many cells, so these
accuracies are associations, not causes.** They point at where to trace.

**Low cells (n ≥ 40):**

| axis | configuration | n | accuracy |
|---|---|---|---|
| curator | positive / complex with set / downstream | 240 | **42.9%** |
| curator | catalyst / set of complexes / root | 30 | 40.0% |
| curator | negative / protein / downstream | 42 | 35.7% |
| curator | catalyst / set of complexes / downstream | 316 | 55.7% |
| curator | positive gene expression / complex with set / downstream | 280 | 58.2% |
| curator | requirement / complex with set / upstream | 46 | 63.0% |
| curator | negative / complex / upstream | 312 | 65.1% |
| curator | negative / complex with set / upstream | 948 | 67.8% |
| experimental | catalyst / set of complexes / downstream | 43 | **20.9%** |
| experimental | negative gene expression / complex with set / upstream | 56 | 37.5% |
| experimental | negative / complex / upstream | 62 | 46.8% |
| experimental | negative / complex with set / upstream | 87 | 48.3% |
| experimental | positive gene expression / complex with set / upstream | 97 | 56.7% |

**High cells**, for contrast: catalyst / set / upstream (curator 91.3%, n = 738),
catalyst / set of complexes / upstream (experimental 96.4%, n = 165),
negative / set / root (curator 100%, n = 60).

## 4. Mismatches between wiring and Reactome's meaning

| # | mismatch | size | status |
|---|---|---|---|
| M1 | Requirement wired as a multiplying positive regulator | 25 annotations, 118 edges | not tested |
| M2 | set-valued catalysts and regulators as separate AND/OR member terms | 1,881 + 56 + 146 groups | specs/033: pooling is neutral-to-positive; draft LNG #99 |
| M3 | catalyst fed by a downstream (recycled) copy | 8,481 edges, 32% of catalysts | `DS_SCC_BREAK_CATALYST` (TP53 +104, held-out −65) and specs/018 roles: not adopted |
| M4 | inputs fed by a downstream copy (RAF F-actin, CNKSR2, Ca2+) | not yet counted | specs/033 result: traced, open |
| M5 | drug-derived regulators propagate | 202 entities, 24 pathways | specs/032: inert not adopted |
| M6 | negative regulators that are complexes (with or without sets) | the low cells in §3 | not traced |
| M7 | gene-expression regulation wired like reaction regulation | 437 annotations | semantics plausible; low experimental cells, not traced |

**Next:** trace one case from each low cell in §3 (M6 and M7 first, because they
are the least explained), and count M4. Then propose a rule set for review.
