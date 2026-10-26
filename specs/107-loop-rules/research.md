# specs/107 results

Code: `src/core/material_loops.jl` (rules 1–3), `leaf_refs` on `ReactionNetwork` and
`parse_leaf_refs` (`src/io/tsv_parser.jl`), applied first in `solve_steady_state` and in
`solved_carrier_network`. Tests: `test/test_material_loops.jl` (44). Generator: `leaf_refs.csv`
(logic-network-generator PR #128, branch `leaf-refs-107`).

Catalogs: the reference table needs the generator change, so both catalogs were copied as overlays
with `leaf_refs.csv` added and nothing else changed: `20261024-2211_7f4b28a_refs107` (benchmark) and
`20261025-1318_7f4b28a_sample101_refs107` (sample). Each overlay's `BUILD.json` records its source build
and how the file was written.

## Implementation notes found while building

- A protein node is a `simple_entity` too; rule 1's "small molecule" therefore means every leaf
  is an `R-ALL-` entity (caught by the first test fixture).
- The specs/106 audit missed cofactor holds; that is specs/106 amendment 1, applied here.

## B1 coverage (specs/106 audit, `loop_coverage.jl` with the rules applied)

Sample (101 pathways, 361 loops):

| arm | dissolved | partly handled | untouched | loop nodes in dissolved loops | material-closed remnants |
|---|---|---|---|---|---|
| rules off | 299 | 20 | 42 | 53% | 20 loops / 15 pathways |
| rule 1 | 303 | 20 | 38 | 53% | 16 / 13 |
| rules 1 + 2 (`strict`) | 304 | 20 | 37 | 53% | 15 / 12 |
| rules 1 + 3 (`wide`) | 334 | 14 | 13 | 83% | 11 / 9 |

Benchmark (92 pathways, 744 loops):

| arm | dissolved | partly handled | untouched | material-closed remnants |
|---|---|---|---|---|
| rules off | 657 | 33 | 54 | 21 / 7 |
| rules 1 + 2 | 660 | 33 | 51 | 21 / 7 |
| rules 1 + 3 | 670 | 27 | 47 | 15 / 4 |

## Predictions

- **P1 fails:** rule 1 dissolves 4 of the 11 sample metabolite cycles (predicted at least 8). The
  rest are two-substrate reactions, which a one-in/one-out conversion does not cover: transaminations
  (GOT2, OAT), transketolase, GPX2 (GSH with H2O2), pyruvate carboxylase (pyruvate with HCO3-), an
  N-acetyltransferase with acetyl-CoA. Cytochrome c oxidised ⇄ reduced is a protein redox pair the
  generator's pool finder does not list; rule 3 reaches it through the reference table.
- **P2 fails:** `strict` changes no benchmark material-closed remnant. Traced
  (`why_not.jl`, per entity node on each remaining cycle):
  - the carriers of PAK (Signaling by Rho GTPases) and L1/microtubule pass through set-variant
    nodes, which `strict` does not walk and `wide` does;
  - tBID:BCL-2 is listed by the generator as an open loop and kept on purpose (BAD displaces tBID,
    specs/077);
  - the DDX58/IFIH1 TRAF3:MAVS cycle is listed as an open scaffold loop and kept on purpose
    (specs/077);
  - cyclin B1:CDK1 between nucleus and cytoplasm changes CDK1's modification (different leaves),
    which only the reference table joins.
- **P3 fails narrowly:** rules 1–3 leave material-closed remnants in 9 of the 15 sample pathways
  (predicted at most 8).

## Behaviour checks B2–B5 (`behaviour.jl`, `score_behaviour.py`)

Perturbations and readouts were chosen once, with every rule on (`select`), and every arm solved the
same ones: the root nodes reaching a node the rules touch (at most 8 per pathway, by stable id), at
doses 0, 0.5, 2, 10, plus the unperturbed solve; readouts are the touched nodes and up to 10 terminal
nodes they reach. Rails: a readout ≥ 99× or ≤ 0.01× under a 0.5× or 2× root. Monotone: folds over
doses 0, 0.5, 1, 2, 10 never rise and fall. Relabel: one perturbation per pathway on a copy with every
uuid renamed and the edge order reversed.

Sample: 14 pathways touched, 99 roots, 410 solves per arm.

| arm | converged | rails | monotone | relabel, max fold difference |
|---|---|---|---|---|
| rules off | 100% | 394 / 2,332 (16.9%) | 1,132 / 1,166 (97.1%) | 7.8e-14 |
| rule 1 | 100% | 385 (16.5%) | 1,132 (97.1%) | 7.8e-14 |
| rules 1 + 2 | 100% | 385 (16.5%) | 1,132 (97.1%) | 7.8e-14 |
| rules 1 + 3 | 100% | 88 (3.8%) | 1,159 (99.4%) | 2.2e-14 |

No sample pathway worsens under any arm. Under rules 1 + 3, Macroautophagy goes from 80 rails to 0,
R-HSA-75893 from 84 to 0 and R-HSA-977225 from 18 to 0; Hedgehog 'off' state from 37 to 64 of 64 pairs monotone.

Benchmark (overlay of the canonical build): 13 pathways touched, 421 solves per arm.

| arm | converged | rails | monotone | relabel |
|---|---|---|---|---|
| rules off | 100% | 767 / 2,372 (32.3%) | 1,061 / 1,186 (89.5%) | 5.6e-14 |
| rule 1 | 100% | 761 (32.1%) | **1,059** (89.3%) | 3.3e-14 |
| rules 1 + 2 | 100% | 670 (28.2%) | **1,057** (89.1%) | 3.3e-14 |
| rules 1 + 3 | 100% | 227 (9.6%) | 1,154 (97.3%) | 1.2e-14 |

**Against the registered criterion** (B2 and B4 not worse, B3 not higher, B5 holds), on both catalogs:

- **Rules 1 and 2 fail on the benchmark (B4), narrowly.** Rule 1 loses two monotone pairs: one readout,
  PI(4,5)P2 in PIP3, under the PI4P and PI5P roots (the second listing of the same readout as a
  terminal node is not counted twice). Rule 2 also loses two in the Fanconi anemia pathway
  (R-HSA-6783310), where it cuts rails from 187 to 112.
- **Rule 3 (with rule 1) passes on both catalogs:** benchmark rails 32.3% → 9.6%, monotone 89.5% → 97.3%.
  Per pathway it is not uniform: cyclin B1 / G2/M (R-HSA-453274) rails 136 → 144 while its monotone pairs go 46 → 122;
  Class I MHC (R-HSA-983169) monotone 95 → 87; PIP3 as under rule 1.

**The PIP3 trace (rule 1).** Rule 1 makes PI(4,5)P2 ⇄ PI(3,4,5)P3 a pool (PI3K forward through 372
reaction copies, PTEN back). With PI4P at 2×, the pool's total rises to 1.33× (supply: PIP5K at 2×,
PIP4K at 1×), as conservation requires, but the split moves to PIP3 (2.40×) and PIP2 falls to 0.26×:
PI3K's catalysts lie downstream of PIP3, a catalyst-closed positive loop, so more PIP3 speeds PI3K and
drains PIP2. At 0.5× PI4P the loop slows and PIP2 rises to 1.17×. PIP3 itself is monotone (0.50,
1, 2.40); with the rule off it rails to 12.2× at 2× PI4P. The non-monotone readout is the substrate
of an enzyme under positive feedback, not a numerical fault, but the registered criterion counts it.

## No-harm guard and reported axes (benchmark arms)

Catalog `20261024-2211_7f4b28a_refs107` (canonical build plus `leaf_refs.csv`). Control `ctl107`
(solver 04b88f2, defaults) is identical to `dflt102` on the canonical build, row for row on both
axes. Arms at solver 6fe789a, each read against `ctl107` with
`holdout_report.py --pair-on-valid --set-aside RAF_MAP_kinase_cascade`. Each log's tally confirms its
rule ran (32 metabolite pool-solves in every arm; `strict` added 36–48 loop-row-solves per axis,
`wide` 484–3,224).

| arm | curator held-out (guard) | tuning | experimental | guard |
|---|---|---|---|---|
| `m107` rule 1 | 0 (0 / 0) | 0 | 0 | passes |
| `ms107` rules 1 + 2 | **+6 (6 / 0), p 0.031**, 1 of 70 pathways | 0 | 0 | passes |
| `mw107` rules 1 + 3 | **−50 (6 / 56), p < 1e-4**, 3 of 70 pathways | −35 (5 / 40) | +1 (3 / 2, p 1) | **fails** |

RAF changes nothing in any row (tuning with RAF set aside equals tuning).

**Traced (`wide_rows.jl`, the rows each mode adds per pathway, entities read in Neo4j):**

- **Rule 2's +6 is free ubiquitin.** Its only row in the Fanconi pathway is free Ub (R-HSA-68524,
  a set of UBB/UBC/UBA52/RPS27A, protein leaves, so `strict` accepts it) released by FANCD2
  deubiquitination by USP1:WDR48 (R-HSA-6786171). Closing it holds free ubiquitin, and the 6 USP1
  cases are exactly the specs/098 gain (`DS_UBIQUITIN_MODE=inert`: held-out +6, USP1). It is not a
  multi-step protein carrier gain.
- **Rule 3's TGF-beta losses (−48 of the curator total, SMAD2/3, TGFBR1/2) are a modification cycle
  read as a carrier.** The reference table matches SMAD2 [cytosol] with p-S465,S467-SMAD2 [cytosol]
  (and SMAD3 likewise), so phosphorylation in the receptor complex, release (R-HSA-170850), and
  MTMR4 dephosphorylation (R-HSA-2187401) form a closed "carrier" loop of four forms. The closed rule
  then holds p-SMAD2/3, the signal itself, at baseline: receptor knockouts and SMAD overexpression
  read NORMAL. A modification changes what a protein does; it is not a carrier released unchanged.
- **Rule 3's Mitotic G1 losses (SRC, JAK2, ABL1, PTK6 knockouts, 8 each) are free ubiquitin again,**
  here Ub [nucleoplasm] closed through proteasomal degradation of ubiquitinated p27/p21
  (R-HSA-187574), so the p27 degradation route these kinases act through no longer moves. This is
  the specs/098 loss class (CCNE1 there).

**Outcome against the registered gate:** no arm passes every criterion. Rule 1 passes the guard
(no class change) and fails B4 on the benchmark by two pairs (PIP2 under PI3K feedback, traced
above). Rule 2 passes the guard and fails B4 by two pairs (Fanconi); its accuracy effect is the
ubiquitin hold specs/098 already measured and Adam declined. Rule 3 passes B2–B5 on both
catalogs and fails the no-harm guard. Defaults stay off; adoption is Adam's.

**What the traces say about the rules, not the data.** Both rule 3 mechanisms come from the
reference match, not from the multi-step walk: matching across modifications admits
phosphorylation cycles, and free ubiquitin is a protein set that every ubiquitination releases.
A variant that matches references only across compartments (same modifications) and leaves
free ubiquitin to `DS_UBIQUITIN_MODE` would answer both, but it is chosen from curator held-out
traces, so it would need its own pre-registration and Adam's decision before it is run.

## Amendment 1: rule 3b (`DS_CARRIER_MULTISTEP=compartment`)

Generator: `leaf_refs.csv` gains `modifications` (LNG #129, branch `leaf-mods-107`). New overlays
`20261024-2211_7f4b28a_refs107m` and `20261025-1318_7f4b28a_sample101_refs107m` (hard-linked copies
plus the three-column table, written by LNG e9483e9; `BUILD.json` records it). Solver 6fadb58.

**Rows (A1, rows part):** rule 3b adds no row in TGF-beta (the SMAD2/3 rows are gone) or Mitotic G1
(free ubiquitin is gone). In Fanconi, the free-ubiquitin row is gone too, and FAN1 [nucleoplasm], released
unchanged after unhooking the crosslink, is listed as a carrier (it was listed under `wide` as well).
Cyclin B1:CDK1 (G2/M, R-HSA-453274) is no longer touched, as predicted.

**Behaviour checks** (the registered selections, rule 1 off):

| catalog | arm | converged | rails | monotone | relabel |
|---|---|---|---|---|---|
| sample | rules off | 100% | 394 / 2,332 (16.9%) | 1,132 / 1,166 | 7.8e-14 |
| sample | rule 3b | 100% | 350 (15.0%) | 1,152 | 2.2e-14 |
| benchmark | rules off | 100% | 767 / 2,372 (32.3%) | 1,061 / 1,186 | 5.6e-14 |
| benchmark | rule 3b | 100% | 471 (19.9%) | 1,074 | 5.6e-14 |

B2–B5 pass on both catalogs. Per pathway it is not uniform: Macroautophagy (R-HSA-1632852) monotone
240 → 233, Class I MHC (R-HSA-983169) 95 → 87 (as under `wide`).

**No-harm guard (`mc107` against `ctl107`): fails.** Held-out **−10 (0 / 10), p 0.002**, one pathway;
tuning 0; experimental 0.

All ten are knockouts in Class I MHC antigen processing: nine E3 ligases (CBLB, FBXO11, FBXW7, KEAP1,
RNF213, SOCS1, TRAF7, VHL, ZBTB16) at antigen peptide bound class I MHC [plasma membrane] (fold
0.667 → 0.917, DOWN → NORMAL), and HLA-A at the same complex [Golgi membrane] (0.0028 → 0.857). Curators
say DOWN. Rule 3b lists two carriers there: UBA7 (E1, released by "Transfer of ubiquitin from E1 to E2")
and the BAP31 oligomer (released by "Transport of MHC heterotrimer to ER exit site"). Adding them one
at a time (`mhc_iso.jl`) shows that **BAP31 alone moves every case** and UBA7 changes nothing. Rules
off, BAP31's loop iterates: a single ligase knockout drives Golgi MHC to 0.0015× and an HLA-A
knockout to 0.0028×. Those are rails (B3 counts them), and the curator DOWN was matched through them.
With BAP31 held as a closed carrier, the same knockouts still lower both readouts (0.75×, 0.92×,
0.86×), but the plasma-membrane values fall inside NORMAL. The direction agrees with the curators; the
class change comes from the rail's magnitude going away. BAP31 is the carrier this rule was written
for (an ER chaperone escorting MHC to exit sites and returning), so this loss is the rule working as
designed, not a new wrong join.

**Predictions:** A1 fails (rows as predicted, but the guard fails, at the edge of the predicted −10 to
+10 range); A2 holds on the benchmark (19.9% ≤ 20%), fails on the sample (15.0%, predicted ≤ 10%: most
of `wide`'s sample gain came through modified-form matching); A3 holds (0).

**Outcome:** rule 3b fails the registered no-harm guard, so it stays off. Whether the ten cases count
as justified disagreements (curator DOWN reached only through a railed loop) is Adam's call.

**Decision (Adam, 2026-10-26): not adopted.** Rule 3b improves every behaviour check but no accuracy
axis (held-out −10, tuning 0, experimental 0), and under specs/108 only one of the ten losses (HLA-A at
the Golgi readout) is justified. A "dose-insensitive response" artefact class, which might have excused the
nine ligase cases, is **not** registered: it would be defined after seeing the cases it rescues. Defaults
for rules 1–3b stay off. Work moves to positive loops closed by a catalyst (specs/106).

**Why the behaviour gain does not show in accuracy (analysis after the decision, 2026-10-26).** Rule 3b
removes 296 rails among the behaviour probes at the mild doses, 192 of them in NHEJ (R-HSA-5693571), 48 in
R-HSA-74752, and 18 in Class I MHC. Against `ctl107`, only 85 valid curator values move and no experimental value
moves. Of the scored curator cases it takes off a rail, 32 in all, **30 keep the same class and were right
before and after**. One goes from right to wrong (HLA-A at the Golgi readout), and one was wrong and stays
wrong. In NHEJ, all 32 valid curator cases are right and all 32 are railed under the control; 6 come off the
rail and stay right. The three-class score gives no credit when a railed answer pointing the right way becomes
a moderate one, so a rule whose main effect is removing rails can lose on these axes but hardly gain. The nine
ligase losses are the only class changes not on a rail: a readout near the 0.85 cutoff (0.667 → 0.917).
