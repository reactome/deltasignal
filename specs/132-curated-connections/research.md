# specs/132 research

## R1. The curated handoffs (2026-11-03, no case data read)

`bench/analysis/curated_links_132.py` (Neo4j, Release 97) lists, for each of the 92 catalog pathways P, every reaction
r inside P whose curated `precedingEvent` s lies outside P, keeping the pair only where an entity (same stId) is an
output of s and an input, catalyst or regulator of r. Output: `~/deltasignal-catalogs/analysis/132/links.tsv`.

- **463 handoff pairs** (566 pair-entity rows) into **64 of the 92 pathways**. Roles: input 359, catalyst 190, positive
  regulator 14, negative regulator 3. 180 rows (32%) also pass L2 (s lies in a pathway P's diagrams encapsulate).
- Most rows: DSB repair-adjacent and TP53 pathways (R-HSA-75153 75, R-HSA-3700989 53), Rho GTPase/BCR, PIP3 (29),
  S phase (28).
- The entities handed off most are **active signalling forms**: TP53 tetramer, p-S1981,Ac-K3016-ATM, p-T,p-S-AKT,
  PI(3,4,5)P3, RAC1:GTP, p-T,Y MAPK dimers, CALM1:4xCa2+, cytosolic Ca2+ and IP3, G-beta:gamma, active IKK, active
  caspase-8, p-STAT3, p-STAT5. None of the most frequent is a free ligand, the specs/056 loss class (P2 is checked
  on the bridges themselves in R2).
- Each s's source network Q is its lowest ancestor pathway with a diagram (`s_q.tsv`): **124 pathways**, 40 already
  in the frozen build. The other 84 are built with the frozen generator (8b6527e) into the variant build
  `20261103-1700_8b6527e_ctx132` (the 92 plus the 84; a first attempt at 15:47 was stopped at Adam's request and
  deleted unfinished). All 176 built; no node id is shared between pathways, so a bundle cannot merge nodes by id.

## R2. The bridges (census, 2026-11-03, no case data read)

`bench/context_build.py --curated links.tsv s_q.tsv --members members.tsv --centres bench/catalog_pathways.tsv
[--l2]`, base `20261103-1700_8b6527e_ctx132`; per-row census in each bundle's `curated_census.tsv`.

Two matcher bugs were found by tracing the first census, before any arm ran, and fixed with tests (dfbe6c1):
- a catalyst or regulator edge points at the **reaction node** and carries no `edge_reaction_id` (only input and
  output edges do), so the first census missed every catalyst handoff (190 rows read "entity not at r");
- a handed-off **set** has no node of its own once the generator expands it into members (13 entities, 68 rows).
  `bench/analysis/set_members_132.py` lists their members from Neo4j (hasMember / hasCandidate, up to 3 levels), and
  a set handoff bridges each member to the same member.

| | L1 | L2 |
|---|---|---|
| handoff rows | 566 | 180 |
| bridged | 117 | 59 |
| P's copy not a root | 447 | 119 |
| entity not at r in P | 1 | 1 |
| held (cofactor, drug, pathogen) | 1 | 1 |
| pathways with a bridge | 27 | 12 |
| bridge edges | 74 | 31 |

**Why most handoffs make no bridge.** The registered rule bridges only into a copy that is a root in P, by the
bench's own root definition (in-edges from `logic_network.csv` and `boundary_edges.csv`, depletion ignored). Of the
447: about 150 are copies P produces itself (an `output` in-edge), so P already computes them; about 214 are complexes
P takes from outside but the generator rebuilds from their components (specs/030 boundary decomposition), so they have
`assembly` in-edges and are not roots. Bridging those would need a combination rule: under the default
`max(AND, OR)` a bridge OR-joined to the assembly would pass Q's increases and hide its decreases. That is a
different rule, not registered here; it is recorded as the main limit of L1 as registered.

**P2 (few free ligands or released forms).** Bridged L1 rows by entity: p-T,p-S-AKT as catalyst 25, p-S1981,Ac-ATM as
catalyst 21, p-CHEK2 4, p-CHEK1 3, mismatch-repair proteins (MSH2, MLH1, PMS2) 8, p14ARF 4, PI4P/PI5P 6, then single
rows (p-PLK1, p-IRS1/2, BID(1-195), CDC25C, MDM2, ATF3, ...). Free extracellular ligands: 2 rows (IL4 and IL21 into
the IL-4/13 and IL-2 family pathways, from *Affinity selection of immunoglobulins*). Small molecules: H2O2, GSSG,
reduced cytochrome c, PA, L-Gln (into TP53 transcription, from metabolic pathways). **P2 is supported**: the bridges
are mostly activated kinases handed to their substrates' pathways; the specs/056 loss class (free EGF, KDM1A) is
almost absent.

Pathways receiving most: Transcriptional Regulation by TP53 (40 rows, 14 from DNA Damage/Telomere Stress Induced
Senescence, 8 from PIP3), Cell Cycle Checkpoints 12, Intrinsic Pathway for Apoptosis 9, PIP3 8, Mismatch Repair 4,
Mitotic G1 4 and S Phase 5 (p-AKT from PIP3), Signaling by ERBB2 5.

Bundles: `.c0` (empty link file: each pathway alone, the control), `.l1`, `.l2`. Arms `s132_c0`, `s132_l1`,
`s132_l2` at solver dfbe6c1 on curator, experimental, korkut and decryptM (hand target map,
`analysis/114/cases_frozen92_hand.tsv`). The korkut and decryptM benches now read only P's own nodes in a bundle
(0ed9ad0), as the curator bench has since specs/056; before that fix a bundle would have scored readouts on its
neighbours' nodes.

**Deviation from the registered order.** Step 1 (ceiling count, no solving) needs each case's pins and readouts
resolved in the bundle, which the benches do. The count is therefore read from the arms' own case files (pins,
readouts and the bundle's reachability; the bench's `no_path` and korkut/decryptM reach columns), and it is computed
and written down before any accuracy figure is read. The stop rule only stops; it chooses nothing, so reading it from
the same runs cannot bias a gate.

## R3. The arms (2026-11-03): no rule passes

Arms at solver 0102aaa (bundles rebuilt with `rnas.csv`, which the first bundles lacked; the bench refused them under
`DS_RNA_TWIN=on`), scored by `bench/analysis/combine_132.py` (dcabf4f, tested). Results under
`builds/20261103-1700_8b6527e_ctx132.{c0,l1,l2}/results/0102aaa/s132_{c0,l1,l2}`; scores in
`~/deltasignal-catalogs/analysis/132/score_{l1,l2}.json`.

**Control.** `s132_c0` reproduces the frozen model (`cf_off`, d210ae3; decryptM `s114_f92_hand`) on every case:
experimental 849, curator 24,100, korkut 11,626, decryptM 1,862, 0 differ. Rebuilding the 92 networks in the
variant build changes nothing, so every difference below is the bridges.

**Ceiling (read first).** Control's wrong cases with no route in P alone but a route through the rule:

| | experimental (RAF aside) | korkut | decryptM | curator held-out |
|---|---|---|---|---|
| L1 | 15 / 194 (7.7%) | 106 / 3,945 (2.7%) | 0 / 671 | 12 / 1,703 (0.7%) |
| L2 | 15 / 194 (7.7%) | 81 / 3,945 (2.1%) | 0 / 671 | 12 / 1,703 (0.7%) |

Above 2% on experimental and korkut, so the stop rule does not apply; the arms are read.

**Candidate gate** (pooled three-class balanced accuracy, experimental RAF aside + korkut, arm − control, 95% interval
resampling perturbations within dataset, 10,000, seed 132):

| arm | pooled gain | 95% interval | experimental | korkut | gate |
|---|---|---|---|---|---|
| L1 | +0.0023 | [−0.0107, +0.0159] | +0.0057 (14 fixed / 8 broken) | +0.0084 (72 / 534) | fails |
| L1+F | +0.0114 | [−0.0011, +0.0252] | +0.0047 (9 / 4) | +0.0174 (69 / 322) | fails |
| L2 | +0.0001 | [−0.0113, +0.0122] | +0.0057 (14 / 8) | +0.0040 (53 / 435) | fails |
| L2+F | +0.0084 | [−0.0018, +0.0200] | +0.0047 (9 / 4) | +0.0121 (50 / 249) | fails |

Reported, not gated: decryptM balanced −0.009 (L1: 2 fixed / 42 broken), −0.005 under F; curator held-out L1 −80
(14 / 94), L1+F −58 (12 / 70), L2 and L2+F +6 (6 / 0). **No rule passes, so nothing goes to stage 5.**

**Traced: experimental (L1 = L2: every experimental change comes through an L2 handoff).**
- Gains: PTEN knockout and over-expression in *Transcriptional Regulation by TP53* (9 fixed, 5 wrong both ways),
  through the p-AKT handoff from *PIP3 activates AKT signaling*; TP53's own network has no route from PTEN. AKT1/AKT2
  knockouts and BRCA1 knockout there gain 5 more (they had a route, but read NORMAL alone).
- Losses: ATM knockout and over-expression in *Cell Cycle Checkpoints* (4). The control was right without a route
  that reachability sees: the route is specs/084's solver-only component-twin edges from ATM's starting copies. A
  bridge into p-ATM's root copy gives it an in-edge, it stops being a starting copy, and the twin edges go. **Bridging
  changes which nodes are roots, and every root-based solver rule (twins, pins) moves with it.** The same blind spot
  makes F read these cases as "no route" and take the arm's answer. SRC over-expression in *Mitotic G1* (3): SRC is
  pinned in the neighbour too; F keeps the control there.

**Traced: curator held-out (guard).** 90 of L1's 94 losses are in *Interleukin-4 and Interleukin-13 signaling*
(IRF4 38, STAT6 32, STAT3 20), all through its one bridge: IL4 handed from *Affinity selection of immunoglobulins*,
where IRF4, STAT6 and STAT3 drive IL4 production, so every IL-4 readout moves where the curators said NORMAL. The
registered sample of 20 broken cases (seed 132) is all IL-4/13; all 20 are carried by the curated connection, as every
change in a feed-in bundle is by construction, so that classification separates nothing. The useful finding is what
the connection is: **IL4 is secreted by T follicular helper cells and acts on B cells: a curated precedingEvent
between two cell types**, the free-ligand class P2 warned of. A one-cell model wires one cell's output to another
cell's input through it. The other 4: TP53 in intrinsic apoptosis (2), SRC in RUNX2 (2).

**Traced: korkut.** 687 changes; 522 of the broken are NORMAL → a change call on readouts measured "not", mostly in
*Transcriptional Regulation by TP53* (251) and *Mitotic G1* (129); 71 correct "down" calls are gained. Balanced
accuracy rises a little while raw accuracy falls by 462: the trade specs/128–131 kept finding, here in the other
direction (more change calls).

**Predictions.**
- P1 (L1 reaches more wrong cases per bridge than specs/056): R4.
- P2 (few free ligands among the bridges): **supported by count, but the one free-ligand bridge (IL4) causes 90 of 94
  curator losses.** Few is not enough; one cross-cell ligand is enough to do the damage.
- P3 (gains concentrate in cell-cycle pathways): **not supported.** The gains are in TP53 transcription (PTEN, AKT);
  the cell-cycle pathways lose (ATM in checkpoints, SRC in Mitotic G1).
- P4 (F matters more on korkut than experimental): **supported.** F doubles korkut's gain (+0.008 → +0.017) and
  barely moves experimental (+0.006 → +0.005).
