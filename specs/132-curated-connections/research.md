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
