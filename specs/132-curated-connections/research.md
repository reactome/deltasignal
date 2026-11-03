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
  `20261103-1547_8b6527e_ctx132` (the 92 plus the 84).
