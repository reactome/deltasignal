# specs/113: Reactome interaction edges (census before any rule)

**Status:** analysis only, 2026-10-27. No rule is registered and nothing has been measured.

Adam, 2026-10-26: use the database's interaction records to add connections. Reactome 97 holds 78,863
`UndirectedInteraction` nodes (IntAct), each joining two or more `ReferenceGeneProduct`s through `interactor`. They
are unsigned and undirected, and they attach to proteins, not to the curated forms (EWAS, complexes) that the networks
are built from. HPN-DREAM is FIT after specs/112, so any rule would need a fresh test. The candidate test is LINCS P100,
which is TEST and still unseen. This census reads only P100's analyte table (probe → UniProt, site) and the graph.

## LINCS P100 readouts and what interaction edges could reach

| | count |
|---|---|
| probes / distinct probe proteins | 96 / 86 |
| probes whose phospho-site is a phosphorylated EWAS in Reactome 97 (specs/111) | 7 |
| probes whose site occurs in the 92 benchmark networks (specs/111) | 4 |
| probe proteins present in some form in the networks (`leaf_refs.csv`, frozen build) | 38 |
| probe proteins absent from the networks | 48 |
| … of which a `ReferenceGeneProduct` in Reactome | 43 |
| … of which have any EWAS in Reactome | 17 |
| … of which have at least one IntAct partner in the networks | 39 |

Interactions touching the networks:
- 47,942 of the 78,863 interactions touch at least one of the 3,956 network proteins.
- 14,812 join two network proteins.

## Reading

- **The readout is the bottleneck, not connectivity.** A P100 probe is one phospho-site. Interaction edges can make a
  protein reachable, but they cannot create a curated phospho-form of it. Only 7 of the 96 sites exist in Reactome at
  all, whatever is added.
- **What edges could reach.** They could reach most of the 48 missing proteins (39 have a partner in the networks).
  Scoring those proteins would need a readout at the protein level ("any form of the protein moves"), not at the
  site. That is a change of protocol, decided before any value is read.
- **No direction or sign.** The edges are unsigned and undirected. They could support a reachability-style score
  (does a perturbation reach the protein), not a signed prediction. specs/112 found reachability does most of the
  work anyway, so this is not useless. But it would measure the edges, not the solver.
- **Scale.** Adding 14,812 protein–protein links among network proteins (47,942 if partners outside are imported)
  to networks of 3,956 proteins would couple most of each pathway to most of the rest. Added edges have cost the curator
  axis every time they were tried (specs/016, 019, 053, 056, 085), and over-coupling is the main error class.

## The drug → target map (amendment 1 build, before any value was read)

Built by `bench/p100_targets.py` from ChEMBL 37. Responses are cached under
`~/deltasignal-catalogs/external/lincs_p100/chembl37_cache/`. The map is in `targets_20261027/` beside it, not in
the repo (CC BY-SA). SHA-256 of `p100_targets_1uM.tsv`:
`a804079ab4b5f6fea764758311813feebef9afdac6502f9f3d4204589db03b67`.

| source (first that gives a target) | compounds |
|---|---|
| kinobeads (R4, Kd ≤ 1 µM) | 17 |
| ChEMBL mechanism | 20 |
| ChEMBL binding, pChEMBL ≥ 6 | 43 |
| none (not scored) | 9 |

- **Identification.** 69 compounds by exact InChIKey, 13 by connectivity block, 4 by name. 764 compound–target
  rows; median 3 targets per compound, maximum 374 (a pan-kinase binder).
- **No target.** BMS-906024, GSK-J4, KN-93, SMER-3, bafilomycin A1, ginkgetin, methylstat, salermide and zebularine.
- **Sensitivity maps.** 100 nM: 722 rows. 10 µM: 896 rows. Both have the same 9 compounds without a target. Neither
  is chosen on outcome.
- **Interactions.** `intact_r97/interactions.csv` is the Reactome 97 export (interaction, accession, taxon), with its
  SHA-256 in `SHA256SUMS` there.

## Predictions (before any value was read)

`p100_bench.py predict` at 2ebbf29 used the frozen catalog `20261026-1349_8b6527e`. Output is
`results/2ebbf29/p100/p100_predictions.tsv` under that build (SHA-256 `a265b29263420a57e2ae0167386d17eb3116bbcf7264a9ba35872fb068f21ce2`). 80 compounds × 96 probes = 7,680 pairs.
74 of the 80 compounds have a target in the networks.

| arm | score 0 | score 1 | score 2 | compounds with any non-zero probe |
|---|---|---|---|---|
| A, curated reach | 7,512 | 168 | – | 58 |
| B, plus one interaction at each end | 2,041 | 5,471 | 168 | 77 |
| C, interactions alone | 4,871 | 2,643 | 166 | 76 |

- **Arm B.** It marks 73% of pairs as reached, as research.md expected: one unsigned hop at each end reaches most
  proteins. It can only help where the 27% it leaves at 0 are mostly unchanged probes.
- **Arm A.** It is non-zero for 2% of pairs. Within a context it can only separate a few probes from the rest.

## Scoring

- **First attempt (9e155ee).** Stopped by the instance guard before anything was computed. The instance table lists
  1,713 samples; the Level 4 matrix holds 1,684, and every one of them is listed. The guard required the two counts
  to be equal. It now requires only that every matrix column is a listed instance. The 29 extra listings were left
  out upstream; they are most likely QC removals (not checked).
