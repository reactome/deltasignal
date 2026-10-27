# specs/114 research

## Part 1: the build (2026-10-27)

`scripts/catalog.sh build --variant top114 --list specs/114-whole-reactome/top_level_pathways.tsv`, generator `8b6527e`
(`freeze-2026-10-26`), Reactome 97. Build `20261027-1346_8b6527e_top114`, 4,232 s.

- **Complete:** 28 of 28 pathways built; none missing, empty or truncated (`BUILD.json`). `current` not moved.
- **Size:** 826,677 nodes and 1,756,544 edges summed over the 28 networks (the frozen 92-pathway catalog: 209,480 nodes). The
  largest: Signal Transduction 109,395 nodes, Cell Cycle 100,117, Immune System 90,072, Cellular responses to stimuli
  87,872, Metabolism of proteins 82,544, Metabolism 68,565, Gene expression (Transcription) 61,363. The largest network of the frozen
  catalog has 25,692 (Signaling by Rho GTPases).
- **Generator warnings:** 537 variant-cap hits (a set combination over the cap bundled into one node, as in every
  build since specs/036) and 849 sets that resolved to no leaves. Neither stops a build; neither was counted per pathway.

## Part 1: coverage (graph only)

`bench/analysis/coverage_114.py` (tested in `bench/test_coverage_114.py`). No measured value is read: decryptM's
`all_sites.tsv` is read for its first two columns only. Definitions as in specs/113: a readout's protein is *in* if a
node's `member_leaves` holds a leaf whose `leaf_refs.csv` reference is its UniProt accession; its site is *in* if a node
holds a human EWAS of that protein phosphorylated (S/T/Y) at that coordinate in Release 97. Gene symbols (korkut) go
through every human `ReferenceGeneProduct`'s gene names (`gene_accession_r97.csv`, beside the phospho export, outside
the repo).

| dataset | readouts | build | protein in | site curated in Reactome 97 | site in |
|---|---|---|---|---|---|
| LINCS P100 | 96 probes | 92 pathways | 42 | 7 | 4 |
| | | 28 top-level | 61 | 7 | **7** |
| korkut2015 | 73 antibodies | 92 pathways | 68 | 29 | 22 |
| | | 28 top-level | **73** | 29 | **29** |
| decryptM | 73,369 measured sites | 92 pathways | 25,298 | 532 | 377 |
| | | 28 top-level | 45,946 | 532 | 506 |

The P100 row reproduces specs/113 on the 92-pathway catalog (4 of 96 sites in). korkut's 73 include 44 total-protein
readouts, which have no site. decryptM's row is every site measured, not the scored unit (a site under a drug whose
target is in the same network); the scored-unit counts come with part 2's join.

**Reading.** The whole of Reactome holds every curated readout site of P100 and korkut, and 506 of decryptM's 532.
Pathway choice no longer limits site coverage; what Reactome curates does. For P100 that ceiling is 7 of 96 sites, so
the frozen model still cannot be scored there in any useful way. The 26 decryptM sites curated but in no network were
not traced (candidates: the Disease top-level pathway, left out, and EWAS in no reaction).
