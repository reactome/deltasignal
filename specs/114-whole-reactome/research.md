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

## Part 2 groundwork: can the frozen solver handle these networks? (2026-10-27)

A feasibility probe, not a score: the API at `src` of `6d2b71b` (solver unchanged since `freeze-2026-10-26`), dev
compose environment, serving the top-level build; one solve per network with an arbitrary root node pinned to 0. No
dataset was read.

| network | nodes | parse | solve | converged | iterations |
|---|---|---|---|---|---|
| Cell Cycle | 100,117 | 31 s | 34 s | yes (residual 4e-15) | 14 |
| Signal Transduction | 109,395 | 22 s | 65 s | yes (residual 9e-15) | 50 |

The container peaked at 2.3 GiB. The first solve includes Julia compilation. The largest unjoined network solves in
about a minute, so per-network solving is feasible; a joined network of all 28 (826,677 nodes before merging) was not
tried.

## Part 2: scores (2026-10-27)

Run with `scripts/run_arm.sh`, code defaults, `src/` unchanged since `freeze-2026-10-26`. Results under each build's
`results/<code>/`: frozen 92-pathway build `20261026-1349_8b6527e` at `a47898d` (`s114_f92_kb`, `s114_f92_hand`);
top-level build `20261027-1346_8b6527e_top114` at `a47898d` (`s114_top_kb`, korkut and decryptM kinobeads) and at
`a0a6a3d` (`s114_top_hand`; that commit only added `bench/analysis/summary_114.py` and its test, so the bench and
solver are the same). Figures over all cases come from `bench/analysis/summary_114.py`; split figures from the bench's
own `summary.json`. Every result below is a development result on FIT data.

### korkut2015

| build | scored cases | pathways | solves converged | model | reachability | always-NOT | regulated right (model / reach, of) | McNemar (model only / reach only) |
|---|---|---|---|---|---|---|---|---|
| 92 pathways | 11,062 | 54 | 1,493 / 1,493 | 0.643 | 0.566 | **0.908** | 450 / 499 of 1,013 | 908 / 56 |
| 28 top-level | 17,836 | 21 | 793 / 814 | 0.723 | 0.631 | **0.913** | 444 / 601 of 1,545 | 1,873 / 244 |

78 drug sets on both builds. Units (drug set, readout): 3,480 on the 92 pathways, 1,125 of them called differently by
different networks; 5,398 on the top-level build, 2,261 called differently (42%).

### decryptM

Kinobeads map (primary), held-out split (on the top-level build every case is held-out):

| build | scored | regulated: model / reach / always DOWN (of) | McNemar regulated (model only / reach only) | by-drug interval, model − reach | unchanged: model / reach (of) | macro-F1 model / reach |
|---|---|---|---|---|---|---|
| 92 pathways, held-out | 1,133 | 97 / 112 / 248 (278) | p 6.1e-05 | [−0.093, 0] | 632 / 592 (855) | 0.368 / 0.364 |
| 28 top-level | 5,096 | 187 / 301 / 606 (725) | 4 / 118, p 3e-30 | [−0.238, −0.100] | 3,316 / 2,834 (4,371) | 0.363 / 0.330 |

Hand map: on the 92 pathways, held-out regulated 54 / 67 of 190, unchanged 485 / 452 of 621, macro-F1 0.365 / 0.364.
On the top-level build, regulated 147 / 272 of 649 (interval [−0.276, −0.126]), unchanged 2,854 / 2,425 of 3,751, macro-F1 0.358 / 0.332.

All cases (both splits on the 92 pathways), kinobeads map: 2,334 → 5,096 scored. Accuracy is 0.635 → 0.687 for the
model, 0.596 → 0.615 for reachability and 0.778 → 0.858 for always-NOT. Units (drug, site): 1,374 → 3,127, of which
196 → 450 get different calls from different networks.

### Coverage count 3: reachable from a pinned target in the same network

A unit is counted if, in at least one network, the bench's reachability call is not NOT.

| dataset | build | units | reachable | readouts reachable |
|---|---|---|---|---|
| korkut | 92 pathways | 3,480 | 2,170 | 53 of 73 |
| | 28 top-level | 5,398 | 3,797 | 60 of 73 |
| decryptM (kinobeads) | 92 pathways | 1,374 | 701 | 197 sites |
| | 28 top-level | 3,127 | 1,621 | 348 sites |

### LINCS P100

Score A (curated reach AUROC over contexts) on the top-level build is **0.498** [0.493, 0.503]; on the 92 pathways it is
0.499. Predictions cover 11,429 network proteins, and 79 of the 80 compounds have a target in a network.

### Hypotheses

- **H1 holds.** More cases are scored on the top-level build: korkut 11,062 → 17,836; decryptM 2,334 → 5,096
  (kinobeads) and 1,802 → 4,400 (hand).
- **H2 fails.** On the regulated decryptM cases, the model gets 187 of 725 right and reachability 301, with an interval
  of [−0.238, −0.100] that excludes 0. The loss is larger than on the 92 pathways (97 vs 112 of 278).
- **H3 fails.** On korkut the model's accuracy (0.723) beats reachability's (0.631) but not always-NOT's (0.913).
- **H4 fails.** P100 score A is 0.498, and its interval includes 0.5.

**Reading.** The whole of Reactome adds cases but does not change the pattern from the 92 pathways. On cases measured
unchanged, the model beats reachability because it predicts NOT for most reachable sites. On the regulated sites the
model loses, by more than it did on the 92 pathways: on top-level decryptM, 118 regulated cases are right only under
reachability and 4 only under the model. Both benches stay below always-NOT. Larger networks also give more
conflicting calls: 42% of korkut units get different calls from different top-level networks. These networks are
larger and overlap, and 21 of 814 korkut solves did not converge, against none on the 92 pathways.

**Amendment 1 (2026-10-27, after the results): what the 118 decryptM cases are.** The first reading above said the
model calls NOT on reachable regulated sites. That was an inference from the totals, written before any case was
looked at, and it is mostly wrong. Of the 118 regulated cases right only under reachability, **73 are inversions**
(measured DOWN, model UP) and 45 are NOT calls. The inversions are mostly ERBB-family inhibitors (afatinib 19,
gefitinib 18, lapatinib 12) and SHP099 (9), 68 of 73 in Signal Transduction, and they include ERK1/2 activation-loop sites
(MAPK3 T202/Y204, MAPK1 T185/Y187). Over all 725 regulated cases: 391 are not reachable in the network they are scored in (305 measured
DOWN, 86 UP), so the model cannot get them there; reachability calls only DOWN, so it gets none of the 119
measured UP (corrected 2026-10-28 from 111: 86 unreachable + 33 reachable). No case was traced through a network here.

**Amendment 2 (2026-10-27, post hoc): balanced accuracy.** Plain accuracy is dominated by NOT: 91% of korkut and
86% of top-level decryptM cases are measured NOT, so always-NOT wins on it while predicting nothing. Balanced accuracy
(the mean of recall on measured DOWN, NOT and UP; always-NOT = 0.333) and macro-F1 over the three classes, over all
scored cases. These measures were chosen after the results; specs/114 registered plain accuracy for korkut.

| case set | model | reachability | always-NOT | macro-F1 (model / reach / always-NOT) | model change calls right |
|---|---|---|---|---|---|
| korkut, 92 pathways | 0.413 | 0.394 | 0.333 | 0.337 / 0.295 / 0.317 | 11.5% of 3,922 |
| korkut, top-level | 0.396 | 0.387 | 0.333 | 0.360 / 0.309 / 0.318 | 10.0% of 4,418 |
| decryptM kinobeads, 92 pathways | 0.394 | 0.387 | 0.333 | 0.374 / 0.355 / 0.292 | 28.6% of 811 |
| decryptM kinobeads, top-level | 0.365 | 0.382 | 0.333 | 0.363 / 0.330 / 0.308 | 14.1% of 1,323 |

On both measures the model is above always-NOT and about level with reachability. When it calls a change on korkut,
it is right about one time in ten.

**Amendment 3 (2026-10-28): P100 re-scored with the corrected reach (specs/113 amendment 1).** Score A on the
top-level build is 0.503 [0.496, 0.509] (was 0.498 [0.493, 0.503]); B 0.511, B − A +0.008 [0.003, 0.012]. The interval
still includes 0.5, so H4 still fails. Output `results/833abd6/p100/` under `20261027-1346_8b6527e_top114`.
