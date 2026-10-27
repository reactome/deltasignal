# specs/114: the frozen model on all of Reactome (decryptM, korkut2015, LINCS P100)

**Status:** part 1 (build and coverage) registered 2026-10-27, before anything was built; result in research.md.
Part 2 (solve protocol and scores) registered 2026-10-27, before any case table was built or any solve run.

Adam, 2026-10-27: the external tests used 92 pathways (13 for HPN-DREAM), so coverage was bound to be low; "this
needs to be done again with all the reactome pathways". Not HPN-DREAM: its 13 pathways were already every pathway
with a reaction on MTOR. This is for the other three datasets.

**What this is and is not.** decryptM, korkut2015 and LINCS P100 are FIT (each scored once). A score here is a
development result about the network, not a validation of the model, and it is reported next to the first score of
each dataset, never in its place. The model (generator and solver at `freeze-2026-10-26`) is not changed.

## What was used before

| dataset | pathways | how |
|---|---|---|
| decryptM (specs/021) | the 92 benchmark pathways | separately; a site scored only in a pathway holding a target |
| korkut2015 (specs/062) | the 92 benchmark pathways | separately, per drug set and pathway |
| LINCS P100 (specs/113) | the 92 benchmark pathways | separately (reachability only) |

The one earlier look across pathways (specs/021, 2026-10-13, graph only, 92 pathways joined by node identity):
148 of decryptM's 154 unreachable regulated sites became reachable, but so did 96% of unchanged sites. No solve was
run.

## Part 1: build and coverage

**Pathways** (`top_level_pathways.tsv`, 28): every human `TopLevelPathway` with `isInDisease` false, Reactome 97
(Adam, 2026-10-27: "only generating top level pathways"). The generator walks `hasEvent*`, so each network holds its
whole hierarchy: every reaction is in one network (a reaction listed under two top-level pathways is in both), and
sub-pathways inside one top-level pathway are already joined. Only links between the 28 networks need a join (part
2). The largest are Signal Transduction (2,589 reactions), Metabolism (2,380) and Immune System (1,900), several
times any network built so far. The Disease top-level pathway is left out (as specs/112): mutant and pathogen states
the cell lines in these datasets do not have.

*Superseded first list:* every non-disease pathway holding a reaction directly (1,728, commit 75c7a30). It
duplicated each reaction across the hierarchy; its build (`20261027-1337_8b6527e_all114`) was stopped after 45
pathways, is incomplete, and is not used.

**Build:** `scripts/catalog.sh build --variant top114 --list specs/114-whole-reactome/top_level_pathways.tsv` from
the generator at `8b6527e` (`freeze-2026-10-26`). A variant build never moves `current`. A pathway that fails to
build (or exceeds memory or time) is listed and left out, not retried with other settings; such a failure is
itself a result about the frozen generator, and Adam decides what follows.

**Coverage (graph only, no measured value read),** for each dataset, three counts of its readouts (decryptM: scored
sites; korkut: the 73 used antibodies; P100: the 96 probes):
1. in a network of the 92-pathway frozen catalog (as scored);
2. in a network of the all-pathway build;
3. reachable from a node holding a perturbed target in the joined network (part 2's join).

## Part 2: solve protocol and scores (registered 2026-10-27, before any case table is built or any solve run)

**Protocol (Adam, 2026-10-27: option A).** Each of the 28 top-level networks is solved on its own, exactly as the
92 benchmark pathways are: no join, no new rule. Sub-pathways are already joined inside each top-level network. Option
B (all 28 merged into one network by node identity) is not run; earlier joins lost accuracy (specs/056, the July 2026
stitch). The 28 networks are never joined, so coverage count 3 is reported per network: a readout reachable from a
pinned target in the same top-level network.

**Model.** Solver and bench at this branch's HEAD when the arms start: `src/` unchanged since `freeze-2026-10-26`,
code defaults, served by `scripts/run_arm.sh --catalog`. No setting is changed in response to any result.

**Two builds per dataset.** Each dataset is scored on the top-level build `20261027-1346_8b6527e_top114` **and** on the
frozen 92-pathway build `20261026-1349_8b6527e`, with the same code. decryptM and korkut were first scored with
earlier models, so the frozen model has no score on the 92 pathways yet. Without that control a difference would mix
the model change with the network change. `--max-edges` is raised above the largest top-level network (147,774 edges),
so that no network is skipped.

**decryptM** (`bench/decryptm_bench.py`, the specs/021 protocol unchanged):
- Case tables from `bench/analysis/decryptm_coverage.py` on each build, with `all_sites.tsv` and `shp099_sites.tsv`,
  the Release 97 phospho export (`phospho_ewas_r97.csv`) and `target_ewas.csv` re-exported from Release 97 by
  specs/051's query. The primary target map is the bench's default, the kinobeads map (specs/051 R4); the hand map
  (`DRUG_TARGETS=hand`, specs/021's) is reported as well.
- A case is (drug, network, site). The top-level pathways are not in `TUNING_PATHWAYS`, so every case is in the
  `held-out` split; on the 92-pathway build the held-out split is reported as before, and all cases are reported too.
- Reported: regulated cases right (model, reachability, always DOWN), unchanged cases right, macro-F1, McNemar by
  drug and the by-drug 95% interval of model minus reachability, as `summarize` computes them.

**korkut2015** (`bench/korkut_bench.py`, the specs/062 protocol unchanged): a case is (drug set, network, readout).
Reported per build: scored cases; accuracy of the model, of reachability (the bench's `reach_predicted`) and of
always-NOT; regulated cases right; McNemar of model against reachability over scored cases.

**LINCS P100** (`bench/p100_bench.py`): the model is not scored, since only 7 of the 96 sites are in any network.
Score A (curated reach, protein level, specs/113) is computed on the top-level build, with the registered AUROC over
contexts, and set beside the 92-pathway value (0.499). The script also computes B and C; they are not part of this
test and are not reported as results.

**Cases in more than one network.** A reaction under two top-level pathways is in both networks, so one (drug, site)
can be scored more than once. The case unit above is kept for comparability with the original benches. Also reported:
the number of distinct (drug, site) units, and how many of them get different model calls from different networks.

**What would count.**
- H1: the top-level build scores more decryptM and korkut cases than the 92-pathway build (expected from part 1).
- H2: on the top-level build, the frozen model gets more regulated decryptM cases right than reachability does
  (kinobeads map). On the 92 pathways an earlier model lost (86 vs 115 of 278).
- H3: on korkut, the frozen model's accuracy on the top-level build is above both reachability's and always-NOT's.
- H4: P100 score A on the top-level build is above 0.5, with a bootstrap interval over compounds that excludes 0.5.
Each is reported whichever way it falls. These are development results on FIT data (spec preamble), not a validation.
