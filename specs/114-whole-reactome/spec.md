# specs/114: the frozen model on all of Reactome (decryptM, korkut2015, LINCS P100)

**Status:** part 1 (build and coverage) registered 2026-10-27, before anything was built. Part 2 (the join and the
scoring rules) is added before any solve on the joined network.

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

## Part 2

The join, the solve protocol, the scores and the predictions are registered here before any solve.
