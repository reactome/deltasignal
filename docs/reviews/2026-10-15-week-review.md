# Week review, 2026-10-15

An adversarial review of everything behind the canonical numbers changed
since 2026-10-08: deltasignal `5979e48..adopt-061-defaults` (94 files) and the
generator `6990015..01206ca` (LNG #109–#111, which built the canonical
catalog). Seven reviewers, one per dimension; every finding was checked by
two independent skeptics told to refute it; a completeness critic looked for
what nobody examined. 54 agents. 20 findings upheld by both skeptics, 2 by
one, 1 refuted by both.

## Findings and what was done

| # | severity | finding | disposition |
|---|---|---|---|
| 1 | high | The specs/061 source gate (adopted that day) is a switch read on the live iterate; inside a cyclic component it toggles between sweeps and the solve never converges. SMAD7 KO in TGF-beta swings 0.3×–38× within ±5 iterations of the cap. | Checked every discordant case of the adoption at caps 500–4000 (specs/061 R5): the held-out +32 is real (33/33 fixes converged and stable); the tuning +112 and experimental +7 are mostly iteration-cap artifacts (111/113 and 6/7 fixes are unconverged TP53 solves). Amendment 3 registered and implemented: the gate is decided once from a first ungated solve. Arm `rtwin58g4`. |
| 2 | high | Generator: an output complex with one bound and one shared slot enumerates the bound slot, then overwrites it: 2,905 duplicate reaction nodes in 36 pathways; duplicates weigh 2–6× in OR-means and pushed 4 reactions over the variant cap. | LNG PR #114 (dedupe choices before the cap; test fails on the old code). Catalog variant build `_dedup` and re-measurement. |
| 3 | medium | A carrier shared by several pools (the 12 RAS GAPs in RAF) is written by each pool from the other pools' release flux; RAF does not converge. | specs/049 amendment 7 (exclude every pool's release steps); test fails on the old code. Arm `carrierfix`. |
| 4 | medium | Non-converged solves are scored like converged ones; only a total is printed. | Per-case `converged` and `residual` columns in the case dumps; API returns `residual`. PR #123. |
| 5 | medium | Uploaded bundles cannot carry their tables, so an upload solves without the drug, pathogen, carrier, pool, containment and cofactor rules. | Upload accepts the tables as optional parts; documented in API.md; test. |
| 6 | medium | Neither API nor bench records the own-product or twin rule, now defaults. | API reports `depletion_own_product` and `own_source_depletions`; bench prints a `Depletion:` line (in ARM.json) and refuses a server that does not report them. PR #123. |
| 7 | medium | `holdout_report.py` leaves cases that lose coverage out of the net figure without saying so. | It names them, with how many the baseline had right, and breaks down every unpaired case; test. PR #123. |
| 8 | low | Under the new pin default, a depletion-fed dead end could be pinned as a root (0 on the canonical build). | Requires an out-edge; test. PR #123. |
| 9 | low | Bench POSTed payload omitted `pools` and `node_leaves`. | Added. |
| 10 | low | decryptM protocol line omits the pin setting; korkut predict and score axes overwrite one file; `cycle_carriers_held` counts a shared carrier per pool. | Fixed. PR #123. |
| 11 | low | Pool modifiers read activators at baseline, which would disable the live source gate there (latent). | Moot under amendment 3: a gated pair is gated for the whole solve. |
| 12 | low | Generator `_sig_cache` persists across pathways and can read prefetch-limited caches. | Open: clearing it per pathway may change the catalog, so it needs its own build and measurement. |
| 13 | low | Generator docstring of `_leaf_multiset` contradicts the code. | Fixed (LNG #114). |
| 14 | low | Docs: API.md drug default; spec status lines for 058–061; `LNG_CAP_POOLS` default; asserting-file count. | Fixed. |
| critic | low | The determinism suite's relabel test drops every bundle table. | Relabel test under the twin and source rules with containment kept. Pools/carriers already had relabel tests in `test_cycles.jl`. |

## Gaps the critic named and not yet examined

- `LNG_BIND_STOICH` (default on): the claim that it moves no prediction was not re-checked by regeneration.
- LNG `pathway_reaction_ids` now drawn from the full connection table.
- `benchmark_vs_mpbiopath.py` defaulted to `--max-edges 20000` while catalog.sh and run_arm.sh pass 100000, so a direct run skipped large pathways. Fixed: the default is 100000.
- The curator pathway `_NEW_ROCK_signaling_regulates_MRLC_phosphorylation` is not in `catalog_pathways.tsv` (skipped as "no network").
- CLI parse → solve with the new JSON fields is not executed by any test.
- LNG `neo4j_connector.py` was reviewed only through its callers.
- A catalog-wide count of scored cases from non-converged solves (possible now that the case dumps carry convergence).
