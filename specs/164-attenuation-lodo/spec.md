# Feature Specification: Attenuation, leave-one-dataset-out

**Feature Branch**: `specs/164-attenuation-lodo`

**Created**: 2026-11-16

**Status**: Draft, registered BEFORE any arm runs

**Input**: Adam, 2026-11-16 ("It doesn't hurt"), after specs/163. Every magnitude lever has failed at current data, but each attenuation fit so far was judged within one or two datasets (specs/128–131). This is the decisive version: learn how strongly signals carry along edges on two drug-perturbation datasets, and test on the third, which it never saw.

## The lever (no new code)

The specs/128 per-edge-class response curves (`DS_EDGE_CURVES`): activating slots read their input's fold through g(f) = f^h (1 + K^h) / (f^h + K^h). The identity (no table) is today's model. The grid is specs/129's, fixed there: K ∈ {∞, 30, 9, 3} × h ∈ {0.5, 0.75, 1}, minus the identity point, so 11 tables (`bench/analysis/129/cf_K*_h*.tsv`) plus the identity control `cf_off`. Nothing new is fitted beyond choosing one grid point.

## The datasets (all FIT; each scored as its own bench writes it)

- **korkut** (`korkut_bench.py`, scored cases).
- **decryptM** (`decryptm_bench.py`, scored cases; the kinobeads case table `analysis/114/cases_frozen92_kinobeads.tsv` with the Release-98 phospho export `phospho_ewas_r98.csv`, as specs/153 ran its subset on 98).
- **CPPA** (`cppa_bench.py`, single-compound samples × mapped antibodies; the prediction for a (drug, readout) is the mean clipped log2 fold over pathways, called at log2 0.85 / 1.15; measured call at ±log2 1.5). Numbers stay outside the repository (CPPA terms).

**Exclusions (the scorecard rule)**: a drug's own direct markers on korkut (specs/138) and CPPA (specs/139); on decryptM, cases whose phosphosite gene is one of the drug's targets.

## Runs

The served Release 98 build `20261113-1311_baf608b_r98fix` (the gene and site lookups use the Release 98 Neo4j on 7687, so the build must be 98; amended before any run), at one commit, with one replay server per grid point carrying exactly `DS_EDGE_CURVES=<table>` (or none for `cf_off`). Each server's curve rule is verified POSITIVELY from the solve response (`edge_curve_rule`) before any bench runs. All three benches run against each server.

## Leave-one-dataset-out (`bench/analysis/164/score_164.py`, fixture-checked before use)

For each held-out dataset D:
1. **Selection**: on the OTHER two datasets, choose the grid point with the highest mean three-class balanced accuracy. Ties go to the identity, then the larger K, then the larger h.
2. **Test on D**, comparing the chosen point with the identity, under the 2026-11-16 false-alarm gate rule: among the units whose call moves, fixed > broken by a one-sided sign test at p < 0.05; AND reachable sensitivity (reachable meaning the identity's fold is not 1) drops by at most 0.05. Balanced accuracy, specificity, and precision are reported beside it.

**Verdict**: attenuation **generalizes** only if the gate passes on ALL THREE held-out datasets; a pass on one or two is recorded as partial. If the selection picks the identity for a held-out dataset, that fold is "no attenuation chosen" (no test).

**Not an adoption test**: specs/129 already recorded that every contracting grid point loses curator held-out balanced accuracy. Even a full pass here would need a curator answer before any default changes (Adam's call). This spec answers whether attenuation learned on drug data transfers across drug datasets.

## Non-goals

- No new curve family, gradient fitting, or per-edge weights (specs/131 found per-edge fitting data-starved).
- No change to any default or served behaviour.
