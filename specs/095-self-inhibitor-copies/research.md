# specs/095 research: arm sic

Canonical `20261022-2011_233665c`, solver `dfb7160` with `--server DS_SELF_INHIBITOR_COPIES=on`,
against `results/8162f5f/dflt091` (main's defaults). Results `results/dfb7160/sic`. Read with
`holdout_report.py --set-aside RAF_MAP_kinase_cascade --pair-on-valid`; cases with
`specs/089-ub-ligase-unit/diff_cases.py`. The solves added 2,822 (curator) and 1,262
(experimental) copy pair-solves.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | 0 (0 / 0) | 1 |
| curator tuning | **+14 (14 / 0)** | 0.0001 |
| experimental (RAF set aside) | **+9 (9 / 0)** | 0.0039 |

Every changed case is PTEN overexpression in PIP3 activates AKT signaling: 83.5× → 0.17× (one
readout 100× → 0×), DOWN as the curators and the experiments say. The other pathways in the
census (R-HSA-446728, TP53 transcription, R-HSA-2559583, RAF) move no class. Converged solves
unchanged (1,724 / 1,725, 244 / 244).

**The gate passes** (held-out 0, experimental +9, no pathway below 0).

## Predictions

1. At least 6 of the 9 PTEN overexpression experimental cases fixed: **holds, 9 of 9**.
2. TP53 in PIP3: no prediction; its 4 experimental cases do not move.
3. No prediction for the nets: tuning +14, experimental +9, held-out 0.

## Reading

Recommended for adoption. The rule is the specs/022 reading (an inhibitor built from its own
reaction's input counts that input at weight w) extended to another node copy of the same
entity, which the generator emits when a species is both a root and a product (PTEN mRNA). It
fixes the largest class of experimental errors on which every reference agrees, and breaks
nothing on either axis. The 115 slots it adds elsewhere change no class.

## Decision

Adam, 2026-10-23: "yes", adopt. `DS_SELF_INHIBITOR_COPIES=on` is the default; `off` restores
the earlier pair test.
