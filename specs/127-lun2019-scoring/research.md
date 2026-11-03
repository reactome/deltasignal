# specs/127 research: values-blind predictions, committed before the seal is broken

Nothing here has read the Lun 2019 data. The archive stays sealed (mode 000); `bench/lun_bench.py` never opens it.
These are step 1 of the spec's pre-scoring list for the frozen model; the learned model's predictions are added when
specs/125 reaches stage 5.

## Frozen model predictions (2026-11-03)

Solver commit `0d16eb9` (`run_arm.sh --axis lun_predict`, protocol of record: `DS_PIN_SCOPE=root_cycle`,
`DS_PIN_ROOT_IGNORE=depletion`, over-expression pin 80×). Its `src/` differs from tag `freeze-2026-10-26` only by
specs/122's `edge_redirects` (an API request field, unused here) and specs/124's `DS_BASELINE` at its default 0.01, whose
control arm reproduced the frozen run case for case.

| build | file | sha256 |
|---|---|---|
| `20261027-1346_8b6527e_top114` (primary, defines P) | `results/0d16eb9/lun127_top/lun_predict_cases.tsv` | `0c0ec0339c5e25eea347d1bc8346cf5913e94fab48d1e5a9cb59bc1c48515d59` |
| `20261026-1349_8b6527e` (frozen 92 pathways, also reported) | `results/0d16eb9/lun127_frozen/lun_predict_cases.tsv` | `3b0ccfda7d7030657caa1c8bda66cf818a9145ad36e060939c29980b250298bc` |

Combined over networks with `lun_score.combine_networks` (non-self pairs with status `scored`):

| | top114 (population P) | frozen 92 |
|---|---|---|
| (protein, readout) pairs | **6,953** | 1,659 |
| proteins | 428 | 241 |
| frozen model: NOT / UP / DOWN / conflict | 5,792 / 1,027 / 123 / 11 | 1,137 / 455 / 58 / 9 |
| reachability: NOT / UP | 5,271 / 1,682 | 1,125 / 534 |
| solves converged | 960 of 984 | 773 of 773 |
| case rows on a non-converged solve | 478 of 10,920 | 0 |

Population P is the 6,953 top114 pairs. The 478 rows from the 24 non-converged solves are kept, as every other axis
keeps them; their pairs are listed by the `converged` column, and P is reported with and without them at scoring.

The scoring code (`bench/lun_score.py`, `bench/test_lun_score.py`) is on main at `76d74b7` (PR #209); with this file
merged, steps 1–3 are done for the frozen model. The seal is not broken until stage 5.
