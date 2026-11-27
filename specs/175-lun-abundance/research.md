# specs/175 research (2026-11-26)

Arm `lun175_abund`: code `356cb2b2`, build `20261027-1346_8b6527e_top114`, `run_arm.sh --axis lun_predict`; protocol line
identical to specs/127 apart from `DS_TOTAL_READOUT=abundance_a3`. Gene lookup on Neo4j Release 97 (checksum
501466562143, the build's), reached by container address because the usual container's network would not attach.
Output `results/356cb2b2/lun175_abund/lun_predict_cases.tsv` (sha256 `21e287c3…`); 960 of 984 solves converged, as in
specs/127. Scoring: `bench/analysis/175/score_175.py` at `1aa46476`; output outside the repo
(`~/deltasignal-catalogs/analysis/175/score_a.json`). Half B was not read.

## Control

**0 differences** on phospho-readout rows: the same rows, the same `predicted` and `reach_predicted`, and `pred_ui`
within 1e-6. HEAD's solver reproduces `0d16eb9` on these networks. Both files have the same 11,028 (pathway,
protein, readout) rows.

## Gate (registered, half A)

| | |
|---|---|
| pairs scored | 3,127 (same set as specs/173) |
| moved | 42, all total β-catenin |
| fixed / broken | **35 / 0**, sign test p = 2.9e-11 |
| recall of measured changes | 0.194 → 0.194 |
| **gate** | **passes** (paired gate and sensitivity guard) |

The 42 moves:
- 35 UP → NOT where nothing moved (fixed);
- 4 UP → NOT where β-catenin fell (wrong before and after);
- 3 UP → DOWN where nothing moved (wrong before and after).

β-catenin calls: before, 44 UP and 98 NOT; after, 2 UP, 3 DOWN and 137 NOT. Measured: 1 UP, 9 DOWN, 132 NOT.

**Reported (three-class balanced accuracy):**

| | old rule | abundance |
|---|---|---|
| model, all pairs | 0.379 | **0.384** |
| reachability, all pairs | 0.372 | 0.380 |
| model, β-catenin | 0.232 | 0.321 |
| reachability, β-catenin | 0.182 | 0.321 |

## Reading

The readout defect specs/174 traced is fixed: the abundance rule removes 35 false β-catenin calls and breaks none. It
does not make the model better than reachability: 0.384 against 0.380. Reachability gains as much as the model,
because most β-catenin pairs lose their route to an abundance copy. Per the registration, `abundance` stays the Lun
bench default, as on korkut and CPPA. Half B stays sealed for the full candidate protocol.
