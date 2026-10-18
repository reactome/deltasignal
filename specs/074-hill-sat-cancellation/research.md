# specs/074 research: results of arms hsc and hscd

Solver `501e480`, run 2026-10-18 with `scripts/run_arm.sh`, no overrides.

| arm | build | against | results |
|---|---|---|---|
| hsc | `20261017-1246_5e663d9` (canonical) | `results/6fb4329` | `results/501e480/hsc` |
| hscd | `20261018-0815_ddbd828_dsub` | `results/40adace/dsub` | `results/501e480/hscd` |

## Classes: nothing moves

`holdout_report.py`, with and without `--set-aside RAF_MAP_kinase_cascade`:
held-out +0, tuning +0, experimental +0 on canonical. Per case
(scratch `cmp074.py`), no class changes on either build or axis.

**The registered gate passes** (held-out net 0, experimental 0, no pathway
below 0). It passes because no class changes, which is what a numerical fix
in a band where every value already reads UP should do.

## Values and convergence

| | valid | `pred_ui` moved > 1e-3 | converged, before → after |
|---|---|---|---|
| canonical curator | 23,511 | 141 (TP53 54, Mitotic G1 21, RUNX2 19, PIP3 13) | 22,930 → 23,117 |
| canonical experimental | 845 | 24 (PIP3 13, TP53 7, Mitotic G1 4) | 804 → 824 |
| faithful curator | 23,511 | 617 (TP53 540) | 22,436 → 22,616 |
| faithful experimental | 845 | 68 (TP53 51) | 753 → 767 |

- Canonical curator: 209 cases become converged and 22 stop converging. The
  22 are one perturbation, CBFB overexpression in RUNX2: its residual moves
  from 6.9e-7 to 1.1e-6, across the 1e-6 tolerance, with the same classes.
- specs/066 fallback, canonical curator: 64 → 58 component-solves above
  tolerance, rescued 33 → 47. Experimental: 17 → 12, all 12 rescued.

## Predictions

1. At most 235 class changes: **holds** (0).
2. Converged solves do not decrease; non-converged valid cases do not
   increase: **holds** (581 → 394 curator, 41 → 21 experimental).
3. Class changes concentrated in non-converged solves: **vacuous**, there are
   none.
4. Readouts in the band keep their class: **holds**.

## Reading

The fix moves no class and makes 187 more curator cases and 20 more
experimental cases converge on canonical. It is what specs/073 needs: the
traced TP53 component converges under `DS_NONCONV_POOLS=joint` only with it
(specs/073 spec, probe).
