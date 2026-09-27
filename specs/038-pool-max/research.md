# specs/038 — Set-valued regulators and catalysts as one OR node, read as `max`

**Status:** pre-registered 2026-09-27, before the build or any arm.

## Why (Adam, 2026-09-27)

"A reaction should be made unique based on its inputs and outputs, not the
regulator. The regulator should be merged into an OR node, and then have that
regulate the reaction, because any of them fulfil the regulation of the
reaction. We just have to figure out the math."

The OR node is the set pool of specs/033 (`LNG_SET_POOL`, draft LNG #99). Its
math was measured under `product`, `extreme`, `geomean` and `mean`, but not under
the literal OR: **`max`**, the strongest member fulfils the role.

- A knocked-out member is covered by its unperturbed peers.
- A raised member lifts the pool.

The ground truth leans the other way on knockouts (specs/033: a single set
member is expected to change the readout in 68% of curator cases). So this is a
real question, not a formality.

**Scope:** pools for catalyst and regulator sets only (`LNG_SET_POOL=1`,
`LNG_CAP_POOLS=0`). Inputs keep their virtual-reaction expansion, as Adam
specifies.

## Arms

All at one solver commit, through `scripts/run_arm.sh`:

| arm | build | settings |
|---|---|---|
| `ctrl037` | `20260927-0833_bddfe78_ctrl037` (the fresh control, pools off) | code defaults; already run at 6f261e2 and re-run at this commit |
| `pmax` | a new build at LNG `bddfe78`, `LNG_SET_POOL=1`, `LNG_CAP_POOLS=0` | `DS_SET_POOL_MODE=max` |
| `pprod` | the same build | `DS_SET_POOL_MODE=product` (the specs/033 reference) |

## Adopt `max` only if all hold (against `ctrl037`)

- curator held-out net ≥ −15 and not significantly negative;
- experimental net ≥ 0;
- no pathway, tuning or held-out, loses more than 10 cases;
- any gain spans at least 2 pathways and 5 perturbations (or is reported as
  concentrated).

**Predictions:**
- `max` loses single-member knockout cases (those were 68% "expected change");
- `max` gains on up-regulation;
- the net depends on which dominates;
- `pprod` reproduces specs/033's +6 / +10 / +3 within the regeneration noise
  floor.

**Also reported:**
- the set-member-case table of specs/033 for each arm (KO and UP separately);
- per-pathway net;
- convergence counts.
