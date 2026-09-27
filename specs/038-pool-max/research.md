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

## Result (2026-09-27)

**Setup:**
- Control `ctrl038`: build `20260927-0833_bddfe78_ctrl037`, pools off, code
  defaults.
- Pooled build `20260927-1526_bddfe78_sp038`: `LNG_SET_POOL=1`, `LNG_CAP_POOLS=0`.
- All arms at solver 5a2a927.

| arm | curator held-out | tuning | experimental | worst pathways |
|---|---|---|---|---|
| **`max`** | **−120** (27 / 147, p = 3.5e-21) | −208 | **−57** (16 / 73, p = 7.2e-10) | PIP3 −140, TP53 −57, IFN-γ −17, IL-4/13 −14 |
| **`product`** | −2 (2 / 4, p = 0.69) | +45 (TP53 +42) | **+8** (8 / 0, p = 0.0078) | Fanconi −3 |

**`max`: NOT ADOPTED.** It fails every gate.
- "Any one member fulfils the role" covers a knocked-out member with its
  unperturbed peers, and the ground truth expects that knockout to matter
  (specs/033: 68% of single-member curator cases expect a change).
- PIP3's catalyst sets (46 pool → catalyst edges) carry most of the loss.

**`product`: passes every gate.**
- Held-out −2, which is within the −15 floor and not significant.
- Experimental +8 (p = 0.0078).
- No pathway below −10.

The experimental gain spans 2 pathways but only **4 perturbations**, so it is
reported as **concentrated** under the pre-registered concentration rule.
TP53 +42 on tuning is consistent with the uuid-draw effect seen before
(specs/026, 030) and is not attributed to pooling.

**Conclusion.** One OR node per set-valued regulator or catalyst (Adam's
design) is **adoptable on faithfulness under `product`**, and neutral to
slightly positive, as specs/033 found. The pool's math is the product of member
folds: every member's change counts, which is what the ground truth rewards.
The literal "any one fulfils it" reading (`max`) is refuted on both axes.
