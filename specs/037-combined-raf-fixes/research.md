# specs/037 — Combined structural fixes (from the RAF trace), measured catalog-wide

**Status:** pre-registered 2026-09-27, before the build or any arm.

## Why

Tracing RAF/MAP kinase iteration by iteration (specs/034 §8, specs/036) found
compounding amplifiers. Each one alone was neutral or harmful catalog-wide, or
did not move RAF:

| # | amplifier | fix | alone |
|---|---|---|---|
| 1 | variant-capped set inputs, with every alternative required at once | `LNG_CAP_POOLS=1` (specs/036) | RAF still collapses |
| 2 | set-valued catalysts/regulators as separate AND terms | `LNG_SET_POOL=1` (specs/033) | held-out +6, experimental +3 |
| 3 | drug-bound complexes inhibiting the kinase steps | `DS_DRUG_MODE=inert` (specs/032) | held-out +30, experimental −5 |
| 4 | catalyst = substrate squaring | `DS_DEDUP_ACTIVATORS=1` (specs/012) | held-out −15 |

On RAF alone, 1 + 3 + 4 with `mean` pools turned HRAS KO from "no signal" into
0.016x (correct). HRAS up read 0.83x (still wrong, but no longer a collapse).

## Arms

All through `scripts/run_arm.sh` at one solver commit:

| arm | build | settings |
|---|---|---|
| `base` | canonical `20260926-1221_590301c` | code defaults |
| **`combined`** | a new build from LNG `feat/combined-037` (`LNG_SET_POOL` default on, `LNG_CAP_POOLS=1`, ships `drugs.csv`) | `DS_SET_POOL_MODE=mean`, `DS_DRUG_MODE=inert`, `DS_DEDUP_ACTIVATORS=1` |

Attribution arms, on the same new build, are **reported, not gated**:

| arm | settings |
|---|---|
| `pools_mean` | `DS_SET_POOL_MODE=mean` |
| `pools_drugs` | + `DS_DRUG_MODE=inert` |
| `pools_dedup` | + `DS_DEDUP_ACTIVATORS=1` |

## Pre-registration (`combined` against `base`)

**Adopt only if all hold:**
- experimental net ≥ +15 with p < 0.05 (the axis this targets: we trail
  MP-BioPath's published result by 8.6pp);
- curator held-out net ≥ −15 and not significantly negative;
- no held-out pathway loses more than 10 cases.

**Predictions:**
- RAF experimental rises (at least +10 of 49);
- PIP3 is not harmed (the specs/035 failure mode);
- the attribution arms show which component carries RAF.

`mean` is fixed in advance from the RAF trace. It is not chosen on the
evaluation set; `extreme` still collapsed there.

**Also reported:**
- tuning curator;
- per-pathway net;
- concentration;
- convergence counts;
- `drugs_held`;
- cyclic node count excluding pools;
- the RAF iteration trace on the new build;
- the paired gap against MP-BioPath's supplementary table.
