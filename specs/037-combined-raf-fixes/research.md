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

## Revision after review (2026-09-27, before any arm on the combined build)

An adversarial review of both integration branches found a blocking generator
bug and weaknesses in this pre-registration. The first combined build was
stopped part-way, and `20260927-0821_d98d2f8_combined037` is incomplete and
unused.

**Generator (LNG bddfe78):**
- Cap pools fired in reactions the cap never touched: 410 pools in PIP3, which
  has 0 cap hits. The miR-93 RISC input's four subunits were averaged, so a
  miR-93 knockout read 0.75x instead of 0.
- Now:
  - pooling happens only in capped reactions;
  - leaves are grouped by alternative, a complex alternative is an AND unit, and
    a pool needs at least 2 alternatives;
  - shared leaves stay direct;
  - the curated stoichiometry is kept;
  - depletion onto a pool is retargeted to its members.
- Verified by regeneration: RAF's binding reaction has its 4 curated inputs,
  and PIP3 has no input pools (its 11 pools are the specs/033 catalyst and
  regulator pools).

**Pre-registration, revised:**
1. **Control is a fresh regeneration**, not the canonical build. It is built
   from the same LNG commit with `LNG_SET_POOL=0` and `LNG_CAP_POOLS=0`. `base`
   and `combined` were otherwise different uuid draws, and the gates sit at the
   regeneration noise floor.
2. **The experimental gate excludes RAF.** `mean` and dedup were chosen from
   the RAF HRAS trace, and RAF is an experimental pathway, so gating on it would
   select on the evaluation set. The gate is experimental net ≥ +15, p < 0.05,
   **excluding RAF/MAP kinase**. RAF is reported separately as the traced case.
3. **Concentration gate:** the curator held-out gain or loss, and the
   experimental gain, must span at least 2 pathways and at least 5 distinct
   perturbations. Otherwise McNemar is not valid and the result is reported as
   such.
4. **The tuning pathways are gated too:** no pathway, held-out or tuning, may
   lose more than 10 cases. PIP3 (specs/035's failure) is covered.
5. Corrections to the table above:
   - the specs/033 +6 / +3 was measured under `product`, not `mean`, so it does
     not carry over;
   - this spec replaces specs/036's plan (its `cap_*` arms and its choice of
     mode on the tuning split) with `mean` fixed in advance.
6. **Known limitations, stated:**
   - cap pools are minted per reaction copy, not shared;
   - pool `member_leaves` in `nodes.csv` lists every leaf of the set;
   - de-duplication does not see through a pool.

## Result — NOT ADOPTED (2026-09-27)

**Setup:**
- Control: `20260927-0833_bddfe78_ctrl037`, a fresh regeneration with pools off.
  Its structure matches the canonical build exactly (72,025 nodes, 225,147
  edges).
- Pooled build: `20260927-0845_bddfe78_comb037` (72,348 nodes, 214,234 edges).
- All arms at solver 6f261e2.

| arm | curator held-out | curator tuning | experimental (all) | experimental excl. RAF | RAF experimental |
|---|---|---|---|---|---|
| control | — | — | — | 555 / 796 | 14 / 49 |
| **combined** (pools + drugs + dedup, `mean`) | **−162** (53/215, p = 2.6e-24) | −106 | −3 | **−12** | 23 / 49 |
| pools only (`mean`) | −35 (p = 6.9e-7) | −29 | +6 | −2 | 22 / 49 |
| pools + drugs | −4 (p = 0.74) | −36 | +7 | −2 | 23 / 49 |
| pools + dedup | −193 (p = 1.7e-40) | −104 | −6 | −12 | 20 / 49 |

**Gates for `combined`:**
- experimental excluding RAF ≥ +15: it is −12, **FAIL**;
- held-out: −162, **FAIL**;
- no pathway below −10: **FAIL** (IFN α/β −112, TP53 −80, IFN-γ −34, PIP3 −28).

**Not adopted.**

**What the attribution says:**
- **De-duplication is the damage.** pools + dedup reproduces the combined
  losses: IFN α/β −112, TP53 −80, IFN-γ −34, PIP3 −28. That confirms and enlarges
  specs/012's negative (−15 there). Squaring catalyst-and-substrate is
  load-bearing on these networks, even though it amplifies RAF's loop.
- **Pools with `mean`** cost held-out −35 (DSB −13, NOTCH1 −7, TP53 −37 tuning).
  They move RAF 14 → 22.
- **Pools + drugs** are neutral on held-out (−4) and experimental (+7; −2
  excluding RAF), with the specs/032 held-out gains (MET +11, ROCKs +8, VEGF +6,
  KIT +6) offset by the pool losses.
- **RAF gains in every pooled arm (14 → 20–23 of 49).** The cap pools partly
  fix the traced case, but nothing generalises beyond RAF on the experimental
  axis.

**Record:**
- The generator's cap pools stay behind `LNG_CAP_POOLS` (default off). They are
  correct, and they fix the structural defect (RAF's binding reaction has its 4
  curated inputs), but they are not adoptable on accuracy with `mean`.
- A `product`-mode cap-pool arm was not run. With cap pools under `product`, the
  pool re-creates an AND over alternatives, which is the defect itself.
