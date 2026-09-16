# Solver defaults, re-measured on the wide curator set

**Created**: 2026-09-16
**Status**: Landed
**Supersedes**: the default choices in `specs/002-upregulation-propagation` R5

## Why this exists

The code defaults and `docker-compose.dev.yml` disagreed on four variables.
Production (`docker-compose.prod.yml`) sets no solver variables, so it ran the
code defaults; every benchmark ran through the dev container. `specs/002`
(2026-09-10, `17bc6be`) wrote its winning values into the code and never touched
the compose file, which still carried values from 2026-07-16 (`c9805b2`).

The tests pinned the code defaults, so `test_config_validation.jl` **failed
inside the project's own dev container** — and because a failing top-level
`@testset` aborts the file, 141 of its 151 assertions silently never ran.

## Method

Wide curator set, Release97, one shared catalog build (`cat_fresh`), paired and
conditioned. Each arm changes ONE variable from the config that won the
head-to-head, and is compared against it. Every arm verified before use:
container env dumped, catalog mount asserted, and raw `pred_ui` checked to
differ from the baseline — a config change that produces bit-identical output
has not applied.

Two earlier attempts at this comparison were discarded, for the record:

1. A broken health-check (`read -t 2 < /dev/zero` returns instantly) let the
   benchmark start against a booting server.
2. `--force-recreate` drops `PATHWAY_CATALOG`, silently swapping the catalog the
   server serves. Two builds of one pathway share NO uuids, so every observation
   matched nothing and every readout reported baseline — a complete 23,908-case
   run with ONE distinct predicted value, scored and reported as meaningful.

Both now have guards: PR #27 (benchmark aborts on catalog mismatch) and PR #29
(the API reports or rejects observations naming absent nodes).

## Head-to-head

| config | scored | correct | accuracy | macro-F1 |
|---|---|---|---|---|
| code defaults (`hill_sat`, assembly off, eps 1e-12) | 22,902 | 18,931 | 0.8266 | 0.7939 |
| dev compose (`hill_log`, assembly on, eps 1e-3) | 22,902 | 19,179 | 0.8374 | 0.8031 |

768 predictions changed, net **+248** for the dev-compose config. Raw values
differed on 57.3% of cases, so the arms genuinely differ.

## Attribution — one variable at a time

From the dev-compose config, moving each variable toward the code defaults:

| variable moved | changed | net |
|---|---|---|
| `DS_ASSEMBLY_LIMITING` 1 → 0 | 613 | **−196** |
| `DS_AND_MODE` `hill_log` → `hill_sat` | 171 | **−90** |
| `DS_INHIBITOR_EPS` 1e-3 → 1e-12 | 269 | **+11** |

Sum −275 against −248 measured jointly; they do not compose additively.

**Assembly-limiting is the driver.** This restores the result already on file
from 2026-07-16, where assembly-limiting was the first real curator win.

## The epsilon sweep

All at the winning config, paired against eps 1e-3:

| `DS_INHIBITOR_EPS` | correct | accuracy | macro-F1 | net |
|---|---|---|---|---|
| 1e-3 | 19,278 | 0.8374 | 0.8035 | — |
| 1e-4 | 19,292 | 0.8380 | 0.8045 | +14 |
| 1e-6 | 19,289 | 0.8379 | 0.8041 | +11 |
| 1e-12 | 19,289 | 0.8379 | 0.8041 | +11 |

**1e-6 and 1e-12 are bit-identical** — same 269 changes, same +109/−98. At or
below 1e-6 the epsilon has saturated into a pure divide-by-zero guard, which is
exactly what `specs/006` FR-002 required and nobody had confirmed.

1e-4 scores 3 better in 23,022. It is **not** taken: at 1e-4 the value still
changes results, so it is still acting as a model parameter, and choosing it
optimises a guard constant against the evaluation set.

## The experimental ground truth agrees

`specs/002` chose on the experimental set, so overturning it on the curator set
alone would be picking the benchmark that agrees. Measured directly, same
catalog, epsilon held at 1e-12 in both arms:

| defaults | scored | correct | accuracy | macro-F1 |
|---|---|---|---|---|
| old (`hill_sat`, assembly off) | 846 | 597 | 0.7057 | 0.6116 |
| new (`hill_log`, assembly on) | 846 | **616** | **0.7281** | **0.6521** |

Net **+19** (+26/−7); per pathway −4 TP53, −1 WNT, +1 ERBB2. `propagator_missed`
falls from 89 to 84.

So the change is not a ground-truth artifact: it wins on both. That does not
retroactively make `specs/002` wrong on its own evidence — it ran a different
comparison on 564 cases against a different baseline and a catalog that no
longer exists (its own research.md records the crash). It does mean the choice
does not survive re-measurement on either set today.

## Decision

    DS_AND_MODE          = hill_log   (was hill_sat)
    DS_ASSEMBLY_LIMITING = true       (was false)
    DS_INHIBITOR_EPS     = 1e-12      (unchanged; compose corrected from 1e-3)
    DS_HILL_SAT_EPS      = 1e-5       (unchanged; inert under hill_log)

Code and compose now agree, so the test suite passes in the dev container and
production computes what the benchmark measures.

## The unresolved tension

`specs/002` R5 was not wrong about the curve. `hill_sat` really does implement
the stated AND intent — 10x10 reads 99.98 against `hill_log`'s 74.06, and a lone
node at UI 50 reads 50.00 against 41.43. `hill_log` compresses throughout the
operating range, which by that argument is a defect.

It also scores 90 cases better on 23,022 curator cases.

So the compression is doing useful work that "correct" multiplication does not,
and **nobody has explained why**. Candidate: compression damps the multi-hop
amplification that drives false-change calls, which is the dominant error mode
(see `project_over_coupling_is_the_target`, and the reach/false-change
table in `specs/008-cycle-handling`). That is a hypothesis, not a finding.

Until it is explained this is an **empirical default, not a principled one**.
`test/test_and_curves.jl` names the testset accordingly rather than claiming the
defaults implement the design intent.

## What this says about the 564-case set

`specs/002` chose `hill_sat` and assembly-off on +31 over 564 experimental
cases. On 23,022 curator cases those choices cost 286 between them. Both
measurements are real; they answer different questions against different ground
truths. FR-008 exists because the small set has now reversed at scale
repeatedly, and this is the clearest instance yet — it silently set the
repository's defaults for six days.
