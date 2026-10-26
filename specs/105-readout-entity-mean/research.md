# specs/105 research

## Arm `entmean` (canonical `20261024-2211_7f4b28a`, results `results/54a1761/entmean`; control `results/23a526c/dflt102`)

`scripts/run_arm.sh entmean --bench DS_KO_AGG=mean_entity`. Solver unchanged; the protocol line
reads `DS_KO_AGG=mean_entity`.

| axis | net | fixed / broke | p | RAF set aside |
|---|---|---|---|---|
| curator held-out | 0 | 0 / 0 | 1 | 0 |
| curator tuning | +6 | 6 / 0 | 0.03 | +6 |
| experimental | **+2** | 2 / 0 | 0.5 | +2 |

Experimental 641 → **643/849** (75.74%; MP-BioPath 643). Curator all 89.20% → 89.22% (with the
specs/102 supplement). Every fixed case is in Signaling by ERBB2: EGFR and ERBB3 knockout at the
ERBB2 heterodimer readouts (R-HSA-1963580, 1963585, 1963588), 0.96 / 0.88 → 0.75. Nothing broken.
**Gate passes.**

Values changed in 10 curator readouts (56 cases) across 8 pathways (ERBB2, RAF, IL-2, IL-3/5/GM-CSF,
FGFR1–4), all set-valued readouts read through proxies or members; only ERBB2's change class.

### Predictions

1. **Partly holds.** 2 of the 4 agreed EGFR knockout experimental errors are fixed (1963585, 1963588,
   0.75). The other two (PI(3,4,5)P3 179838, RAS:GTP 109783) resolve to their own nodes, so the rule
   does not touch them. EGFR overexpression stays UP (25.75×).
2. **Holds** by construction: a direct-route readout is one group, identical to `mean`.
3. **Holds.** Experimental +2, 643/849.

### Reading

This changes the measurement, not the model. Each curated member of a set readout (and each entity a
proxy reaction copy produces) weighs once, which is the rule the case-level benchmark already applies
to set readouts and specs/091 applied to pool resting copies. ERBB2's p-heterodimer readouts read 1/4
EGFR, 1/4 ERBB3 and 1/2 ERBB4 (two ERBB4 entities), where the plain mean read 9/225 EGFR.

Circularity: the rule was found tracing EGFR, an agreed error. The 4 ERBB3 knockout curator cases and
the 2 EGFR ones at 1963580 were not used to choose it. If adopted, headlines are reported under both
readings, as specs/102's are.

**Decision (Adam, 2026-10-26, at the freeze, specs/111): adopted as the bench default** (`DS_KO_AGG=mean_entity`). Headlines are reported with it and with `mean`.
