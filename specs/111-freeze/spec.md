# specs/111: the freeze

**Status:** in progress, 2026-10-26.

Adam, 2026-10-26: "go ahead with your plan" (freeze under specs/096 §3, then one external dataset under §4).

## Why now: the stop rule (specs/096 §1)

Condition 2, "four consecutive pre-registered arms fail their gate", was met in specs/107 and not flagged at
the time:

| arm | spec | gate |
|---|---|---|
| `m107` (metabolite pools) | 107 rule 1 | fails (B4 −2 pairs) |
| `ms107` (strict multi-step carriers) | 107 rule 2 | fails (B4 −2 pairs) |
| `mw107` (wide multi-step carriers) | 107 rule 3 | fails (held-out −50) |
| `mc107` (compartment copies) | 107 amendment 1 | fails (held-out −10) |
| `cat109`, `cat109b` (catalyst loops) | 109 | passes; no class change on any axis |
| `dedup110` (self-modification counted once) | 110 | fails (held-out −60) |

Every loop rule since specs/106 improved the behaviour checks and none improved accuracy; specs/110 showed
that current correct calls lean on a structural double count. More rules now mostly risk overfitting.

## Decisions at the freeze (Adam, 2026-10-26)

- specs/105 `DS_KO_AGG=mean_entity`: **adopted** as the bench default (gate passed: experimental +2, tuning +6,
  held-out 0). Headlines are reported with it and with `mean`.
- specs/104 `DS_POOL_DRIVE_TWINS=hold`, specs/109 `DS_SELF_FED_MODE=catalyst_loop`: **not adopted** (no class
  change on any axis).
- Every other measured-but-not-adopted mode stays at its default.

## What is frozen

- **Generator:** logic-network-generator `8b6527e`, tag `freeze-2026-10-26`, every `LNG_*` at its code default
  plus the build environment recorded in the catalog's `BUILD.json`.
- **Solver and bench:** deltasignal at the commit tagged `freeze-2026-10-26` (the merge of this spec), every
  `DS_*` at its code default as listed in CLAUDE.md "Solver Configuration".
- **Catalog:** one canonical build from the tagged generator against Reactome Release 97 (recorded below once
  built and verified).
- **Benchmark protocol:** `DS_PIN_SCOPE=root_cycle`, `DS_PIN_ROOT_IGNORE=depletion`, knockout pin 0
  (`DS_PERTURB_UI_DOWN=0`), overexpression pin 80× (`DS_PERTURB_UI_UP=80`), `DS_KO_AGG=mean_entity`; calls
  DOWN < 0.85 ≤ NORMAL < 1.15 ≤ UP.

No rule, default or threshold changes after the freeze for the external test.

## Backlog, for after the external test (recorded, not run)

- `LNG_POOL_REGEN=protein` (exists since specs/048 lever B2, never measured on its own): a pool path counts its
  enzyme regenerated when the same proteins come out in another complex. Would make the TGF-beta SMAD2/3,
  insulin receptor and macroautophagy cycles pools (specs/110 traces).
- A flow-direction rule for large catalyst-only loops (EGFR, WNT, HOX activation, DSB/DNA repair; specs/109).
- Ku-type carrier cycles (NHEJ, Fanconi) under the multi-step carrier rules (specs/107), with the specs/109 gate.

## External dataset (specs/096 §4)

A shortlist is drawn from descriptions only and given to Adam. Nothing is downloaded, and no account or terms
are accepted, before he approves one.
