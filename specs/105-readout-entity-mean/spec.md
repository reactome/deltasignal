# specs/105: a readout weighs each curated entity once, not each variant copy (pre-registration)

Adam, 2026-10-25: "I want you to figure out the last two". Found tracing EGFR knockout in Signaling by
ERBB2 (4 experimental cases where curators and experiments agree, and MP-BioPath is right).

## Why (traced, canonical `20261024-2211_7f4b28a`)

The readout p-6Y-ERBB2 heterodimers (R-HSA-1963585) is a DefinedSet of three curated members: the
EGFR, ERBB3 and ERBB4 heterodimers. It has no node of its own, so the benchmark reads its producing
reaction, "Trans-autophosphorylation of ERBB2 heterodimers" (R-HSA-1963582), through the generator's
proxy table: 225 reaction copies, one per variant combination. By the entity each copy produces:

| output entity | copies |
|---|---|
| ERBB4 heterodimers (R-HSA-1250316) | 126 |
| ERBB4 heterodimers (R-HSA-1250329) | 63 |
| ERBB3 heterodimer (R-HSA-1963583) | 27 |
| EGFR heterodimer (R-HSA-1963593) | 9 |

`DS_KO_AGG=mean` averages all 225, so EGFR weighs 9/225 and its knockout reads 0.96 (NORMAL) at
1963585, 1963588, 109783 and 179838. Its overexpression reads 5–100×. The weight is how many ligand
and isoform combinations the generator enumerated, not biology. `bench/benchmark_mpbiopath_cases.py`
already states the principle for set readouts ("a member split into three positional occurrences
outweighs one split into two purely by how finely decomposition happened to cut it, which is an
artifact and not biology"), and specs/091 applied it to pool resting copies. The benchmark of record
(`benchmark_vs_mpbiopath.py`) does not apply it.

## The rule (benchmark side)

New value `DS_KO_AGG=mean_entity` (default stays `mean`). The readout's uuids are grouped by the
curated entity each one carries; the prediction is the mean over groups of each group's mean.
- Readout resolved to its own nodes (the direct route): one group. Identical to `mean`.
- Resolved to producing- or consuming-reaction proxies: grouped by each reaction copy's non-cofactor
  output entities (`diagram_entity_id` in `nodes.csv`).
- A set resolved to its members' nodes: grouped by member.

Also: an unknown `DS_KO_AGG` value becomes an error. It used to fall through to `max` silently.

This changes what is measured, not the model. If adopted, headlines are reported under both readings,
as specs/102 is.

## Gate (as specs/077)

Curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10; read with RAF set aside too.

## Predictions

1. EGFR knockout in ERBB2: the 4 agreed experimental errors read ≤ 0.80 (DOWN) and are fixed; EGFR
   overexpression stays UP.
2. Only readouts resolved through proxies or set members move; direct-route readouts are bit-identical.
3. Experimental net ≥ +2 (643/849).
