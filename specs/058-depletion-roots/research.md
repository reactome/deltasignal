# specs/058 research

## R1. Arm `deproot`: +9 held-out, not significant; the gate fails

`scripts/run_arm.sh deproot --bench DS_PIN_ROOT_IGNORE=depletion` at DS
`4a2261b` (solver identical to `64d23cb`), canonical build, results
`builds/20261012-1401_01206ca/results/4a2261b/deproot`. Protocol line:
`DS_PIN_SCOPE=root_cycle DS_PIN_ROOT_IGNORE=depletion`; 1,940 pins over 864
perturbations, 1,016 of them roots. Scored against `results/64d23cb` with
`--pair-on-valid --set-aside RAF_MAP_kinase_cascade`; no case valid in one
arm only.

| axis | net | fixed / broke | p | pathways |
|---|---|---|---|---|
| curator held-out | **+9** | 33 / 24 | 0.29 | 5 / 70 |
| curator tuning | +2 | 2 / 0 (MDM2 OE, Cell Cycle Checkpoints) | 0.5 | 1 / 11 |
| experimental (= RAF set aside) | +1 | 1 / 0 (MDM2 OE) | 1 | 1 / 10 |

**Gate: fails** (held-out p 0.29, not < 0.01; TGF-beta loses 24 > 10).
Not adopted; the setting stays unset.

| perturbation | net | pins before → after |
|---|---|---|
| MET OE (MET) | **+10** (predicted) | 1 → 2 |
| CBL KO+OE (EGFR) | +11 | 1 → 2 |
| PTCH KO+OE (Hedgehog) | +8 | 1 → 4 |
| BRCA1, BARD1 (DSB) | +4 | 4 → 6, 2 → 3 |
| MDM2 OE (Checkpoints, tuning + exp) | +2 / +1 | 1 → 2 |
| SMAD2, SMAD3 KO+OE (TGF-beta) | **−24** | 4 → 2 |

Every loss is one mechanism, the one specs/057 R1 found. Plain SMAD2
[cytosol] (R-HSA-3006353) is fed only by a depletion (SMAD2:SMURF2 from
"SMURF2 ubiquitinates SMAD2"), so under the rule it is a true root, and the
specs/024 fallback, which had pinned the recycled signal-carrying forms
(including p-S465,S467-SMAD2), no longer applies. In this pathway plain
SMAD2 enters only SMURF2 binding; the signaling route takes SMAD2 inside the
receptor complex variants of R-HSA-171180 (produced by R-HSA-170868), which
the gene's own pins never reach. The same happens for SMAD3.

### What this means

- The rule is faithful to MP-BioPath's protocol and gains in five
  perturbations (+33), as predicted for MET. It loses TGF-beta's SMAD2/3
  cases (−24), which score only because the gene has no true root and the
  fallback pins its produced forms.
- So the SMAD2/3 TGF-beta cases depend on the fallback: any change that
  gives SMAD2 a true root loses them (specs/057 R1 and here). Whether the
  curators' SMAD2 means the cytosolic pool or the receptor-bound signal form
  is a protocol question, not a network one.
- Nothing is adopted.
