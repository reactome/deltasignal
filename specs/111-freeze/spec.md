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

### Shortlist (from published descriptions only, 2026-10-26; nothing downloaded)

| dataset | perturbations | readouts | access | fit |
|---|---|---|---|---|
| LINCS P100 (Litichevskiy et al., Cell Systems 2018; PXD017458) | 90 drugs × 6 cell lines in the pilot; a later release covers 128 compounds | ~100 phosphopeptides (reduced-representation phosphoproteome) | public, CC BY 4.0, PanoramaPublic / clue.io | largest and openly licensed; drug perturbations need a target map registered before scoring (decryptM, specs/021, lost to reachability under its first map) |
| HPN-DREAM (Hill et al., Nature Methods 2016) | 4 breast cancer lines × 8 ligand stimuli, a few kinase inhibitors (mTOR, AKT, MEK, ...) | ~45 phosphoproteins by RPPA, time courses | challenge data via the DREAM platform (an account and data terms are needed) | well-curated readouts; few perturbations; stimuli and time courses do not map onto our steady state directly |
| Kinase knockouts in HCT-116 (ERK2, PLK1, PIK3CA, PAK4; 2026, TMT proteome and phosphoproteome) | 4 genetic knockouts | proteome and phosphoproteome | public repository not confirmed | genetic knockouts are the closest match to our protocol, but only 4 genes |
| HMMR knockout in KELLY neuroblastoma cells (Mendeley Data 10.17632/6wr2tj8wr9) | 1 genetic knockout | phosphopeptides | public | one gene: too small to test anything |

Adam chooses. An account, data terms, or any download wait for his approval.

### LINCS P100: approved, then set aside on coverage (2026-10-26, before any value was read)

Adam approved LINCS P100 ("feel free to start with LINCS"). The processed files came from GEO GSE101406
(levels 3 and 4, the analyte, instance, perturbagen and cell tables; SHA-512 checked against the series'
`SHA512SUMS`) into `~/deltasignal-catalogs/external/lincs_p100/raw_20261026`. Only the metadata tables were read:

- 89 compounds × 6 cell lines (A375, A549, MCF7, NPC, PC3, YAPC), 3 h, 0.5–10 µM; about 35 compounds act on
  kinases or other signalling proteins the networks contain.
- Of the 96 phosphopeptide probes, **7** match a phosphorylated EWAS in Reactome 97 (UniProt and residue, from
  the graph: `analysis/phospho_ewas_r97.csv`), and **4** occur in the 92 benchmark networks (CDK1 T161 in 5
  pathways, RPS6KA3 S369 in 2, ABI1 S183 and BRAF S446 in 1 each). P100 was designed to cover many proteins,
  not the canonical signalling sites Reactome curates.

Four readouts cannot test the model, so P100 is not scored. No value file was opened; the dataset is still
unseen and keeps its TEST role in `bench/datasets.tsv`.

## Frozen state (recorded after the tags)

- Generator: logic-network-generator `8b6527e`, tag `freeze-2026-10-26`.
- Solver and bench: deltasignal `05dc867`, tag `freeze-2026-10-26`.
- Catalog: `20261026-1349_8b6527e` (92 of 92 pathways, 209,480 nodes, 492,034 edges; Reactome Release 97),
  served by the API since the tag.
- Benchmark: `builds/20261026-1349_8b6527e/results/05dc867` under the protocol line
  `DS_PIN_SCOPE=root_cycle DS_PIN_ROOT_IGNORE=depletion DS_PERTURB_UI_DOWN=0 DS_PERTURB_UI_UP=80 DS_KO_AGG=mean_entity`.
  Every prediction on both axes is identical to specs/105's `entmean` arm on the previous canonical build.

**Development figures** (development set, never generalisation; specs/096 §6):

| axis | cases | accuracy | macro-F1 |
|---|---|---|---|
| curator, all | 24,100 | 89.22% | 0.8693 |
| curator, held-out (SELECT, specs/051) | 19,000 | 89.77% | 0.8687 |
| curator, tuning | 5,100 | 87.18% | 0.8656 |
| experimental (849 cases) | 849 | 75.74% (643; MP-BioPath 643) | balanced accuracy 0.6857 |

Converged: 1,724 of 1,725 curator solves, 244 of 244 experimental.
