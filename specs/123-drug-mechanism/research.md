# specs/123 research: a drug acts by its mechanism

## Preview by composition (2026-10-30), before the confirming runs finish

CPPA, decryptM and HPN-DREAM solve each drug on its own, so a `mechanism` prediction is exactly the specs/122
`knockout` or `activity` prediction for that drug; the composition below uses those runs (CPPA `8276c28/cppa_frozen`,
`ff81b6a`+`597463d/cppa_activity`; decryptM `a47898d/s114_top_kb`, `f491284/a122_activity_dm`; korkut
`a47898d/s114_top_kb`, `597463d/a122_activity`). The HPN-DREAM `mechanism` run (`ac2db81/hpn112_mechanism`) gives
predictions identical to `ff81b6a/hpn112_activity` file for file, as composition says. The CPPA and korkut `mechanism`
runs at `c00308e` are in progress and will replace this section.

| check | knockout | mechanism | mechanism − knockout (95%) | verdict |
|---|---|---|---|---|
| M1: CPPA training per-drug correlation | – | – | interval includes 0 | **fails** (numbers private) |
| M1, `DS_DRUG_MECHANISM_UNKNOWN=activity` (reported) | – | – | interval above 0 | would pass |
| M2: decryptM balanced accuracy (18 drugs) | 0.365 | 0.359 | [−0.010, +0.002] | no harm |
| M2: HPN-DREAM mean AUROC | 0.6930 | 0.6925 | [−0.0010, −0.00002] | **fails** (AZD8055 is labelled activity, so this is specs/122's result) |
| korkut, drug sets of one mechanism only (preview) | 0.392 | 0.400 | [+0.002, +0.013] | – (mixed sets running) |

On the 46 CPPA training drugs with a known label, `mechanism` and `activity` score almost the same; the gain `activity`
has over `knockout` on all 76 drugs is mostly in anonymised compounds the labels leave `unknown`.

## Traced: drugs whose `activity` loss the mechanism reading did not explain (CPPA training)

- **Fulvestrant** (a degrader, labelled knockout): ER-alpha falls; the label fixes it.
- **Topotecan**: TOP1 is its most-lowered readout (camptothecins trigger TOP1 degradation); labelled activity.
- **Dasatinib**: its own readouts include Src's autophosphorylation and its CSK site, which fall; the knockout's
  exaggerated fall of the target's own forms carried the per-sample correlation.
- **Palbociclib**: a defect in the activity rule. In Cell Cycle (R-HSA-1640170) cutting CDK4/6's catalyst and
  regulator edges predicts p-Rb S807/S811 **33× up**; knockout predicts it down, as measured. Under trace.
- **AMG900**: four samples; knockout's upward calls on mitotic markers (FoxM1, cyclin B1) happen to match.
