# specs/123 research: a drug acts by its mechanism

## Preview by composition (2026-10-30), before the confirming runs finish

CPPA, decryptM and HPN-DREAM solve each drug on its own, so a `mechanism` prediction is exactly the specs/122
`knockout` or `activity` prediction for that drug; the composition below uses those runs (CPPA `8276c28/cppa_frozen`,
`ff81b6a`+`597463d/cppa_activity`; decryptM `a47898d/s114_top_kb`, `f491284/a122_activity_dm`; korkut
`a47898d/s114_top_kb`, `597463d/a122_activity`). The HPN-DREAM `mechanism` run (`ac2db81/hpn112_mechanism`) gives
predictions identical to `ff81b6a/hpn112_activity` file for file, as composition says. The CPPA `mechanism` run at `c00308e`
completed shard 2 of 4 before a memory reap stopped every run (2026-10-30); its 10,735 cases are identical to the
composed predictions for the same pathways (status, value and call), which checks the composition. The korkut run,
needed for the drug pairs that mix the two actions, was stopped too; it was re-run alone and is reported below.

| check | knockout | mechanism | mechanism − knockout (95%) | verdict |
|---|---|---|---|---|
| M1: CPPA training per-drug correlation | – | – | interval includes 0 | **fails** (numbers private) |
| M1, `DS_DRUG_MECHANISM_UNKNOWN=activity` (reported) | – | – | interval above 0 | would pass |
| M2: decryptM balanced accuracy (18 drugs) | 0.365 | 0.359 | [−0.010, +0.002] | no harm |
| M2: HPN-DREAM mean AUROC | 0.6930 | 0.6925 | [−0.0010, −0.00002] | **fails** (AZD8055 is labelled activity, so this is specs/122's result) |
| M2: korkut balanced accuracy, each arm's scored cases | 0.396 | 0.402 | [+0.002, +0.010] | no harm |

On the 46 CPPA training drugs with a known label, `mechanism` and `activity` score almost the same; the gain `activity`
has over `knockout` on all 76 drugs is mostly in anonymised compounds the labels leave `unknown`.

## korkut, the real run (2026-10-30)

`dab1132/a123_mechanism_korkut` on `20261027-1346_8b6527e_top114`, one server, every other setting the default. Scored
as specs/122 (balanced accuracy over each arm's scored cases, paired bootstrap over drugs, seed 123, 2,000 resamples);
the same script reproduces specs/122's activity figures exactly (0.396 → 0.379, [−0.035, −0.004]).

| korkut subset | cases (knockout / mechanism) | knockout | mechanism | mechanism − knockout (95%) |
|---|---|---|---|---|
| all sets, each arm's scored cases (**M2**) | 17,836 / 17,860 | 0.396 | 0.402 | +0.006 [+0.002, +0.010] |
| all sets, cases scored in all three arms | 17,836 | 0.396 | 0.402 | +0.006 [+0.002, +0.010] |
| sets of one action (both drugs knockout, or both activity) | 10,021 / 10,037 | 0.393 | 0.400 | +0.007 [+0.002, +0.013] |
| mixed sets (one drug of each; 32 pairs) | 7,815 / 7,823 | 0.402 | 0.407 | +0.004 [+0.001, +0.007] |

Against `activity` (each arm's cases): 0.379 → 0.402, +0.023 [+0.012, +0.041]. The knockout-labelled AKTi 1/2 (AK)
and SB216763 (ST) keep their pins, so their own targets' readouts stay unscored as under knockout (24 more cases than
knockout, 336 fewer than activity). Converged: knockout 16,639 of 17,836, mechanism 16,663 of 17,860.

**Verdict on M2: korkut no harm (positive), decryptM no harm, HPN-DREAM fails by a hair** (its one drug, AZD8055, is
labelled activity, so its result is specs/122's). With M1 failing as registered, `mechanism` is **not adoptable as
registered**; the decision is Adam's.

## Traced: drugs whose `activity` loss the mechanism reading did not explain (CPPA training)

- **Fulvestrant** (a degrader, labelled knockout): ER-alpha falls; the label fixes it.
- **Topotecan**: TOP1 is its most-lowered readout (camptothecins trigger TOP1 degradation); labelled activity.
- **Dasatinib**: its own readouts include Src's autophosphorylation and its CSK site, which fall; the knockout's
  exaggerated fall of the target's own forms carried the per-sample correlation.
- **Palbociclib**: a defect in the activity rule. In Cell Cycle (R-HSA-1640170) cutting CDK4/6's catalyst and
  regulator edges predicts p-Rb S807/S811 **33× up**; knockout predicts it down, as measured. Traced in Mitotic G1
  (frozen build, same rails: hyperphosphorylated RB1 55×, p-S807/p-S811-RB1 82×). Reactome curates the cyclin D:CDK4/6
  set (R-HSA-8949896) both as the catalyst of RB1 phosphorylation and as a **negative regulator** of "RB1 binds and
  inhibits E2F1/2/3:DP1/2" (the regulation is the consequence of that phosphorylation). The activity rule cuts both: 25
  negative-regulator edges are redirected to blocks held at 0, which under `divide` de-represses their reactions to the
  10× ceiling, and the RB1 / E2F / cyclin E:CDK2 loop rails up (cyclin E:CDK2 makes hyperphosphorylated RB1 without
  CDK4/6). Cutting a target's *negative* regulations treats "the drug removes the target's inhibition of X" as a
  full release; a knockout pin does the same, but it also removes the target's other forms, which kept the loop low.
  The exact chain around the rail is not yet followed node by node.
- **AMG900**: four samples; knockout's upward calls on mitotic markers (FoxM1, cyclin B1) happen to match.
