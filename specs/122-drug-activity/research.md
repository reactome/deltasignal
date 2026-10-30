# specs/122 research: a drug blocks its target's activity

Arms on `20261027-1346_8b6527e_top114` (HPN-DREAM on `20261026-1643_8b6527e_hpn112`), every other setting the default.
- `knockout` (control): korkut and decryptM `s114_top_kb` (code `a47898d`); CPPA `cppa_frozen` (code `8276c28`);
  HPN-DREAM `hpn112_amend4` (code `cb71cef`).
- `activity`: korkut `597463d/a122_activity`; decryptM `f491284/a122_activity_dm` (the first decryptM run read every
  case `no_pin`, a bench bug fixed in f491284: under activity nothing is pinned, so validity now comes from the target's
  roots); CPPA `ff81b6a/cppa_activity` shards 0-2 and `597463d/cppa_activity` shard 3; HPN-DREAM `ff81b6a/hpn112_activity`.
  Code ff81b6a POSTs the redirected network per solve; 597463d onward sends the cut edges by id (`/api/solve`
  `edge_redirects`), which gives identical predictions (all 143 PD0325901 readouts in Developmental Biology, and
  test_api_errors.jl) and is much faster on the largest networks.

## Result

| check | knockout | activity | activity − knockout (95%) | verdict |
|---|---|---|---|---|
| A1: CPPA training per-drug correlation (fold) | – | – | interval above 0 | **pass** (numbers private, CPPA terms) |
| A2: korkut balanced accuracy, each arm's scored cases | 0.396 | 0.379 | [−0.035, −0.004] | **fails** |
| korkut, the 17,836 cases scored in both | 0.396 | 0.391 | [−0.013, +0.003] | no harm |
| A2: decryptM balanced accuracy | 0.365 | 0.362 | [−0.010, +0.008] | no harm |
| decryptM, the 5,096 cases scored in both | 0.365 | 0.364 | [−0.008, +0.009] | no harm |
| A2: HPN-DREAM mean AUROC (30 official contexts) | 0.6930 | 0.6925 | [−0.0010, −0.00002] | **fails** (13 contexts worse, 4 better) |

**A1 passes; A2 fails on korkut and, by a hair, on HPN-DREAM. Not adoptable as registered.**

## Trace (korkut, public)

- **The cases activity adds are the drugs' own targets' readouts**, which a knockout pin excluded: 360 cases, of which
  175 are measured DOWN and activity predicts every one NOT. Several korkut drugs lower their target's own
  phosphorylation: the allosteric AKT inhibitor AKTi 1/2 stops AKT reaching the membrane, so p-AKT falls. These cases
  are most of the loss (each arm's own cases −0.018, common cases −0.005).
- **On common cases** activity removes false alarms (463 DOWN → NOT and 115 UP → NOT where nothing changed) and loses
  some right DOWN calls (36), which weigh more under balanced accuracy since few readouts change.
- **The drug pairs that lose most** contain AKTi 1/2 (AK), the CDK4 inhibitor (RY) or the MEK inhibitor (901).

## Reading

Neither rule fits every drug. An inhibitor of a kinase's catalysis (the MEK inhibitors; most ATP-competitive drugs)
leaves its target's phosphorylation, and can raise it when the target's feedback onto its activators is lost: activity
is right there, and on CPPA's 76 training drugs it helps clearly. An inhibitor that stops its target being activated
(allosteric AKT inhibitors) lowers the target's own phosphorylation as a knockout does. Which one a drug is depends on
its mechanism, not on the network. That is a per-drug property, a natural parameter for specs/121's phase 3 (learned
from the training drugs, with the drug's class as a feature) rather than a fixed rule. The default stays `knockout`;
adoption is Adam's decision.
