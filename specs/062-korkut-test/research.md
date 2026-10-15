# specs/062 research

## R1. Step 1, values-blind (no korkut2015 value read): informative by the registered rule

`scripts/run_arm.sh kp_ctrl --axis korkut_predict` and `scripts/run_arm.sh
kp_rtwin58g2 --axis korkut_predict` (the specs/061 overrides) at `93f78f7` on
the canonical build; results `results/93f78f7/kp_{ctrl,rtwin58g2}`. The bench
ran with `--predict-only` (protocol line: "PREDICT ONLY (values unread)").

- 11,626 cases in each arm (identical sets), 11,062 scored in both; 1,449 of
  1,493 solves converged in each.
- **178 cases change predicted class (≥ 10): step 1 passes**, so step 2 runs
  as registered.
- Caveat, recorded before the values are read: the 178 are **16 distinct
  (pathway, readout, change) units from 3 drugs in 2 pathways**, repeated
  across the combinations that contain the drug:
  - AKTi (AK) in Transcriptional Regulation by TP53 (R-HSA-3700989): p21 and
    ATR UP → DOWN, AMPK pT172 NOT → UP;
  - temsirolimus (Tm) in TP53: p53, p21, TSC2, TSC2 pT1462, ATR, AMPK pT172
    NOT → UP;
  - Nutlin (NT) in R-HSA-8878159 (Transcriptional regulation by RUNX3): p53,
    p21, BIM, cyclin D1, SMAD3, SMAD3 pS423, β-catenin UP → NOT, c-Myc
    DOWN → NOT.
  The case-level McNemar treats them as independent; they are not. The
  registered clustered-by-drug-set report is the honest reading, and with
  three drugs it can hardly reach significance.

## R2. Step 2, scored once: **inconclusive**

`scripts/run_arm.sh ks_ctrl --axis korkut` and `ks_rtwin58g2` (the specs/061
overrides), both at `8de3c7f`, results `results/8de3c7f/ks_{ctrl,rtwin58g2}`.
(A first `ks_rtwin58g2` start refused: `verdict.py` was uncommitted; it was
moved aside so both arms ran the same commit.) `expr.csv` checked before
scoring: 89 × 99, all numeric, −5.08 to +1.90, about half negative,
consistent with log2 ratios, so no amendment. Verdict by
`bench/analysis/062/verdict.py`. 11,062 cases scored in both arms (no
conflicts at the primary threshold); 1,449 of 1,493 solves converged in each.

| | primary (|log2| ≥ log2 1.5) | secondary (log2 1.15), not the verdict |
|---|---|---|
| measured down / not / up | 820 / 10,049 / 193 | 2,752 / 6,443 / 1,866 |
| control accuracy, macro-F1 | 0.638, 0.333 | 0.499, 0.372 |
| arm accuracy, macro-F1 | 0.640, 0.332 | 0.501, 0.374 |
| always-DOWN | 0.074, 0.046 | 0.249, 0.133 |
| discordant: arm right / control right / both wrong | 86 / 68 / 24 | 80 / 51 / 47 |
| net, McNemar p | +18, p 0.17 | +29, p 0.014 |
| by drug set: positive / negative / zero | 11 / 18 / 0 (p 0.26) | 14 / 8 / 6 (p 0.29) |

**Registered verdict: inconclusive** (net > 0 but p ≥ 0.05). The secondary
threshold's "confirms" is reported, not adopted as the verdict.

What the cases say (16 units, step 1):
- At the primary threshold 91% of measured calls are NOT, so the arm gains
  mostly where it stops predicting change (Nutlin in RUNX3: the control's
  UP/DOWN calls on p53, p21, BIM, cyclin D1, SMAD3, β-catenin, c-Myc become
  NOT, and most are measured NOT) and loses where it starts predicting change
  (temsirolimus and AKTi in TP53: AMPK, ATR, TSC2, p21 become UP or DOWN,
  measured NOT).
- At the secondary threshold one unit agrees with the arm's new call:
  temsirolimus → p53 UP is measured UP in 9 of 10 drug sets.
- Not part of either verdict, and recorded because it frames everything
  else: at the primary threshold **both arms are far below always-NOT**
  (0.908 accuracy; the model predicts change in most cases, the data shows
  little). The registered baselines were always-DOWN and reachability; the
  reachability baseline was **not computed** (the bench records pin counts,
  not pin sets), an omission against the registration.

### What this means

- The independent check neither confirms nor contradicts specs/061's
  candidate. It touches 3 drugs in 2 pathways, none of them the perturbations
  the candidate was built for (SMAD2/3, MDM2 KO): its evidence on the
  candidate's own mechanism is close to nil, as the registration predicted.
- korkut2015 is spent as a TEST set and is re-registered FIT
  (`bench/datasets.tsv`), as specs/051 requires.
- The adoption decision rests on the SELECT and experimental evidence
  (specs/061 R3) and on Adam's judgement of the post-hoc risk; the next
  Reactome release remains the clean replication when it exists.
