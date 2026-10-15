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
