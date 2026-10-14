# specs/055 research: result

Arm `dinert` = build `20261012-1401_01206ca`, commit 4683388 (src/test
identical to a807851), `DS_DRUG_MODE=inert`; results
`results/4683388/dinert`. Control `cj0` (`results/a807851/cj0`). Both arms
converged on the same solves (curator 1,652 of 1,725, experimental 215 of
244). Drugs held: 15,086 curator and 6,648 experimental node-solves (D-P4).

## R1. Numbers (`holdout_report.py --compare ... --set-aside RAF_MAP_kinase_cascade`)

| split | cases | net | fixed / broke | p | pathways moved | perturbations |
|---|---|---|---|---|---|---|
| curator held-out | 19,000 | **+8** | 8 / 0 | 0.0078 | 1 / 71 | 1 (KIT up) |
| curator tuning, all | 5,100 | +13 | 23 / 10 | 0.035 | 2 / 11 | 5 |
| curator tuning, RAF set aside | 5,016 | +3 | 3 / 0 | 0.25 | 1 / 10 | 1 (CDK4 KO) |
| experimental, all | 849 | +2 | 12 / 10 | 0.83 | 2 / 10 | 5 |
| **experimental, RAF set aside** | 800 | **+2** | 2 / 0 | 0.50 | 1 / 9 | 1 (CDK4 KO) |

Outside RAF **nothing broke** on either axis. Every change is:
- **KIT overexpression, Signaling by SCF-KIT** (8 held-out readouts, DOWN →
  UP, right). Reactome curates KIT binding TKIs: "KIT:type II TKIs"
  (lenvatinib, masitinib, sunitinib), "KIT:ripretinib", and "KIT:sSCF
  dimer:KIT:type I KIT binding TKIs". With drugs propagating, raising KIT
  raises those inhibitor-bound complexes, which pull the readouts down.
- **CDK4 knockout, Mitotic G1** (3 curator, 2 experimental, UP → DOWN,
  right). "CDK4 inhibitors bind CDK4" (R-HSA-9754130) makes
  p-T172-CDK4:CCND:CDK4 inhibitors. With drugs propagating, the knockout
  lowers that inhibitor-bound complex and de-represses the readouts.

Both are the specs/032 mechanism: a drug-bound complex follows its protein,
in a cell that has no drug. Neither case was traced by solve; the class is
confirmed in Neo4j.

**RAF (set aside, reported):** curator +10 (20 fixed / 10 broke),
experimental 0 (10 / 10). The switch falls to its low state: ARAF, BRAF, HRAS
and NRAS knockouts read DOWN (right), and HRAS/NRAS overexpressions read
DOWN at about 1e-6 (wrong). This is the collapse of specs/054 T3, not signal (D-P3).

specs/032 measured +30 held-out over 4 pathways (MET, ROCKs, KIT, VEGF).
Only KIT remains on the canonical build; why the others no longer move is
not traced.

## R2. Gate (registered, RAF set aside)

| criterion | result |
|---|---|
| curator held-out net > 0 at p < 0.01 | +8, p 0.0078: **pass** |
| experimental loss not significant (RAF aside) | +2, no loss: **pass** |
| no pathway losing > 10 (RAF aside) | none loses: **pass** |
| every experimental loss outside RAF traced | there are none: **pass** |

Predictions: D-P1 (+5 to +40) held at +8; D-P2 (|net| ≤ 3) held at +2;
D-P3 and D-P4 held.

**The gate passes as registered. The evidence behind it is thin.** The
held-out gain is one perturbation (KIT up) read at 8 places, so its McNemar
p overstates independence; the tool flags it as such. Outside RAF the change
touches 2 perturbations and breaks nothing. Per the rule, the
recommendation is to adopt `DS_DRUG_MODE=inert` as the default, on
faithfulness (the benchmark cells have no drug) more than on numbers. This
is Adam's call. If adopted, RAF cases move from the switch's high state to
its low state: RAF experimental stays 0 net, curator +10.
