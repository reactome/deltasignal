# specs/104 research

## Arm `twins` (canonical `20261024-2211_7f4b28a`, results `results/f122706/twins`; control `results/23a526c/dflt102`)

`scripts/run_arm.sh twins --server DS_POOL_DRIVE_TWINS=hold`.

| axis | net | fixed / broke |
|---|---|---|
| curator held-out | 0 | 0 / 0 |
| curator tuning | 0 | 0 / 0 |
| experimental | 0 | 0 / 0 |

No case changes class on either axis. Values moving more than 5%: 31 curator and 4 experimental,
all in Signaling by WNT (APC 15, AMER1 15, CTNNB1 1). The other eight pathways that have twin copies
(EGFR, IFN, DNA repair, RAF, …) do not move. Converged 1,724 / 1,725, as before. **The gate passes
trivially.**

### Predictions

1. **Holds.** APC and AMER1 overexpression read the readout R-HSA-2130284 at 1.11 (was 0.014),
   R-HSA-195297 likewise. NORMAL, where curators and experiments say UP: still wrong. Free β-catenin
   (R-HSA-448839) falls to 0.011 (was 0.0001).
2. **Holds.** Experimental net 0.
3. **Holds, more narrowly than registered.** Only WNT moves, not all 9 pathways.

### Reading

A faithfulness fix: an 80× artefact on the pool's exit drive is gone. It changes no class, like
specs/083 and specs/091. Adoption on faithfulness is Adam's call. It does not close the experimental
gap (641/849). The remaining distance on APC/AMER1 is the release pool's conservation: the curated
exit to SCF(β-TrCP) is read from an intermediate and is not a pool transition, so capture rate cannot
raise the ubiquitination flux above the release flux.
