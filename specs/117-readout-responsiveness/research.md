# specs/117 research

Run 2026-10-27: `bench/analysis/train_117.py` (tested in `bench/test_train_117.py`) on the specs/114 top-level case
logs (build `20261027-1346_8b6527e_top114`, code `a47898d`, `s114_top_kb`); output
`~/deltasignal-catalogs/analysis/117/train_117.json`. Leave one drug out; balanced accuracy (bal, registered) and
macro-F1 (mF1).

| arm | korkut bal | korkut mF1 | decryptM bal | decryptM mF1 |
|---|---|---|---|---|
| R2: responsiveness + DeltaSignal's call | 0.676 | 0.440 | **0.573** | **0.456** |
| R1: responsiveness only | 0.672 | 0.433 | 0.467 | 0.398 |
| R3: responsiveness + reachability's call | **0.688** | **0.448** | 0.497 | 0.411 |
| frozen model | 0.396 | 0.360 | 0.365 | 0.363 |
| reachability | 0.387 | 0.309 | 0.382 | 0.330 |
| always-NOT | 0.333 | 0.318 | 0.333 | 0.308 |
| specs/116 trained layer (fold, hops) | 0.389 | | 0.431 | |

korkut: 33,488 test predictions (a combination once per drug it contains), 12 drugs; decryptM: 5,096, 18 drugs.

95% bootstrap intervals over held-out drugs:

| difference | korkut | decryptM |
|---|---|---|
| R2 − R1 | [+0.001, +0.006] | [+0.019, +0.166] |
| R2 − R3 | [−0.033, +0.002] | [+0.018, +0.164] |
| R1 − frozen | [+0.234, +0.329] | [−0.008, +0.228] |
| R2 − frozen | [+0.238, +0.331] | [+0.102, +0.333] |

**Hypotheses.**
- **H1 passes.** R2 beats R1 on both datasets. On korkut the gain is small (+0.004); on decryptM it is +0.106.
- **H2 fails.** On korkut, reachability's call adds at least as much as DeltaSignal's (R3 0.688 vs R2 0.676). On decryptM,
  R2 beats R3.
- **H3 fails.** Responsiveness alone beats the frozen model on korkut by a wide margin. On decryptM its interval includes 0.

**Where R2 differs from R1 on decryptM (traced).** 370 cases fixed and 143 broken:
- **Most of the gain is in one place.** Sites with a moderate changed share in other drugs (0.1-0.25), mostly down. R1
  calls them DOWN; R2 calls NOT where the model's call is NOT or UP. This fixes 265 unchanged cases.
- **The UP calls are an artefact of the rule.** In the cell (changed share > 0.5, mostly down, model NOT), R2 calls UP.
  UP is rare (119 of 5,096), so a few UP cases in that cell outweigh the rest under the balanced-accuracy rule. This
  fixes 41 UP cases and breaks 94 DOWN cases.

Recall by class, R1 → R2: unchanged 2,957 → 3,192, down 393 → 344, up 9 → 50. Macro-F1, which this artefact moves
less, also rises (0.398 → 0.456).

**Reading.**
- **Responsiveness is the information training was missing.** Learned only from other drugs, it lifts balanced
  accuracy from about 0.40 to 0.67 on korkut and to 0.47 on decryptM.
- **DeltaSignal's call adds information on decryptM, beyond both responsiveness and reachability.** The model's NOT
  marks which responsive sites a given drug leaves alone.
- **On korkut, DeltaSignal adds almost nothing once responsiveness is known**, and less than reachability does. korkut's
  72 antibodies are mostly central signalling readouts that respond to most drugs. Which drug was given matters less
  there than which readout was measured.

Both datasets are FIT, and the readout tables are fitted on them; these are development results. A readout never seen
in training gets no responsiveness, so the layer only helps on panels that overlap readouts measured before.
