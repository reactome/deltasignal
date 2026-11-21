# specs/169 research: the measurement (2026-11-21)

Scripts at `06b71066`. The registration was `55f97aa7`, plus two harness fixes before any truth was read: a pathway
whose only target is the perturbed gene itself sent an empty readouts query, and the gene lookup was made one pass
(equal to the bench's lookup on a check set). A dedicated server served `20261113-1311_baf608b_r98fix`. Outputs
(outside the repo): `~/deltasignal-catalogs/analysis/169/` (predictions, truth, score_169.json).

**Scope:** 1,240 split-A perturbed genes pinnable somewhere; 3,383 solves; 30,031 (pathway, gene, target)
predictions. 1,028 genes were engaged (knockdown ≤ 0.5, ≥ 25 cells), giving 24,054 scored cases. Split B was not read.

**Truth check:** for each engaged gene, its own transcript's fold computed by `truth_169.py` matches the dataset's
`fold_expr`: median 0.166 vs 0.165, correlation 0.989, n = 864.

## Result

| truth label | reached / cases | changed, of reached | changed, of unreached | call rate | precision | lift [95%] | declined that changed |
|---|---|---|---|---|---|---|---|
| **1.5× (primary)** | 5,814 / 24,054 | **0.057** | 0.057 | 0.77 | 0.060 | +0.003 [−0.001, +0.008] | 0.045 |
| 15% (secondary) | 5,814 / 24,054 | 0.334 | 0.334 | 0.77 | 0.348 | +0.014 [+0.006, +0.023] | 0.286 |
| 1.5×, model called at 1.5× (descriptive) | | 0.057 | | 0.65 | 0.063 | +0.006 [+0.001, +0.012] | 0.045 |

Among calls on readouts that really changed, the direction is right 60% (1.5×) and 47% (15%) of the time.

## The reading (registered): base rate among reached ≤ 0.25, comparable to the drug panels

On an unselected genetic panel, a reached target changed 5.7% of the time, against korkut's 14–23%. The gene axes'
88–94% (specs/167) comes from how those records were made, not from genetic perturbations carrying through curated
routes. specs/165's "accurate for gene perturbations" does not hold as a property of the model.

The panel is weaker than korkut in one further way. Reach is not enriched for change at all here: 0.057 reached
against 0.057 unreached (korkut: 0.142 against 0.049).

## Post hoc (described, not registered): do the curated regulations themselves hold in K562?

Pairs where the knocked-down gene is a component of a curated regulator of the target's expression event (Reactome
98, one hop): 104 pairs, 5 changed ≥ 1.5× (0.048), against 0.060 for all other pairs. Of the 5, 4 moved in the curated
direction. Most curated transcription regulations do not act in unstimulated K562 within the screen's window.
Possible reasons: cell-type or stimulus dependence, redundancy, small effects. That is a context problem that no
propagation rule can fix, and it agrees with the failed cell-context attempts (specs/121, 130, 139).

## Caveats

- **One cell line.** K562 is a leukaemia line under no stimulus, so many curated programmes (TP53 in a TP53-null
  line, immune and developmental transcription) are expected to be off.
- **Readout type.** The readouts are transcripts; the model reads the curated product node.
- **Truth noise.** The 15% label is close to pseudo-bulk sampling noise (33% of unreached readouts "change").
- **Pins.** Knockdown was pinned to 0, though the median measured knockdown is 0.17.

Split B (`replogle2022_b`) stays sealed for any model or rule this result leads to.
