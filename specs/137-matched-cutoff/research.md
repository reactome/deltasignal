# specs/137 research

## R1. korkut, matched cutoff (2026-11-05)

`bench/analysis/score_137.py` on `builds/20261026-1349_8b6527e/results/e4809b8/k134_abund_a2/korkut_cases.tsv` (code
`e4809b8`, the default since specs/134 amendment 2), t = 1.5. Output `~/deltasignal-catalogs/analysis/137/korkut_M.json`.
Same 10,469 cases in both arms (9,905 with RAF set aside).

| | balanced acc. | accuracy | macro-F1 | change calls (wrong) |
|---|---|---|---|---|
| C (0.85 / 1.15) | 0.4469 | 0.8117 | 0.4136 | 1,788 (1,408) |
| M (1.5) | 0.4486 | 0.8251 | 0.4201 | 1,638 (1,263) |
| M − C | +0.0017 [−0.0022, +0.0046] | +0.013 | +0.007 | 150 moved: 145 fixed / 5 broken |
| M − C, RAF set aside | +0.0019 [−0.0022, +0.0048] | +0.014 | +0.007 | the same 150 |

**The gate fails** (interval includes 0, both figures). Per the registration, the CPPA training check runs only after
a pass, so it was not run.

Every move is a change call turning NOT: DOWN → NOT 115 (111 measured NOT), UP → NOT 35 (34 measured NOT). Half are
total readouts, half phospho, spread over many pathways (Signaling by EGFR 44 the most).

**Why so little moves: the change calls are not near the cutoff.** Change calls by predicted fold:

| predicted fold | change calls | wrong |
|---|---|---|
| ≤ 0.01 (read through to zero) | 853 | 607 |
| 0.01 – 0.1 | 61 | 38 |
| 0.1 – 0.667 | 431 | 328 |
| 0.667 – 0.85 (M moves these) | 115 | 111 |
| 1.15 – 1.5 (M moves these) | 35 | 34 |
| 1.5 – 10 | 145 | 144 |
| ≥ 10 | 148 | 146 |

Almost every call between the two cutoffs is wrong (145 of 150), so matching does what it should there, but only 8%
of the change calls sit in that band. Half of all change calls (853) are predicted at 1% of baseline or less, and 71%
of those are wrong: a knockout read straight through to the readout, which no cutoff reaches.

**Threshold-free.** AUROC of the predicted log fold: DOWN against the rest 0.679 (RAF set aside 0.666); UP against the
rest **0.427** (0.422), below chance. Of 193 measured UP cases the model calls 4 UP, 55 DOWN, 134 NOT; of its 293 UP
calls, 290 are wrong. Upward changes in korkut (drug-induced rises, mostly feedback relief across pathways) are not in
the model's single-pathway networks.

**S (truth at 1.15 instead).** Balanced 0.4031, accuracy 0.586 on 10,468 cases: a 15% truth turns many measured noise
changes into calls the model does not make. Not a gain either way.

**CV diagnostic.** Leave-one-drug-set-out picks 1.5 for every drug set; CV gain over 1.15 +0.0017. In-sample balanced
accuracy by t: 1.15 0.4469, 1.25 0.4485, 1.5 0.4486, 2 0.4436, 3 0.4361. The data would choose the registered threshold
too, and it gains nothing worth having.

**Predictions.**
- Q1 (accuracy above 0.85) — wrong: 0.825.
- Q2 (balanced about +0.01, interval above 0) — wrong: +0.0017, interval includes 0.
- Q3 (fewer than 50 of the 409 single-copy gene-expression total calls move) — not checked on that exact set; 75 of all
  527 total change calls move, and 209 of the 527 are predicted at ≤ 0.01, which no cutoff reaches.
- Q4 (PIP3 switch not fixed) — holds: none of the 153 PIP3 DOWN calls measured NOT moves.
- Q5 (CV picks ≥ 1.25) — holds: 1.5.
- Q6 (DOWN AUROC above 0.5, well below 1) — holds, 0.679. Not predicted: UP is below chance.
- Q7 (CPPA) — not run (gate failed).

## What this says

The cutoff mismatch is real but small. The calls it touches are nearly all wrong, but there are 150 of them out of
1,788. korkut's errors are dominated by readouts the model drives to zero (853 calls) or to a rail (148 at ≥ 10×), not
by modest predicted changes, and by measured rises the model has no route for. specs/135's and specs/136's cases
mostly sit at those extremes (SRC knockout's AKT pS473 at 0.21 is below 0.667 too). A matched cutoff is coherent and
costs nothing, so it could be adopted as protocol on faithfulness, but it does not move the gate.
