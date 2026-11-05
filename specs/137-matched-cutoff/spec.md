# specs/137: call a change at the size the dataset measures one

**Status:** registered 2026-11-05, before any rescoring (Adam: "the miss match should probably have its own spec and be
looked into"). Follows specs/135 and specs/136, which both ended at the mismatch.

## The mismatch

The model calls a readout changed when its predicted fold leaves [0.85, 1.15) (`bm.classify`: UI < 0.85 DOWN, ≥ 1.15
UP). Those cutoffs are MP-BioPath's published convention, inherited for the curator and experimental axes, whose truth
is MP-BioPath's qualitative calls (specs/004 checked they are not tuned: the curve is smooth and 0.85 / 1.15 is not its
peak).

Two measured axes call a change at a different size:

| axis | truth calls a change at | model calls a change at |
|---|---|---|
| korkut2015 (FIT) | \|log2 fold\| ≥ log2 1.5, any dose of the drug set (specs/062) | 0.85 / 1.15 |
| CPPA training, secondary (balanced accuracy) | \|Δ\| ≥ log2 1.5 (specs/121, specs/134 Q6) | 0.85 / 1.15 |

So a predicted 20% fall is scored as a DOWN call against a truth that would only call 50%. specs/135 found 197 total
calls that would sit at 0.848 under a 24 h stock rule, and specs/136 found SRC knockout's AKT pS473 at 0.80 under PL:
both are changes the model thinks are modest, scored as calls a 1.5× truth cannot share.

Not affected: CPPA's primary (a continuous fold, no cutoff), HPN-DREAM (AUROC), decryptM (truth is a curve-fit
regulated call, not a fold threshold), curator and experimental (truth has no size; MP-BioPath's convention stays).

## The rule (no fitting)

**M (matched).** On an axis whose truth calls a change at a stated fold threshold t, the model's call uses the same t,
symmetric in log as the truth's own `call()`: DOWN if log2(pred) ≤ −log2 t, UP if ≥ log2 t, otherwise NOT. For korkut
and CPPA secondary, t = 1.5 (UI ≤ 0.6667 or ≥ 1.5). The threshold is the dataset's registered one, not chosen here.

Rescored offline from the stored `pred_ui` of existing runs; no solver or network change, so no new solves.

## Arms

- **C (control):** korkut `builds/20261026-1349_8b6527e/results/e4809b8/k134_abund_a2` (the default since specs/134
  amendment 2), predicted at 0.85 / 1.15, measured at 1.5.
- **M:** the same cases, predicted at 1.5.
- **Reported only, not gated:**
  - **S** — predicted at 0.85 / 1.15, measured at 1.15 (korkut's existing `measured_secondary`). It aligns the truth to
    the model instead; the truth changes, so it is not compared to C as a gain.
  - **Threshold-free** — AUROC of the predicted log fold for DOWN against the rest and UP against the rest (no cutoff
    moves it); says how much of the error is the cutoff and how much is ranking.
  - **CV diagnostic** — a symmetric cutoff chosen from {1.15, 1.25, 1.5, 2, 3} by leave-drug-set-out CV on korkut
    balanced accuracy, with its CV gain over 1.15. If it lands near 1.5, matching is what the data would choose
    anyway; if not, that is reported, and M is still the rule (it is not replaced by the fitted value).

## Gate

korkut, M minus C on the same cases: three-class balanced accuracy, 95% interval from resampling drug sets
(`score_134.py`'s bootstrap); the interval must be above 0, read with RAF set aside, both figures reported. Raw
accuracy is reported beside it (fewer change calls raises it almost by construction, since always-NOT scores 0.913).
Then CPPA training secondary under M against C: must not fall (interval), numbers private (CPPA terms); CPPA primary
must be identical. Adoption is Adam's; it would make t a per-axis bench setting.

**Honest limit.** The question came from korkut's errors and korkut is FIT, so a korkut pass is development evidence,
even though nothing is fitted. A confirmation needs data none of this has touched (the CPPA test portion or Lun 2019,
both sealed; opening either is Adam's decision).

## Predictions (before any rescoring)

- **Q1.** Accuracy rises (fewer change calls): korkut 0.812 → above 0.85.
- **Q2.** Balanced accuracy (C: 0.4469 on 10,469 cases) rises, by about +0.01; the interval is above 0 with RAF kept. Less sure with RAF set aside.
- **Q3.** M does not fix specs/135's total calls on gene-expression copies: a knocked-out required regulator reads 0,
  below any cutoff. Fewer than 50 of the 409 single-copy ones move.
- **Q4.** M does not fix the PIP3 switch: SRC knockout's AKT pS473 (0.21) stays DOWN; most of the 153 PIP3 DOWN calls
  measured unchanged stay wrong.
- **Q5.** The CV diagnostic picks 1.25 or above.
- **Q6.** Threshold-free AUROC for DOWN is above 0.5 but well below 1 (the ranking is weak too, not only the cutoff).
- **Q7.** CPPA training secondary does not fall; CPPA primary identical.
