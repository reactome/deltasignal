# Phase 0 Research: Long-range upregulation propagation

**Plan**: [plan.md](./plan.md) | **Date**: 2026-09-10

All figures below come from a single reproducible catalog build (`cat`), ten
MP-BioPath evaluation pathways, Reactome Release97, LNG `main` @ 4ff0408,
root-input perturbation, harness defaults.

## R1 — SC0: is the measurement stable? YES. The gate was raised in error.

Two independent catalog builds with `PYTHONHASHSEED=0`, same arm:

| | build 1 | build 2 |
|---|---|---|
| all scored | 334/564 | 335/564 |
| both-converged subset | 233/432 | 233/432 |
| predictions differing | — | 6, all TP53, **0** both-converged |

All ten networks hash identically across builds (canonical form, uuid4
relabelled), so generation is deterministic. Build-to-build noise on the
headline is **±1 case**, and the converged subset is bit-identical.

**The blocking gate in plan.md was wrong and is withdrawn.** It rested on a
334-vs-295 comparison against a catalog destroyed in a machine crash. There
are zero generator changes between the two LNG versions involved (the diff is
logging and credential redaction only), so the discrepancy cannot be
attributed to code, and the old catalog cannot be rebuilt to find out. **295
is not a defensible figure and is withdrawn.**

The lesson is specific: uuid noise was invoked to explain a 39-case gap having
previously been *measured* at ±3–4 cases. A mechanism's known magnitude is
part of the evidence for or against it, and 39 should have counted against
that explanation rather than for it.

**Standing rule adopted:** an A/B is only valid across arms sharing one
catalog build. Comparing across builds is not supported even though the
networks are identical, because the non-converged tail is not.

## R2 — The assembly clamp: real, and one third the size I reported

Same catalog build, only `DS_ASSEMBLY_LIMITING` differs:

| | all 564 | macro-F1 | UP 246 | DOWN 247 |
|---|---|---|---|---|
| clamp ON (default) | 334 | 0.5601 | 119 | 169 |
| clamp OFF | **355** | 0.5653 | **161** | 166 |
| **delta** | **+21** | +0.005 | **+42** | **−3** |

Both-converged subset (382 cases): 207 → 228, the same +21.

**Decision**: the mechanism is confirmed; the magnitude is revised from +73 to
+21. The earlier figure used the withdrawn 295 as its baseline.

**Rationale**: the effect is concentrated exactly where predicted — +42 on
upregulation, −3 on downregulation. `min(elevated, baseline)` is exactly
baseline, so a complex transmits scarcity perfectly and blocks abundance
completely. Structural evidence is unaffected by the baseline error: the
path-exists-but-untouched count and the direction asymmetry both stand.

**Alternatives considered**: per-hop `hill_log` attenuation was eliminated —
`hill_log` is identity to four decimals at both classification cutoffs.
Readout-local masking was eliminated — untouched readouts are 30.0%
multi-producer against 26.1% for correct ones, indistinguishable, and `max`
masks decreases rather than increases.

## R3 — The clamp is not most of the gap

| | vs MP-BioPath (407) |
|---|---|
| clamp ON (334) | −73 |
| clamp OFF (355) | −52 |

Removing the clamp closes **29%** of the gap. My earlier framing — "the single
biggest factor in the low number" — overstates it. It is the largest single
*identified* factor; the majority of the gap is still unexplained.

DeltaSignal also remains behind the structural baseline (393) with the clamp
off, so **SC2 is not met by this change alone**.

Residual failure modes for expected-UP, clamp ON, still to be explained:

| | count |
|---|---|
| no directed path | 55 (22.4%) |
| wrong direction | 42 (17.1%) |
| attenuated into the no-change band | 24 (9.8%) |

The wrong-direction group is the most interesting: 42 cases where an
overexpression produces a *decrease*. That is not attenuation and not
reachability, and no hypothesis currently covers it.

## R4 — TP53 convergence: real, but not a measurement blocker

TP53 is 232 of 564 scored cases (41%) and 122 of them (52.6%) do not converge.
Every other pathway is 0%.

This does **not** destabilise the headline (R1), because the non-converged
cases are stable given a fixed catalog build. It matters for two other
reasons: those cases are scored as though equivalent to converged ones, which
Principle V calls out; and they are the only source of cross-build variation.

**Decision**: report both the full figure and the both-converged subset;
do not gate Stage 2 on it. Prior work characterises the giant SCC as Type II
catalytic recycling rather than genuine feedback, which makes it a
network-structure question for a separate feature.

## R5 — Two departures from the design, not one; together worth +31

Adam identified the AND design intent: multiply the fold-changes with a
Hill-like curve, constrained to 0-100, and "if it is near 1 I want it to be
extremely close to pure multiplication" — 0.5*0.5=0.25, 1*1=1, 2*0.5=1, and
sig(100)*sig(100)=100 rather than 10,000.

Checked against the two implemented curves:

| UI folds | pure x (clamped) | `hill_log` (default) | `hill_sat` |
|---|---|---|---|
| 0.5 x 0.5 | 0.25 | 0.252 | 0.250 |
| 1 x 1 | 1.0 | 1.0 | 1.0 |
| 2 x 0.5 | 1.0 | 1.0 | 1.0 |
| 10 x 10 | 100 | **74.06** | 99.98 |
| **[100] alone** | **100** | **74.07** | **99.99** |
| **[50] alone** | **50** | **41.43** | **50.00** |

`hill_log` tanh-squashes the summed log-fold with `z_max=10`, and e^10 is far
outside the UI range, so the squashing is active THROUGHOUT the operating
range rather than only near the ceiling. A lone node at UI 50 reads 41.4.
`hill_sat` implements the stated intent — its own docstring says "exact
multiplication through the common operating range, with sigmoid transitions
ONLY at the UI=0 and UI=100 boundaries".

**Decision**: `DS_AND_MODE=hill_sat` with `DS_HILL_SAT_EPS=1e-5`, and
`DS_ASSEMBLY_LIMITING=0`. Measured, same catalog build:

| arm | correct | macro-F1 | UP 246 | DOWN 247 | pred UP | non-conv |
|---|---|---|---|---|---|---|
| defaults | 334 | 0.5601 | 119 | 169 | 159 | 122 |
| `hill_sat` only | 339 | 0.5671 | 122 | 171 | 163 | 132 |
| clamp OFF only | 355 | 0.5653 | 161 | 166 | 222 | 164 |
| **both** | **365** | **0.5781** | **168** | **169** | 230 | 179 |

Both-converged subset (362 cases, the comparison immune to non-convergence):
199 -> 216, MP-BioPath 253.

Verified genuine rather than a bias: median predicted UI stays at exactly
1.000, predicted-UP moves toward the true distribution (159 -> 230 against 246
expected) rather than past it, DOWN is unchanged, and macro-F1 improves.

The two effects compose slightly better than they add (+5 and +21 separately,
+31 together), consistent with the mechanism: the clamp destroyed increases
before they grew large enough for the compression to matter, so each masked
the other's cost.

## R6 — NEGATIVE, and a caught false positive: the default `hill_sat` epsilon

`DS_HILL_SAT_EPS` defaults to 0.001, which is LARGER than the internal values
where knockouts live (~0.0006), so it acts as a floor and lifts the whole
network upward. Measured: `hill_sat` at the default epsilon predicts UP on
**317 of 564** cases against 246 actually UP, with median output 1.587 rather
than 1.000. Accuracy 310, macro-F1 0.4848 — the worst arm tested.

I initially reported that arm's UP score (170/246, the best of any arm) as
evidence the curve worked. It was not: it was over-prediction. The tell was in
the same table — DOWN collapsing 169 -> 119 while UP rose. A genuine
improvement does not usually gut the other half. Macro-F1 flagged it
correctly and I quoted the number without letting it inform the conclusion.

**This is a live defect independent of this feature**: anyone selecting
`hill_sat` today gets a systematically upward-biased model.

## Open questions carried into Stage 2

1. What soft-minimum form satisfies FR2 and FR3 together? Untested.
2. Does replacing the clamp preserve the +0.9pp curator win it was adopted
   for? Must be measured against curator ground truth, not only experimental.
3. What explains the 42 wrong-direction upregulation cases?
4. Why is DeltaSignal behind a model-free traversal even with the clamp off?

## R7 — the curator axis does NOT regress; it moves the most (T011/T012)

Open question 2 above is answered, and the answer is the opposite of the
worry that motivated the task. Both arms, same shared catalog build
(`$S/cat`, ten pathways, LNG main 4ff0408, Release97), `--ground-truth
curator`, 3914 cases DeltaSignal scores:

| arm | correct | acc | macro-F1 | F1_DOWN | F1_NO_CHANGE | F1_UP |
|---|---|---|---|---|---|---|
| old defaults (`hill_log`, eps 1e-3, clamp ON) | 2479/3914 | 0.6334 | 0.5798 | 0.5949 | 0.7316 | 0.4130 |
| **new defaults** (`hill_sat`, eps 1e-5, clamp OFF) | **2628/3914** | **0.6714** | **0.6408** | 0.5814 | 0.7550 | **0.5858** |

**+149 cases, +3.81 pp, macro-F1 +0.0609.** Recall by class: UP
343/985 -> **592/985**, DOWN 522/991 -> 523/991, NO_CHANGE 1614/1938 ->
1513/1938. The gain is upregulation and nothing else, which is the
hypothesis this feature was written to test, now confirmed on the ground
truth it was *not* tuned against.

The trade-off T012 was reserved for does not exist on this catalog, so no
decision is being deferred to Adam. Two honest caveats:

- **DOWN F1 drops slightly** (0.5949 -> 0.5814) even though DOWN recall is
  flat, because the arm predicts DOWN 764 -> 808 times. Precision, not
  recall, pays for it.
- **The +0.9pp curator result that justified the clamp was measured on the
  92-pathway catalog, and its named regression (Interferon α/β −75) is not
  in these ten pathways.** So this is not a re-measurement of that exact
  experiment; it is the same axis on the benchmark catalog. Whether the
  clamp still earns its keep on Interferon α/β specifically is untested and
  stays open.

Per-pathway (curator, changed predictions only; 282 of 462 changed with both
arms converged, so this is not the uuid4 sweep-order artifact):

| pathway | + | − | net |
|---|---|---|---|
| PIP3_activates_AKT_signaling | 81 | 9 | +72 |
| Mitotic_G1-G1_S_phases | 44 | 9 | +35 |
| HDR_through_Homologous_Recombination | 41 | 21 | +20 |
| Cell_Cycle_Checkpoints | 18 | 0 | +18 |
| S_Phase | 10 | 0 | +10 |
| Mitotic_Prophase | 9 | 1 | +8 |
| Signaling_by_WNT | 26 | 19 | +7 |
| Signaling_by_ERBB2 | 6 | 6 | 0 |
| Transcriptional_Regulation_by_TP53 | 45 | 66 | **−21** |

TP53 is the only net loss and it is the pathway that does not converge.

## R8 — where this leaves us against the baselines (T013)

Experimental ground truth, identical 564 scored cases:

| model | correct | acc | macro-F1 |
|---|---|---|---|
| DeltaSignal, old defaults | 334/564 | 0.5922 | 0.5601 |
| **DeltaSignal, new defaults** | **365/564** | 0.6472 | 0.5781 |
| `shortest_signed_path` (model-free traversal) | 393/564 | 0.6968 | 0.6177 |
| MP-BioPath | 407/564 | 0.7216 | 0.6707 |

Curator ground truth, identical 3914 scored cases: DeltaSignal 2628 (0.6714),
`shortest_signed_path` 2824 (0.7215), MP-BioPath 3034 (0.7752).

So this feature closes **31 of the 73-case experimental gap to MP-BioPath**
and leaves 42. Open question 4 stands unanswered: **DeltaSignal is still
behind a model-free signed traversal on both axes.** A propagator that loses
to "follow the shortest signed path" is not yet earning its complexity, and
that — not the remaining distance to MP-BioPath — is the sharpest statement
of what is left to fix.

Per-pathway on the experimental axis, new defaults:

| pathway | DS | MP-BioPath | shortest-path |
|---|---|---|---|
| Cell_Cycle_Checkpoints | 35/44 | 42/44 | 38/44 |
| HDR_through_Homologous_Recombination | 24/32 | 28/32 | 24/32 |
| Mitotic_G1-G1_S_phases | 21/73 | 31/73 | 21/73 |
| Mitotic_Prophase | 17/20 | 18/20 | 17/20 |
| PIP3_activates_AKT_signaling | 75/84 | 79/84 | 79/84 |
| RAF_MAP_kinase_cascade | 2/4 | 4/4 | 4/4 |
| S_Phase | 10/15 | 12/15 | 10/15 |
| Signaling_by_ERBB2 | 13/23 | 16/23 | 17/23 |
| Signaling_by_WNT | **24/37** | 23/37 | 28/37 |
| Transcriptional_Regulation_by_TP53 | 144/232 | 154/232 | 155/232 |

DeltaSignal beats MP-BioPath on exactly one pathway (WNT, +1) and loses
everywhere else. `Mitotic_G1-G1_S_phases` at 21/73 is the worst absolute
deficit and is also the **only pathway this feature made worse** on the
experimental axis (−10 net, 15 of its 20 changed predictions converged in
both arms, so it is a real regression rather than noise).

## R9 — the change makes non-convergence worse, as expected (T016)

Non-converged solves on the experimental arm rise **122/564 -> 179/564**.
This is the predicted consequence, not a surprise: removing the clamp lets
larger signals reach the cyclic components, and larger signals circulate.
It is a direct argument that `specs/003-solver-objective` is now *more*
urgent — 32% of scored cases currently report a value that depends on where
the fixed-point iteration was stopped. Out of scope here; do not attempt to
fix it in this feature.
