# Positive feedback loops sit on a knife-edge: loop gain is exactly 1 at baseline

**Created**: 2026-09-18
**Status**: ceiling ε_hi=0.95 gives the first clean loop result — held-out +29 (33/4) over 12 genes, held-out macro-F1 0.8172→0.8191, relabelling churn 73→2 — while the tuning ten lose 61 and the all-pathway headline is flat (83.22%). A correctness win and a held-out win, not a headline win. Combo (ε_lo=0.3): held-out +78 but churn 21 — not adoptable. Default decision left open (§ 'Where the loop program lands').
**Flag**: `DS_LOOP_ELASTICITY` (default `1.0` = off, byte-identical)

## 1. The defect, from the algebra

Everything in DeltaSignal is a fold-change relative to baseline, and AND is
multiplication of folds. So every activator edge has elasticity
`d log(out) / d log(in) = 1`, and the gain around any positive cycle is
`1 × 1 × … × 1 = 1` at baseline.

Take the simplest AND-closed positive loop, `U → A`, `A → r1 → B`, `B → r2 → A`:

```
A = U · B,   B = A     ⟹     A = U · A
```

At `U = 1` every value of A is a solution. Baseline is a fixed point, but so is
everything else — it is **marginal**, not stable. Push `U` to 1.01 and each
pass round the loop multiplies A by ~1.01; push it to 0.99 and A decays. The
loop rails to saturation or collapses to zero, and *which* depends only on the
sign of the leak and how long you iterate.

## 2. Confirmed on the current solver

Fixed-point solver (`DS_SCC_METHOD=fixed_point`, damping 0.5), A reported in
UI, baseline = 1.0:

| U | 50 sweeps | 500 sweeps | 5000 sweeps | converged @500 |
|---|---|---|---|---|
| 1.00 | 1.0 | 1.0 | 1.0 | yes |
| **1.01** | 1.09 | 2.30 | **100.0** | no |
| **0.99** | 0.92 | 0.43 | **0.030** | no |
| 1.10 | 2.31 | 100.0 | 100.0 | yes |
| **0.90** | 0.40 | **0.0026** | 0.0026 | **yes** |
| 80.0 | 100.0 | 100.0 | 100.0 | yes |
| 0.0 | ~0 | ~0 | ~0 | yes |

A 1% perturbation rails the loop given enough sweeps. A 10% drop collapses it
to 0.3% of baseline and **reports converged**. The answer to a small
perturbation is a function of the iteration budget.

## 3. This one fact accounts for every solver anomaly measured on 2026-09-18

| observation (specs/003, specs/013) | why the knife-edge produces it |
|---|---|
| node UUID relabelling moves predictions, 0.47pp spread | on a marginal loop, sweep order decides which tiny drift wins |
| 20× iterations makes divergence *worse* (14 → 25 cases) | more sweeps = more time for `1.0001ⁱ` to reach a rail |
| all-zero is a global minimum of `‖F(x)−x‖²` | `A = U · 0 = 0` is exact for every U |
| false change is 48% of all errors, concentrated in cyclic pathways | loops railing from leaks that reach them from unrelated perturbations |
| the L_SS minimiser lost 72-to-1 at γ=1e-6 | it finds the rails *exactly*; Gauss-Seidel from baseline stops partway, so the 500-sweep cap is accidental regularisation |
| the γ prior loses at every value | it fights the knife-edge globally and drifts everything off baseline instead |
| `hill_log` beat `hill_sat` by 16 held-out at specs/010 time, "unexplained" | tanh-in-log-fold has elasticity < 1, which stabilises loops — while also compressing every acyclic cascade. **Caveat:** on the current stack the sign has flipped (`hill_log` is −26 held-out), so this row is a plausible account of a *past* measurement, not a live effect. |

The last row is suggestive, not evidence: compression may have helped because
of loops and paid for it everywhere else, but the number it explains no longer
holds on today's baseline. Prediction 3 below records the correction.

## 4. The fix: elasticity on cycle-closing edges only

Real positive feedback has gain **< 1 at rest** (the resting state is stable,
noise is damped) and **> 1 only near a threshold** (a real signal flips the
switch). That is a sigmoidal elasticity — and it is the design doc's Section
2.1 per-input sensitivity `α(x)`, which shipped with `s = 0` and has been the
identity ever since.

`DS_LOOP_ELASTICITY = ε` reads the input fold through `fold^ε` on activator
edges whose source and target share an SCC (`activator_in_loop`, set from the
existing Tarjan `comp_id` in `index_reactions`). Nothing else changes. For the
fixture (four in-loop edges):

```
A = U · A^(ε⁴)    ⟹    A = U^(1 / (1 − ε⁴))
```

a **unique, stable positive root**. Zero is still a root but now unstable
(`A^ε ≫ A` near zero pulls the loop back up), so collapse needs a real in-loop
knockout. Same fixture, same solver, `ε = 0.5`:

| U | 50 | 500 | 5000 | analytic `U^(16/15)` | converged |
|---|---|---|---|---|---|
| 1.01 | 1.0106 | 1.0106 | 1.0106 | 1.0106 | yes |
| 0.99 | 0.9894 | 0.9894 | 0.9894 | 0.9893 | yes |
| 1.10 | 1.1069 | 1.1069 | 1.1069 | 1.1069 | yes |
| 0.90 | 0.8938 | 0.8938 | 0.8938 | 0.8937 | yes |
| 80.0 | 100.0 | 100.0 | 100.0 | rail | yes |
| 0.0 | ~0 | ~0 | ~0 | 0 | yes |

Budget-independent, matches the derivation to four figures, a strong drive
still rails, an in-loop knockout still collapses. Acyclic networks are
bit-identical at every ε (no in-loop edges exist to act on).

Pinned by `test/test_loop_elasticity.jl` (66 assertions + 1 `@test_broken`
holding the knife-edge at ε = 1 so a change there announces itself).

## 5. Catalog A/B — predictions stated before the result

Running: `DS_LOOP_ELASTICITY ∈ {0.5, 0.8}` against the fixed-point baseline on
`cat_os`, wide curator set, held-out split, both concentration columns.

Predicted in advance, so they can fail:

1. **Held-out net positive**, distributed over many pathways and readouts (not
   a single-readout cluster).
2. **False change drops** in the loop-heavy pathways specifically
   (TP53, DSB repair, PIP3, WNT).
3. ~~`hill_sat + ε` recovers `hill_log`'s 16-case held-out edge.~~
   **WITHDRAWN before results, as ill-posed.** The +16 came from specs/010,
   measured on an older catalog and config. On the current baseline
   (`ab_onesided`: cat_os, depletion floor, hill_sat eps 1e-9) `hill_log` is
   **−26 held-out vs `hill_sat`** (16,084 vs 16,110 of 18,808). There is no
   edge to recover. The testable remainder: `hill_sat + ε` must not lose
   held-out accuracy the way global compression does — i.e. `hill_sat + ε ≥
   hill_sat` on held-out while `hill_log < hill_sat`. Recorded here so the
   prediction is not quietly rewritten after the fact.
4. **Label-dependence shrinks**: relabelling the catalog should move far fewer
   than 14 cases, because there is no marginal drift left for sweep order to
   decide.

### Results (constant ε)

| | ε=0.5 | ε=0.8 | baseline |
|---|---|---|---|
| accuracy | 82.22% | 83.10% | 83.36% |
| macro-F1 | 0.7762 | 0.7966 | 0.7998 |
| held-out net (fixed/broke, p) | −51 (368/419, 0.07) | **+4** (81/77, 0.81) | — |
| held-out pathways / readouts moved | 26 / 161 | 10 / 60 | — |
| false change, loop-heavy | 637 → **177** | 637 → 546 | 637 |
| false change, other | 800 → 627 | 800 → **802** | 800 |
| hill_sat+ε vs hill_sat, held-out | −51 | +4 | — |

**Prediction 1 — FAILED.** No held-out gain at either ε (−51 n.s.; +4 n.s.).
**Prediction 2 — CONFIRMED.** False change falls, and at ε=0.8 the entire −89
is in the loop-heavy pathways (+2 elsewhere): the effect lands exactly where
the mechanism says it must.
**Prediction 3 (remainder) —** fails at 0.5, passes trivially at 0.8.
**Prediction 4 — CONFIRMED, exactly.** Same network, every UUID relabelled
(verified isomorphic, 92/92), both arms at ε=0.5: **0 of 23,908 predictions
differ**; accuracy identical to the digit (82.22% / 82.22%). Under the
knife-edge (ε=1) the same relabelling moved 14 (all TP53). The solve is now a
function of the graph. This is also the prerequisite the parameter-learning
direction needs: a stable, label-independent `d(prediction)/dθ` (specs/003
research §7, [long-term vision]).

**Why accuracy did not follow:** it is a near one-for-one trade of false
change for **missed change**. Of the broken cases, 98% (ε=0.5) and 87%
(ε=0.8) are "truth is a change, arm now says NORMAL", and the flattened
signals are the small genuine ones — median fixed-point value **1.56 UI** on
the missed cases. TP53 alone is −193 / −73.

**The sharper reason:** "source and target share an SCC" is almost every edge
inside TP53's 836-node component (the giant SCCs are Type II recycling
artifacts, 1–2 reactions exploded by variants — specs/008). A constant ε
therefore compresses *every path through the component* geometrically, not
just the feedback closure. The leaks it should damp sit at ≈1.0×; the signals
it must not damp sit at ≥1.5×. They are separable by magnitude.

That is the case for the design doc's **sigmoidal** elasticity: ε<1 in a
band around baseline, →1 outside it, so drift decays and signal passes. It is
also the first thing to *learn* rather than hand-set — whether a 1.3× signal
through a loop should propagate or be damped is a biological threshold, and
the answer is in the data. Implemented next as `DS_LOOP_ELASTICITY_WIDTH`
(default 0 = constant ε, byte-identical).

## 5b. Sigmoidal elasticity — implemented, fixture-confirmed, catalog pending

`DS_LOOP_ELASTICITY_WIDTH = w` (default 0 = constant, byte-identical) makes the
in-loop elasticity a function of the input's fold:

```
ε(f) = ε_lo + (1 − ε_lo) · tanh(|log f| / w)
```

Inside the band the loop is stable; outside it the edge reads `fold^1` again
and a positive loop becomes a genuine **switch** — above-threshold signal
rails, below-threshold leak decays. Fixture (`ε_lo = 0.5`, A in UI, 2000
sweeps):

| U | constant ε | w=0.1 | w=0.3 | w=1.0 | ε=1 (knife-edge) |
|---|---|---|---|---|---|
| 1.01 | 1.0106 | 1.0107 | 1.0106 | 1.0106 | 27.8 |
| 1.10 | 1.107 | **100.0** | 1.116 | 1.109 | 100.0 |
| 1.50 | 1.541 | 100.0 | **100.0** | 1.632 | 100.0 |
| 2.00 | 2.095 | 100.0 | 100.0 | **100.0** | 100.0 |
| 0.90 | 0.894 | 0.0026 | 0.884 | 0.892 | 0.0026 |
| 0.50 | 0.477 | 0.0003 | 0.0003 | 0.0002 | 0.0003 |

Leaks (±1%) are held at every width and budget-independent (50 vs 5000 sweeps
identical). The switch threshold rises with w: w=0.1 flips at 10%, w=0.3
between 1.1× and 1.5×, w=1.0 between 1.5× and 2×. Width 0 is bit-identical to
the constant case. Acyclic networks untouched. `test/test_loop_elasticity.jl`:
115 assertions + 1 `@test_broken`.

**Choosing w from the benchmark's own band, not by taste.** The classifier is
`DOWN < 0.85`, `UP ≥ 1.15` (`bench/benchmark_vs_mpbiopath.py`), a NORMAL
half-width of ~0.15 in log-fold. For a 3-class score, railing a real
above-band signal is *correct* — only direction counts — so the requirement
is exactly: damp inside the band, permit flipping outside it. w=0.3 puts the
band edge at tanh(0.47) ≈ 0.44 (strongly stabilised) with 1.5× already in the
switch regime; w=0.15 moves the switch to the band edge.

### Predictions for the sigmoid arms (`ε_lo=0.5`, `w ∈ {0.15, 0.3}`), stated before results

5. **False change in loop-heavy pathways stays far below the baseline 637** —
   the leaks that constant ε=0.5 removed (→177) are ≈1.0× and remain inside
   the band.
6. **Missed change is mostly recovered**: the ≥1.5× genuine signals that
   constant ε flattened (median 1.56×) now pass or flip. Held-out `fixed`
   should remain in the hundreds while `broke` falls well below 419.
7. Therefore **held-out net positive** at w=0.3, distributed over many
   pathways and readouts.
8. w=0.15 carries the risk that leaks living at 1.10–1.15× cross the switch
   and rail, partially restoring false change; if so, w=0.15 shows higher
   false change than w=0.3 in the loop-heavy pathways.

### Results (ceiling, ε_lo=0.5, w=0.2)

| | ε_hi=1.0 | **ε_hi=0.95** | ε_hi=0.9 | baseline |
|---|---|---|---|---|
| accuracy | 83.75% | 83.22% | 82.98% | 83.36% |
| macro-F1 (all) | 0.8049 | 0.7987 | 0.7955 | 0.7998 |
| **held-out net** (fixed/broke, p) | +20 (22/2) | **+29 (33/4, <0.0001)** | −7 (33/40, 0.48) | |
| held-out pathways / readouts / genes | 7 / 22 / 8 | **7 / 32 / 12**, no dominance warning | 9 / 44 / — | |
| held-out accuracy / macro-F1 | | **85.66 → 85.81% / 0.8172 → 0.8191** | | |
| tuning net | +74 | **−61** (111/172) | −84 | |
| **relabelling churn** (same network, UUIDs renamed) | **73** | **2** | **1** | 14 |

**P9 — CONFIRMED.** Churn 73 → 2 → 1. The ceiling removes label-dependence
essentially completely, as the algebra said it would.
**P10 — CONFIRMED at 0.95, FAILED at 0.9.** Held-out +29 at 0.95 — larger
than the churny +20, and for the first time *distributed*: 12 genes, no single
perturbation over half. At 0.9 the gain is gone (−7 n.s.), not "modestly
reduced".
**P11 — consistent:** the 0.9 loss is tuning-side missed change.

**The honest reading.** With churn gone, hi=1.0's tuning +74 is revealed as
mostly label noise; the true tuning effect of the sigmoid is **−61** (TP53).
So the clean configuration wins on the 70 held-out pathways by every measure
and loses on the tuning ten, and the all-pathway headline is flat to slightly
down (83.22%, macro-F1 0.7987). By the project's protocol the held-out figure
is the report — but this is not a headline-accuracy improvement and is not
recorded as one. Its unambiguous wins are correctness properties: the solve is
a function of the graph (churn 2), and a loop's answer no longer depends on
the iteration budget.

### Combined configuration ε_lo=0.3 / w=0.2 / ε_hi=0.95

| | value | baseline |
|---|---|---|
| accuracy (all) | 83.46% | 83.36% |
| **held-out net** (fixed/broke, p) | **+78 (113/35, <0.0001)** | |
| held-out pathways / readouts / **genes** | **12 / 66 / 31**, no dominance warning | |
| held-out accuracy / macro-F1 | **85.66 → 86.07% / 0.8172 → 0.8214** | |
| tuning net | −53 (128/181) | |
| **relabelling churn** | **21 (12 held-out)** — fails the single-digit criterion; ε_lo=0.5's ceiling config has 2 | 14 |

The best held-out result of the loop program, and distributed by gene (31, no
perturbation over half) — but **not by pathway in net terms**: DSB repair is
**+74 of the +78**; the other eleven moved pathways roughly cancel (+13 PTK6,
+6, +6, +5, +5, +3 against −18 MHC, −10 RUNX2, −3). The concentration columns
count pathways *moved*, which hid this; the per-pathway net is the honest view
and is now part of how loop results are read.

**A pattern across every loop configuration:** TP53 loses (−80 here) and DSB
repair gains (+74), and both are giant-SCC pathways (836 and 1,127 nodes). The
same global ε helps one loop and hurts the other. That is the concrete case
for *learning* ε per loop or per edge type rather than setting one constant —
which is the design's intent (Section 2.1) and the parameter-learning
direction. The tuning ten lose again; every loop configuration has, and it is TP53
(the 836-node component) each time. Recorded as: held-out +78 and headline +0.10pp, tuning −53, **and churn 21** —
so the +78 carries label noise of the same order as its own concentration in
DSB repair. Not adoptable as it stands.

### Where the loop program lands (2026-09-19)

| configuration | held-out net | churn | headline | tuning |
|---|---|---|---|---|
| original (knife-edge) | — | 14 | 83.36% | — |
| ε_lo=0.5, w=0.2, ε_hi=1.0 | +20 | 73 | 83.75% | +74 (noise) |
| **ε_lo=0.5, w=0.2, ε_hi=0.95** | **+29 (33/4), 12 genes** | **2** | 83.22% | −61 |
| ε_lo=0.3, w=0.2, ε_hi=0.95 | +78 (113/35), 31 genes, DSB +74 | 21 | 83.46% | −53 |

The only configuration that is both held-out-positive and reproducible under
relabelling is **ε_lo=0.5 / w=0.2 / ε_hi=0.95**. It is a correctness win (the
solve is a function of the graph; a loop's answer no longer depends on the
iteration budget) and a held-out win (+29 over 12 genes, macro-F1 +0.0019),
and it is **not** a headline win (−0.14pp overall, tuning −61). Whether to
make it the default is a judgement about which of those to weight, and is
recorded here as an open decision rather than made. The per-loop pattern
(TP53 loses, DSB gains under every ε) says the durable answer is learned
per-loop parameters, not a better constant.

### Results (sigmoid, ε_lo = 0.5)

| | w=0.15 | w=0.3 | baseline |
|---|---|---|---|
| accuracy | **83.70%** | 83.34% | 83.36% |
| macro-F1 | **0.8043** | 0.7991 | 0.7998 |
| ALL net (fixed/broke, p) | **+82** (85 / **3**, <0.0001) | −4 (128/132, 0.85) | |
| held-out net | **+19** (19 / **0**, <0.0001) | +9 (32/23, 0.28) | |
| held-out pathways / readouts / **genes** | 4 / 18 / **5** | 10 / 34 / — | |
| false change, loop-heavy | 616 | 602 | 637 |
| vs hill_log, held-out | +45 | — | −26 |

**P5 — FAILED.** False change moved −49 / −75, nowhere near constant-ε's
−460. The sigmoid does not damp the leaks constant ε damped; whatever it
fixes, it is not by that route.
**P6 — CONFIRMED.** Held-out broke 0 and 23, against constant-ε's 419.
**P7 — not established as stated** (w=0.3 is +9 n.s.). The positive result is
at **w=0.15**, which P8 had flagged as the risky width.
**P8 — directionally true, trivially** (loop-heavy false change 616 vs 602)
and its implied conclusion — that w=0.3 is safer — was wrong.

**What w=0.15 actually fixed** (85 cases, 3 broke):

| transition (truth) | n | reading |
|---|---|---|
| UP → DOWN (DOWN) | **33** | knife-edge railed the loop *up* on a leak and drowned a real downward signal |
| UP → NORM (NORM) | 31 | false change removed |
| DOWN → NORM (NORM) | 19 | false change removed |
| DOWN → UP (UP) | 2 | |

Fixed-point median 1.53× → sigmoid 0.93× on these cases: they had been pushed
just past the 1.15 UP cutoff by loop amplification of a leak. The mechanism
that pays is **wrong-sign correction**, not damping — stabilising the loop
stops amplification from overriding a genuine signal.

**Concentration, stated plainly.** The 19 held-out fixes are **13 from one
gene (DOK1, PTK6) + 3 from one gene (STAT3, MET) + 3 singletons**. 18 distinct
readouts, but effectively two perturbations. The tuning side has the same
shape: **62% of its 69 discordant cases are one perturbation, BRCA1 in TP53**
(the report's new dominant-gene share surfaced this; a plain gene count did
not — TUNING shows "6 genes"). `p < 0.0001` on 19/0 is therefore
overstated: the gain is real and nearly free (3 broke catalog-wide), but the
*independent* held-out evidence is thin. The tuning +54 is TP53 and is
subject to the relabelling check (running). This result is the reason the
report now counts distinct genes beside distinct readouts.

**Width sweep, judged on TUNING** (ε_lo=0.5):

| w | acc | macro-F1 | ALL (fixed/broke) | held-out | TUNING |
|---|---|---|---|---|---|
| 0.1 | 83.45% | 0.8005 | +23 (27/4) | +18 (18/0) | +5 n.s. |
| 0.15 | 83.70% | 0.8043 | +82 (85/3) | +19 (19/0) | +63 |
| **0.2** | **83.75%** | **0.8049** | **+94 (99/5)** | **+20 (22/2)** | **+74** |
| 0.3 | 83.34% | 0.7991 | −4 | +9 n.s. | −13 n.s. |

TUNING chooses **w=0.2**; held-out reports **+20 (22/2)** over 7 pathways / 22
readouts. Baseline 83.36% / 0.7998.

### Label-dependence RETURNS above the band — worse than the knife-edge

Same network, every UUID relabelled, both arms identical config:

| config | predictions that differ | held-out |
|---|---|---|
| ε=1 (original solver) | 14 | 0 |
| constant ε=0.5 | **0** | 0 |
| sigmoid ε_lo=0.5, w=0.15 | **73** | **20** |

The relabelled catalog at w=0.15 scored **83.89% vs 83.70%** — a 0.19pp swing
from node names alone, the same order as the held-out gain. Cause, from the
algebra: above the band the sigmoid returns ε to 1, restoring loop gain 1 and
therefore the knife-edge for every above-band signal — now behind a sharp
boundary that sweep order can tip a case across. Constant ε had no such
boundary and no churn.

**This caps how much of the +82/+94 can be believed.** The sign is consistent
(both labellings beat baseline), so the gain is very likely real, but its
magnitude is uncertain by ~0.2pp and the pre-registered "distributed, p<0.0001"
framing is overstated twice over — by dominance (DOK1 68%, BRCA1 62%) and by
churn.

**Fix under test: a ceiling.** `DS_LOOP_ELASTICITY_HI = ε_hi < 1` (default
1.0, byte-identical) so the sigmoid rises to ε_hi rather than to 1. Loop gain
then stays below 1 at every amplitude — unique root everywhere, so churn should
return to 0 as it was for constant ε — while a real signal still clears the
classifier: 1.5× through ten in-loop edges at ε_hi=0.9 reads ≈1.19 > 1.15.
It is the third loop parameter (ε_lo, w, ε_hi), all three learnable.

**Fixture, ε_lo=0.5, w=0.2** (A in UI; `*` = not converged at 2000 sweeps):

| U | ε_hi=1.0 | ε_hi=0.95 | ε_hi=0.9 | analytic `U^(1/(1−ε_hi⁴))` at 0.9 |
|---|---|---|---|---|
| 1.01 | 1.0106 | 1.0106 | 1.0106 | in band |
| 1.5 | **100.0** | 8.90 | **3.25** | 3.25 |
| 2.0 | 100.0 | 42.0 | 7.50 | 7.50 |
| 3.0 | 100.0 | 100.0 | 24.4 | 24.4 |
| 10, 80 | 100.0 | 100.0 | 100.0 | rail |
| 0.9 | 0.859 | 0.874 | 0.880 | in band |
| 0.0 | ~0 | ~0 | ~0 | 0 |

With the ceiling at 1.0 every above-band signal rails (the knife-edge). With
ε_hi<1 it settles at the analytic root, matched to three figures, and
converges to one value (500 vs 5000 sweeps agree). Direction is preserved,
which is all a 3-class score needs; the amplification exponent `1/(1−ε_hi⁴)`
is the biologically meaningful "loop gain" and the natural thing to learn.
Pinned in `test/test_loop_elasticity.jl`.

### Predictions for the ε_hi arms, stated before results

9. **Relabelling churn returns to 0** (or single digits) at ε_hi ∈ {0.9, 0.95}
   with w=0.2, because no amplitude has loop gain 1.
10. **The held-out gain survives** at ε_hi=0.95 (≥ +15) and is at most
    modestly reduced at 0.9, because the above-band compression is mild
    relative to the 1.15 cutoff.
11. If (10) fails at 0.9 but holds at 0.95, the loss is missed change on long
    in-loop paths, and that is where a *learned* ε_hi earns its keep.

## 6. Where this goes

`ε` is the first **learnable loop parameter**, and it is per-edge-type by
construction (in-loop vs not). The full sigmoidal `α(x)` — low elasticity at
baseline, rising near a threshold — is the path to genuine bistability with
hysteresis, and it is already implemented and inert
(`activator_sensitivity_s = 0`). Negative feedback needs nothing extra: with
divide-form inhibition the loop attenuates (`A = U^(1/(1+ε²))`) as biology
expects.

Non-goals: no change to acyclic evaluation (it is exact and stays exact); no
edge deletion; no entity-level loop classification (specs/008 records why).
