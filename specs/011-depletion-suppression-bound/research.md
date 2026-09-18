# Bounding depletion suppression

**Created**: 2026-09-18
**Status**: Landed
**Related**: `specs/010-and-multiplication-fidelity` (same asymmetry, different
operator); supersedes nothing

## Why this exists

Depletion inhibition was capped above and unbounded below:

```julia
H_dep = clamp(H_dep, zero(T), h_max_dep)    # h_max_dep = 10
```

De-repression could raise a node at most 10x. Suppression could drive it to
**zero**, and multiple depletion edges compound multiplicatively. The model was
therefore free to assert that an abundant complex depletes its own free subunit
without limit. That is not mass conservation — it is an unbounded reciprocal,
and it is the same capped-above / unbounded-below asymmetry `specs/010` found
in the AND operators.

## How it was found — one case, not an aggregate

This matters methodologically and is worth recording, because nine consecutive
levers in this project died as between-pathway confounds while this took under
an hour.

1. Root-caused **all 2,519** cases that MP-BioPath's networks answer and ours
   do not, same propagator and config so only the network differs. Labels are
   mechanical facts about each case, and were checked for internal consistency
   (seven checks, zero violations) before being believed.
2. That gave **431 wrong-SIGN** cases.
3. Split them by path parity: **174** have an odd-negative path in OUR network
   and lose the inversion anyway; 104 have an all-positive route theirs lacks;
   **94 have no inversion in EITHER network yet MP-BioPath is right**, which
   may not be our deficit at all given their solver's known direction-inversion
   defect.
4. Traced one of the 174 node by node: `EPS15` knockout → `p21 RAS:GTP` in
   `Signaling_by_EGFR`, truth UP, ours DOWN.

The inversion **worked** — the negative depletion edge de-repressed correctly
and the signal reached **70x** baseline at the EGFR dimer. It then collapsed to
**0.007** at one reaction whose two AND inputs were 70x and **1e-4**; 70 × 1e-4
is exactly the output, so the AND was right too. Walking the 1e-4 backwards
reached free `GRB2-1 [cytosol]`, whose only inputs are **two negative depletion
edges** from EGFR:CBL complexes the knockout drove to 100x. Each contributed
`bl/x = 1/100`; compounded, 1e-4. With the floor the case reads **8.5 → UP**.

## The fix

`clamp(H_dep, T(config.depletion_h_min), h_max_dep)`, with
`DS_DEPLETION_H_MIN` defaulting to `1/DS_DEPLETION_H_MAX` = 0.1 — symmetric in
log space: depletion may suppress at most as hard as it may de-repress.

Made configurable deliberately. "Symmetric with `h_max`" is principled only
relative to `h_max = 10`, which is itself a tuned constant, so the value had to
be settled by measurement rather than by the symmetry argument.

## Measurement

93-pathway catalog regenerated from LNG main and gated on its own
`src_sha256`; per-arm env verified by `docker exec` before any number was
taken; conditioned pairing; 0 of 23,908 cases dropped. 81 pathways scored — 11
tuning (5,100 cases), 70 held out (18,808).

| split | net | fixed / broke | McNemar p |
|---|---|---|---|
| tuning | +99 | 114 / 15 | < 0.0001 |
| **held-out** | **+28** | **35 / 7** | **< 0.0001** |

macro-F1 **0.7835 → 0.7929**, accuracy 82.42% → 82.95%. This also beats the
`hill_log` default that preceded `specs/010` (0.7846). Positive and significant
on **both** splits, so it is not overfitting.

### The floor value, swept on the TUNING split only

| floor | tuning macro-F1 | TP53 | PIP3 | held-out net | held-out p |
|---|---|---|---|---|---|
| 0.0 (original) | 0.6945 | 1349 | 420 | — | — |
| 0.001 | 0.6940 | 1347 | 420 | not run | |
| 0.01 | 0.7047 | 1387 | **420** | +4 | **0.1250** |
| **0.1 (symmetric)** | **0.7199** | **1459** | 406 | **+28** | **< 0.0001** |

**0.01 is the Pareto-safe option and was rejected on evidence.** It gains 38
TP53 cases with PIP3 completely unharmed, which looks obviously correct — but
its held-out gain is **4 discordant pairs at p = 0.125**, not established. It
buys PIP3 protection by surrendering the only out-of-sample effect the change
has. The Pareto-safe option is not automatically the conservative one: choosing
it would have shipped a change whose only measured effect was on the tuning
set, which is precisely what the held-out protocol exists to catch.

## What it costs, stated plainly

**87% of the gain is one pathway.** TP53 is +110 of the +127 net, and TP53 is
a tuning pathway. That fully explains the tuning-vs-held-out asymmetry
(+1.94pp per case tuning against +0.15pp held out).

**It significantly degrades the best pathway.** `PIP3_activates_AKT_signaling`
goes 420/448 → 406/448, i.e. **93.75% → 90.63%, p = 0.0001**, every loss a
`truth DOWN, was DOWN, now UP`. PIP3 carries the "DeltaSignal matches
MP-BioPath" claim, so this is a real cost and is recorded rather than left to
be discovered.

The held-out +28 is genuinely distributed — EGFR +10, NOTCH1 +7, NOTCH2 +5,
ROBO +3 — and is the out-of-sample evidence, independent of TP53.

**Experimental axis: no detectable effect.** Net −3, 8 fixed / 11 broke,
McNemar p = 0.65. Same directional pattern as curator (TP53 +5, PIP3 −9), too
small to call at n=849. Every experimentally-evidenced pathway is inside the
tuning ten, so there is no held-out empirical split for either tool.

## The dominant mechanism is broader than the traced case

102 of the 149 fixed cases are `truth DOWN, was UP, now DOWN` — not the
UP-recovery the trace found. Unbounded depletion was mostly causing **runaway
de-repression**: suppressing inhibitors toward zero maximally de-repressed
their targets into spurious UP calls. Bounding it keeps them nearer baseline.
The EPS15 trace is real but not representative of the bulk.

## Relation to two earlier negative depletion results

`substrate depletion` measured −14pp and `general consumption` measured
macro-F1 0.663 → 0.529. Both **added** depletion edges. This **bounds** ones
already present. Recorded so the earlier results are not read as contradicting
this one.

## Open

The same bucket still holds **173 more odd/odd wrong-sign cases**, 104 where
our network offers a positive route theirs lacks, and 94 where neither network
has an inversion. Each is a separate trace, and the method now has a record.

Closing the full 2,519-case gap while keeping the 573 cases our networks
already win would land at **95.88%** on the 22,106 cases both network sets
score — so a 95% target sits below an already-scoped ceiling, with "both
wrong" at 4.12% as the practical floor.
