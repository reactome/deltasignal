# Inhibitors that contain their own substrate

**Created**: 2026-09-18
**Status**: NEGATIVE RESULT — mechanism confirmed, the obvious fix refuted
**Flag**: `DS_SKIP_SELF_INH` (benchmark-side diagnostic, default off)

## The mechanism, which is real

An inhibitor whose entity **contains** one of the same reaction's activators is
not an independent variable. It is a partition of the substrate pool, and it
rises *because* the substrate rose. With divide-form inhibition
(`H = baseline / x`) the two cancel, and an elevated substrate yields exactly
baseline flux.

Traced from a real failure in the "path exists, signal never arrives" bucket.
`ATM` over-expression in `Cell_Cycle_Checkpoints`, truth UP, we said NORMAL:

| step | node | fold |
|---|---|---|
| 0 | p-S1981,Ac-K3016-ATM | 80.00 |
| 1–10 | ATM → TP53 → CDKN1A transcription → **CDKN1A mRNA** | **76.83** |
| **11** | **PCBP4 modulates CDKN1A translation** | **1.0000** |
| 12–14 | → Cyclin A:Cdk2:p21 (readout) | 1.0000 |

That reaction has exactly two inputs:

| role | entity | fold |
|---|---|---|
| activator (`pos and input`) | CDKN1A mRNA | 76.83 |
| inhibitor (`neg or regulator`) | **PCBP4:CDKN1A mRNA** | 76.83 |

`76.83 × (0.01 / 0.7683) = 0.9998`. Ten steps of correct propagation annulled
to exactly baseline, and the output lands on baseline rather than being
suppressed — the signature that made it visible.

**Prevalence: 150 reactions, 22.7% of all reactions carrying both an activator
and an inhibitor, across 37 pathways.** Measured from the containment table the
generator ships, so nothing is inferred.

## The fix that does not work

Dropping those inhibitor edges at solve time. Measured on the 93-pathway
catalog, conditioned pairing, per-arm env verified:

| split | net | fixed / broke | McNemar p |
|---|---|---|---|
| tuning | −37 | 30 / 67 | 0.0002 |
| **held-out** | **−61** | **28 / 89** | **< 0.0001** |

macro-F1 0.7929 → 0.7878; accuracy 82.95% → 82.54%.

**It fixes 28 cases and breaks 89 — roughly 3:1 against.** Those edges encode
real curator-asserted negative regulation that mostly does useful work. The
cancellation is a real defect in a minority of the cases where the pattern
appears, and blanket removal costs far more than it recovers.

This is the same shape as the four uuid-silo bridging attempts: a real
structural observation, and an intervention that over-corrects.

## What would be worth trying instead

Not removal but attribution. The inhibitor's rise has two components — the part
explained by the substrate's own rise, and any excess. Only the excess is
independent inhibition. A formulation that inhibits on the excess would leave
genuine regulation intact while removing the self-cancellation.

That is a change to the inhibition operator rather than an edge filter, so it
needs a stronger prior than this result provides. Recorded as a direction, not
a plan.

## Do not re-derive

The mechanism is real and worth knowing; the edge-removal fix is measured and
negative at p < 0.0001. `DS_SKIP_SELF_INH=1` reproduces the arm.
