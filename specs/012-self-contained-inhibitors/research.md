# Structural double-counts: two real mechanisms, both refuted as edge filters

**Created**: 2026-09-18
**Status**: NEGATIVE RESULT — mechanism confirmed, the obvious fix refuted
**Flags**: `DS_SKIP_SELF_INH`, `DS_DEDUP_ACTIVATORS` (benchmark-side diagnostics, both default off)

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

---

# Second mechanism: an entity counted twice as its own activator

## The mechanism, which is also real

An entity that is both the **catalyst** and a **substrate** of the same
reaction contributes TWO activator edges, so the AND product squares its
fold-change. A 5x input yields 25x.

Traced from the false-change bucket — `TP53` knockout in
`Cell_Cycle_Checkpoints`, truth NORMAL, we said UP:

| step | node | fold |
|---|---|---|
| 0 | p-S166,S188-MDM2 dimer | 0.00 (knocked out) |
| 2–6 | de-repressed downstream complexes | 5.00 |
| **7** | **MDM2 ubiquitinates phosphorylated MDM4** | **24.9995** |
| 8 | readout | 24.9995 → UP |

That reaction's activator edges:

| source | role | fold |
|---|---|---|
| Ub | input | 1.00 |
| `R-HSA-6804936` MDM2:MDM4 | **input** | **5.00** |
| `R-HSA-6804936` MDM2:MDM4 | **catalyst** | **5.00** |

5 × 5 = 25. The same entity, counted once per role it plays.

**Prevalence: 7,013 reactions — 15.8% of all reactions with two or more
activator edges — across 78 pathways.** The repeated pair is `catalyst+input`
in 7,042 of 7,045 cases.

## The fix that does not work, now established

Collapsing duplicate activator edges to one per (source, reaction):

| split | net | fixed / broke | McNemar p |
|---|---|---|---|
| tuning | −4 | 20 / 24 | 0.65 |
| **held-out** | **−15** | **3 / 18** | **0.0015** |

accuracy 82.95% → 82.87%.

A "dedup activators" arm was recorded as negative once before, at **p = 0.25 on
the 742-case set** — which is noise, not a refutation, and that set has reversed
decisions in opposite directions twice. On the wide curator set the result is
**significant**: 3 fixed against 18 broken.

**One honest tension:** macro-F1, the project's stated primary metric, moves the
other way — 0.7929 → 0.7933, with DOWN-F1 0.734 → 0.747. The paired held-out
test at p = 0.0015 is much stronger evidence than a +0.0004 macro-F1 move, so
this is recorded as refuted, but the disagreement is real and is why both
numbers are written down.

## What the two results have in common

Both mechanisms are real, verifiable, and affect thousands of reactions. Both
"obvious" corrections lose significantly on held-out cases. Together with the
AND-compression result in `specs/010` (correct arithmetic costs 16 held-out
cases), three of the four structural defects characterised on 2026-09-18 turned
out to be **load-bearing**: the model leans on them, and removing them alone
makes predictions worse.

The one that won — bounding depletion suppression, `specs/011` — did not remove
a mechanism. It **bounded** one that was unbounded in a single direction.

That is the pattern worth carrying forward: on these networks, *constraining* a
runaway operator has paid off, and *deleting* a structurally-wrong term has not.

## Addendum 2026-09-25: two self-contained inhibitors invert the signal

specs/021 traced a dose-response reversal (WNT5A, Signaling by WNT) to a
reaction with **two** self-contained inhibitors. One such inhibitor cancels
(x · 1/x = 1), as recorded above. Two give x · min(x^-2, 10): a knockdown reads
UP until the de-repression ceiling binds, then DOWN. 94 reactions catalog-wide
have two or more (PIP3 61, WNT 16). The full trace and the untested bounded
alternative (combine them by min) are in `specs/021-empirical-holdout-axis/research.md`.
