# specs/144 research: the engagement dial, measured

## Verdict

**H1 fails as registered.** Balanced accuracy over non-marker scored readouts
(n = 9,476; markers excluded from both arms): control 0.4031, conditioned
0.4046, difference +0.0016 with a 10,000-resample by-drug-set bootstrap of
[−0.0116, +0.0116] — the interval includes 0. RAF set aside: +0.0056, same
conclusion. Secondary threshold: −0.0067.

## What the dial does do (the traces, before any claim)

- **Raw accuracy +5.0pp** (0.8545 → 0.9047), **495 fixed / 19 broke** — and
  471 of the 495 fixed (95%) and 19 of 19 broken touch a drug the specs/138
  table marks non- or unknown-engaged. The predicted mechanism is exactly
  what happened: non-engaged drugs stop making false change calls.
- **The fitted dials recover the engagement table without being told it.**
  Single-drug α: SR 0.004, P6 0.055, NT 0.100, AK 0.099, ST 0.145, RY 0.177,
  PLX 0.57, ZS 0.992, 901 1.0 — the specs/138 "did not move its own marker"
  drugs all land near 0; the engaged ones near 1. (HN, RO, Tm fall back to
  α = 1: no marker on the array, or no marker scored by any network.)
- **Why balanced accuracy cannot move**: the dial is pure shrinkage toward
  NORMAL. It can only delete change calls, never create a correct one, so it
  raises NORMAL recall (already the large class) and cannot touch the DOWN/UP
  recalls that dominate balanced accuracy. The arm's raw accuracy (0.9047)
  still sits below always-no-change (0.9385) for the same reason: on this
  dataset the model's bottleneck is **true-change recall under engaged
  drugs**, which no conditioning-by-shrinkage can supply.
- Fallback usage: per-(drug set, pathway) α available for 1,471 cases; the
  pooled per-drug-set α carried 7,421; 584 fell through to α = 1.

## What this redirects

1. **For the research axis**: stage 2 as imagined (a fitted partial-inhibition
   pin with re-solves) is also shrinkage-shaped at the perturbation source;
   it could differ through loops, but the burden of proof is now on
   mechanisms that ADD correct change calls, not ones that mute false ones.
   The record's standing candidates for that: the synthesis-stock class
   (specs/135: 527 wrong total calls on proteins the pathway only makes) and
   network reach for the drug-dose axes (specs/021: 116 of 190 decryptM sites
   unreachable in-pathway).
2. **For the product axis** ("analyzing datasets"): a 26:1 fixed-to-broken
   ratio and −5pp false alarms from one parameter the dataset itself supplies
   is a real usability gain even though it is not an accuracy-gate win — a
   user who uploads measurements should see the model go quiet where their
   own data says the perturbation never engaged. Serving the dial as an
   explicit, labelled conditioning mode (never silently) is a product
   decision, Adam's to make.

## Provenance

Control: build `20261026-1349_8b6527e`, solver `a490cdc`, arm
`k134_abund_a3/korkut_cases.tsv` (11,626 rows; 10,469 scored+measured; raw
0.812 verified before registration). Script: `bench/analysis/144/engagement_dial.py`
(markers imported from `bench/analysis/138/engagement.py`, never re-listed;
measured folds via the bench's own `conditions()` readers; seed 144).
korkut2015 is FIT (specs/051): development, not a generalization claim. The
CPPA training check was not run: the gate it would check already failed.
