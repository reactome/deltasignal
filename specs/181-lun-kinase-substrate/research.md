# specs/181 research (2026-11-27)

Registration `d4467a10`. Predictions: Release 98 (`20261126-2322_d6d026b_top98`, `results/906698aa/lun178_r98`), both
halves pooled unless stated; both are FIT. Table: OmniPath enzyme–substrate, human, downloaded 2026-11-27 (41,505
rows, 39,037 phosphorylation and 1,072 dephosphorylation; sha256 `503f8028…`; kept outside the repo). Output outside
the repo under `analysis/181/`.

## Where the errors are (the starting ledger)

| measured changes (1,437) | count |
|---|---|
| no curated route | 1,001 |
| route, model right | 224 |
| route, model silent | 125 |
| route, wrong direction | 87 |

| model change calls (1,012) | count |
|---|---|
| right | 226 |
| readout measured no change | 699 |
| readout changed the other way | 87 |

## Part K: kinase–substrate edges do not cover the gap

371 of Lun's proteins are enzymes in the table.

| arm | pairs added | no-route changes reached (of 1,001) | added pairs changed | direction right | null 95th pct. | above null |
|---|---|---|---|---|---|---|
| K1, exact site | 136 | 39 (3.9%) | 29% | 24% | 23% | yes, barely |
| K1p, any site of the readout's protein | 295 | 77 (7.7%) | 26% | 20% | 20% | yes, barely |
| K2, one hop into the curated network | 1,603 | 420 (42%) | 26% | 19% | 19% | no |

Base rates among the 5,275 no-route pairs: 19% changed, 13% UP.

- **Registered reading: K1 does not cover the gap** (3.9% against the registered 20%). Per half it is 4.0% and 3.8%,
  and on neither half alone is it above the null.
- 160 of K1's 296 calls on scored pairs already have a curated route: Reactome holds more than half of the known
  direct relations for this panel. The table adds little that is new.
- K2 reaches 42% of the gap but at the relabelled-null rate: it reaches a lot because it calls a lot.
- Balanced accuracy: adding K1 to model + tyrosine-kinase rule changes it by −0.0005 [−0.004, +0.003]; K2 by −0.011
  [−0.025, +0.002].

**The missing changes are not known direct kinase → site relations.** A direct substrate changed 29% of the time,
against 19% for an arbitrary unrouted pair.

## Part M: the route-exists errors

**Silent (125).** From the prediction rows:

| | count |
|---|---|
| route exists, predicted fold exactly 1 | 64 |
| fold moves the right way, under the 15% cutoff | 31 |
| fold moves the wrong way, under the cutoff | 26 |
| networks disagree | 4 |

- 119 of the 125 are reached in one network only, mostly Signal Transduction; p-p70S6K (25), p-AMPKα (16), p-S6 (16).
- Four exact-1 cases traced by hand (ALK → p-S6, ACVR2B → p-S6, ACVR1B → p-AMPKα, CRKL → p-S6): the shortest routes are
  38 to 58 edges long and run through transcription of another ligand, a second receptor, PLCγ, then small-molecule
  second messengers (calcium) into the mTOR branch. These are graph connectivity, not a mechanism anyone would claim;
  the model reading 1 is the defensible answer and the measured change comes from somewhere else.

**Wrong direction (87).** 68 are called UP and measure DOWN; 36 of those are p-AKT, 9 p-ERK1/2, 8 p-GSK3β.
- p-AKT S473 is the readout that falls most often in the whole screen (106 pairs); p-S6, p-AMPKα, p-4EBP1 and
  p-p70S6K follow. Tyrosine kinases are 24% of the proteins behind measured falls.
- The model routes a receptor or kinase through PI3K to AKT and calls UP. The measured fall has no curated cause
  (131 of 435 measured falls have any route; 15 are called right).
- 77 of the 87 wrong-direction calls are strong or railed folds, as are 578 of the 699 false alarms: magnitude does
  not separate right from wrong.

**By route length** (`routes_181.py`; reaction steps on the shortest route; 298 routed pairs the tracer could not
match are left out):

| reaction steps | pairs | measured changed | calls right |
|---|---|---|---|
| 1–5 | 317 | 34% | 29% |
| 6–10 | 267 | 27% | 18% |
| 11–20 | 284 | 26% | 20% |
| 21+ | 298 | 26% | 16% |

Base rate over all pairs: 21%. **Even at one to five reaction steps, two readouts in three do not move, and 176 of
280 calls there are on readouts that did not change.** Long routes are barely above the base rate.

## Reading

1. The 70% of changes without a route are not missing kinase–substrate facts. With specs/176 (protein activity ×
   readout responsiveness, AUROC 0.82; a route adds nothing) the consistent picture is that an active over-expressed
   kinase disturbs the phospho-network broadly, and a small set of responsive readouts take most of it.
2. Where a route exists, the dominant error is not attenuation or direction: it is that the over-expressed protein
   does nothing down a short, correct route two times in three. Reactome records what an *active* protein does; the
   screen over-expresses proteins that are mostly not active.
3. The "silent" class is mostly spurious connectivity and is not a solver defect. The wrong-direction class is mostly
   one unexplained biological effect (p-AKT falls under kinase over-expression).

## What this leaves

- Not worth pursuing on this dataset: more edges of any kind (specs/113, 166–172 and K here), route-length or
  magnitude cutoffs (specs/176 T5 and the table above).
- Worth testing next, from structure alone: **whether the protein needs activating.** Reactome distinguishes a protein
  whose starting form is itself the catalyst from one that must first be modified or bound by something upstream. If
  short-route false alarms concentrate in the second kind, an over-expression query can say "this needs an upstream
  signal" instead of calling a change. Both halves are FIT, so it would be sized here and confirmed on new data.
- The p-AKT fall is a separate, unexplained effect and a candidate for its own trace.
