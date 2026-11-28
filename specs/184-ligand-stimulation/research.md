# specs/184 research: ligand stimulation on the HPN-DREAM training data

Run 2026-11-27. Build `20261126-2322_d6d026b_top98` (Release 98), predictions `results/2e399999/stim184` (54 solves,
51 converged; 356 prediction rows in 9 networks). Scoring output is outside the repo, under `analysis/184/`.
**HPN-DREAM is FIT**: nothing below is a held-out test.

## Verdict

- **G1 fails, narrowly.** Model − reachability is +0.015 [−0.005, +0.036] over the 24 line–ligand conditions
  (by ligand: [−0.008, +0.041]). The model beats never-call by +0.079 [+0.037, +0.116].
- **The route carries more here than on Lun.** Enrichment 2.0 [1.2, 3.6] (Lun 1.5): 23.1% of reachable pairs change
  against 11.5% of unreachable ones. H1 holds.
- **Direction is right.** 76 of 82 measured changes the model calls are in the right direction (93%). H2 holds.
- **Precision is still the reached base rate.** 24.7% of the model's change calls are right; 23.1% of reachable
  pairs change. The specs/167 pattern holds on the fairest perturbation scored so far.

Registered reading: the route carries the signal and the solver adds little. But the traces below show two things the
aggregate hides: most misses are specific, fixable gaps in the inputs and the networks, and most false alarms are
cell-line differences that no context-free model can call.

## The registered figures

Six single ligands × four lines, early window (5–60 min). 700 scored pairs; 374 measured pairs have no prediction
(11 of the 47 readouts have no phospho form in Reactome). Measured: 120 UP, 11 DOWN, 569 no change.

| | Balanced accuracy | Change calls | Right | Precision | Recall | Direction |
|---|---|---|---|---|---|---|
| model | 0.412 | 308 | 76 | 24.7% | 58.0% | 76 of 82 |
| reachability | 0.397 | 438 | 94 | 21.5% | 71.8% | 94 of 101 |
| never | 0.333 | 0 | | | 0 | |

| Ligand | Pairs | Changes | Reachable | Enrichment | Model BA | Reach BA | Model precision |
|---|---|---|---|---|---|---|---|
| EGF | 122 | 24 | 95 | 3.1 | 0.458 | 0.403 | 19 of 67 |
| NRG1 | 114 | 31 | 95 | 1.9 | 0.423 | 0.362 | 21 of 64 |
| HGF | 114 | 34 | 95 | 1.5 | 0.376 | 0.366 | 22 of 68 |
| Insulin | 122 | 11 | 58 | 1.3 | 0.377 | 0.344 | 5 of 46 |
| FGF1 | 114 | 18 | 95 | 1.0 | 0.333 | 0.326 | 9 of 63 |
| IGF1 | 114 | 13 | 0 | – | 0.333 | 0.333 | no call |

| Line | Changes of pairs | Enrichment | Model precision |
|---|---|---|---|
| BT20 | 55 of 184 | 2.2 | 32 of 82 |
| MCF7 | 46 of 172 | 1.1 | 24 of 80 |
| BT549 | 18 of 172 | 5.0 | 13 of 73 |
| UACC812 | 12 of 172 | 6.8 | 7 of 73 |

**By engagement (report-only; receptor readouts left out).** Where the ligand raised its own receptor readout
(7 conditions): precision 42 of 105 (40%). Where it did not (5 conditions): 13 of 76 (17%). No receptor readout on the
array (FGF1, IGF1, insulin): 14 of 109 (13%). The same split as specs/179, in the same direction.

**Other registered variants.**
- Late window (2–4 h): model − reachability −0.017 [−0.053, +0.020]; enrichment 2.0.
- δ doubled: +0.041 [−0.006, +0.080]; enrichment 2.6.
- Serum (4 conditions): +0.003; enrichment 0.9.
- Withdrawal contrast: −0.065 [−0.120, −0.007], direction 48%. Pinning a ligand to 0 flips the PI3K switch in Signal
  Transduction (AKT pS473 reads 57× with HGF removed, specs/136's switch), so this contrast is not usable.

## Traced: where the 131 measured changes go

| | Pairs |
|---|---|
| called in the right direction | 76 |
| no route | 30 |
| silent with a route | 19 |
| wrong direction | 6 |

**No route (30).**
- **IGF1, 13.** IGF1's roots reach 11 reactions and stop at IRS binding and phosphorylation on the IGF1 receptor.
  Nothing connects those to PI3K or RAS in the Signal Transduction network, so IGF1 reaches none of the readouts:
  the continuation is curated under the insulin receptor's entities.
- **p90RSK pT359/S363, 9.** A phospho form exists in four networks and no ligand reaches it in any.
- The rest are spread over single readouts and ligands (insulin 5).

**Silent with a route (19).**
- **p70S6K pT389, 7 (and 4 of the no-route pairs).** An input error in specs/112's site map, not the model. The
  antibody's residue is T389 in isoform 2 numbering and T412 in the canonical sequence. Reactome curates the mTORC1
  target as `P23443-1` phosphorylated at T412 (R-HSA-165714); the map reads only the one `P23443-2` form at T389,
  which no ligand moves. bench/analysis/127's Lun table already carries T412. The same error sat in specs/112's
  scoring of the mTOR inhibitor, where p70S6K is the canonical mTOR target.
- S6 pS235/236 (2 silent, 2 no route) is one node that stays at 1 under every ligand; it sits behind the same kinase.

**Wrong direction (6).** HGF lowers MEK1 pS217/S221 in the model (0.85×, at the cutoff) and raises it in three
lines; NRG1 in UACC812 lowers EGFR, HER2 and SRC phosphorylation where the model raises them.

## Traced: the 232 wrong change calls

- **211 are UP calls on readouts that did not move.** The model makes one prediction per ligand–readout pair for all
  four lines. Of its 82 pairs, none is right in all four lines, 8 are right in three and 36 are right in none
  (of the 73 measured in all four: 29 in none, 20 in one, 14 in two, 8 in three). The lines differ: 30% of pairs
  change in BT20 and 7% in UACC812. A context-free model cannot have high precision here whatever its wiring.
- **15 are DOWN calls on unchanged readouts**, 12 of them c-Met pY1235 under EGF, FGF1 and NRG1 (0.06× in the model:
  other receptors' ligands drain what MET signalling shares). The EGF solve in Signal Transduction did not converge.
- p27 pT157 (20 calls) and STAT3 pY705 (16 calls) are called UP in every reachable condition and change in none.

## What this means for the next step

1. The perturbation's activity matters: enrichment rises from 1.5 (Lun) to 2.0 here and 3.1 for EGF, and precision
   doubles where the receptor demonstrably responded.
2. 31 of the 55 missed changes sit in three specific places (the p70S6K site map, IGF1's handoff, p90RSK; 35 with
   S6, which sits behind p70S6K), not in
   the propagation rules. They are checkable against Reactome one at a time.
3. The false alarms are mostly cell context. The measurement that separates them is on the array itself: the
   receptor's own phosphorylation.

Step 3 (a fresh ligand dataset scored once) should wait until the input errors in point 2 are corrected and
re-measured here as a post-hoc correction; otherwise the fresh test would be spent on a known defect.
