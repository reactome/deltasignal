# specs/177 research (2026-11-26)

Registration `a0968539` (rule, class tables, scorer and tests pushed before half B was opened). The half-B file's
sha256 matched the sealed record (`6ae9507d…`). `score_177.py B` ran once; output outside the repo
(`~/deltasignal-catalogs/analysis/177/score_b.json`). `lun2019_b` is now registered FIT.

## Half B (230 proteins, 3,633 scored pairs; 463 UP, 229 DOWN, 2,941 no change)

Positive control passes: 10 of 13 measured self pairs read UP.

**Gates**

| gate | result | |
|---|---|---|
| **primary:** candidate − reachability (balanced accuracy) | **+0.066 [+0.036, +0.096]** | **passes** |
| secondary: candidate − model | +0.063 [+0.039, +0.085] | passes |
| secondary: added calls right vs wrong | 89 vs 14 (103 calls, 20 kinases), p = 8e-15 | passes |

**Arms**

| | balanced accuracy | change calls | recall | precision |
|---|---|---|---|---|
| candidate (model + rule) | **0.420** | 599 | 0.249 | 0.287 |
| model (frozen, abundance) | 0.358 | 496 | 0.120 | 0.167 |
| reachability | 0.354 | 749 | 0.168 | 0.155 |
| the rule alone | 0.426 | 159 | 0.192 | **0.836** |
| always no change | 0.333 | | | |

With the measured call at 2 δ: candidate 0.487, reachability 0.360 (+0.127 [+0.075, +0.181]). At 3 δ: 0.514 against
0.362 (+0.152 [+0.095, +0.210]).

**The model alone, as specs/173 registered for half B**
- model − reachability +0.004 [−0.013, +0.020]: unresolved;
- model − always no change +0.025 [+0.003, +0.046]: a gain;
- route enrichment 1.50 [1.17, 1.88] (half A: 1.55 under the abundance rule): routes carry signal, short of the
  registered 1.5 bar at the interval's lower end, as on half A;
- lift +0.009 [−0.022, +0.041]: unresolved.

**Replication of specs/176 (no gate)**

| | half A | half B |
|---|---|---|
| tyrosine kinase × tyrosine site: share UP | 72% (306 pairs) | **72%** (234 pairs) |
| tyrosine kinase × S/T site: UP | 14% | 18% |
| other protein × tyrosine site: UP | 14% | 12% |
| other protein × S/T site: UP | 9% | 9% |
| kinases raising at least half of their tyrosine readouts | 27 of 34 | 19 of 26 |
| AUROC: protein rate × readout rate (leave-one-out) | 0.82 | 0.82 |
| AUROC: a route exists | 0.56 | 0.55 |
| AUROC: the model's call | 0.55 | 0.53 |
| AUROC: the candidate's call | 0.61 | 0.60 |
| basal change calls repeating in the EGF well | 71% | 68% |

## Reading

- **The rule holds on unseen proteins.** The primary gate and both secondary gates pass, and the class table is the
  same on both halves. An over-expressed tyrosine-specific kinase raises a tyrosine site 72% of the time, with or
  without a curated route; the rule alone is right on 84% of its calls.
- **The model itself replicates as level with reachability** (+0.004), a small gain over always-no-change, and a route
  enrichment of 1.5. specs/176's structure replicates too: protein activity and readout responsiveness explain the data
  (0.82), a route adds little (0.55).
- **What the rule is not.** It is not signal propagation, and these data cannot tell promiscuous phosphorylation from
  phosphotyrosine antibodies binding the over-expressed kinase's homologous sites. It describes what an
  over-expression experiment with these reagents shows.
- Even with the rule, balanced accuracy is 0.42: recall of real changes is 25%. The rest is what specs/176 listed:
  marginal effects, whether a protein is active at all, and changes with no curated route.

## For Adam's decision

1. Adopt the candidate (abundance + tyrosine-kinase rule) as the Lun bench protocol of record.
2. Serve it as a statement on over-expression queries: "a tyrosine kinase over-expressed at this level raises tyrosine
   sites beyond its curated routes", with the measured 72% / 84%.

Neither is done here.
