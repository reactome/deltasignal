# specs/115 research

Case logs: specs/114's top-level runs (`20261027-1346_8b6527e_top114`, code `a47898d`, `s114_top_kb`). Hop table:
`~/deltasignal-catalogs/analysis/115/korkut_hops.tsv`.

## 1. Where the errors fall

| class | korkut | decryptM (kinobeads) |
|---|---|---|
| unchanged, model says changed | **3,845 (78% of errors)** | **1,055 (66%)** |
| changed, no route from a target | 823 (17%) | 391 (25%) |
| changed, route, model says NOT | 149 (3%) | 66 (4%) |
| changed, route, wrong direction | 129 (3%) | 81 (5%) |
| changed, right | 444 | 187 |
| unchanged, right | 12,446 | 3,316 |

Wrong direction is 3-5% of errors; the specs/114 amendment 1 count (73 of the 118 cases only reachability gets right)
described the cases the model loses to reachability, not the model's errors. The largest class by far is a change
predicted for a readout that did not move. Those calls are not near the cutoff: on korkut, of 3,845, 1,185 predict a
2-10x change and 978 more than 10x (369 read 0).

**A change call carries almost no information about whether the readout changes.** korkut: 8.7% of cases change; of
the model's 4,418 change calls, 10.0% are right (reachability 9.5% of 6,353). decryptM: 14.2% change; 14.1% of the
model's 1,323 change calls are right (reachability 16.1%).

## 2. By distance from the target (korkut)

| hops | cases | changed | model change calls | right | direction agrees where both change |
|---|---|---|---|---|---|
| 1-5 | 1,975 | 14% | 1,532 | 13% | 200 / 263 (76%) |
| 6-10 | 968 | 20% | 767 | 21% | 162 / 187 (87%) |
| 11-20 | 1,670 | 8% | 1,132 | 7% | 77 / 101 (76%) |
| 21-40 | 1,261 | 5% | 738 | 0.7% | 5 / 21 |
| >40 | 479 | 13% | 188 | 0% | 0 / 1 |
| no route | 11,483 | 7% | 61 | 0% | — |

Within 20 hops the model gets the **direction** right in 76-87% of the cases where both it and the measurement call a
change, but it calls a change for 74% of those readouts (3,431 of 4,613) while 13% change: it cannot say **which**
readouts move. Beyond 20 hops it calls UP 712 times and DOWN 214, while measured changes there are mostly DOWN (108 of
128); its calls there are worse than chance (balanced accuracy 0.27 at 21-40 hops, 0.19 beyond). 61 change calls have
no route in the bench's adjacency (not traced).

## 3. Other methods (literature, read-only; 2026-10-27)

No published prior-knowledge-only method reports a held-out readout prediction score on korkut2015, decryptM or P100.
Every score found is from a model fitted to the dataset:
- korkut2015: belief-propagation ODE models with Pathway Commons priors, leave-one-drug-out correlation 0.81-0.84
  (Korkut et al. 2015, eLife); CellBox 0.78 and a linear regression 0.78 leave-one-drug-out (Long et al. 2025, BMC
  Bioinformatics; CellBox figures not read from Yuan et al. 2021 itself). Continuous correlations, not three-class calls.
- decryptM (Zecha et al. 2023, Science): no modelling benchmark (Herbst & White commentary).
- P100: ternary ILP on a 35-node network fitted to P100, about 84% leave-one-out (Ji et al. 2017, BMC Syst Biol).
- Prior-knowledge-only, nearest analogue: CausalPath (Babur et al. 2021, Patterns), 29 of 32 checkable relations in the
  predicted direction on breast RPPA; target-inhibition effects explained for 14 of 31 drugs on PC3 phosphoproteomics.
- Coverage: under 1% of measured phosphosites map to Reactome (Huckstep et al. 2021, PeerJ); under 5% of sites have a
  known kinase (Needham et al. 2019, Sci Signal; abstract only). Adaptive feedback: BRAF inhibition reactivates EGFR
  and ERK (Prahallad et al. 2012, Nature; abstract only).

## Reading

The deficiency is in what a generic network plus a full-strength perturbation can say, not in the sign rules:
1. **Over-propagation.** Every curated route is taken as fully effective in every cell, so most readouts near a target
   are predicted to move; in the cell few do. The model has no cell context (expression, activity, mutations such as
   BRAF V600E in SKMEL-133), pins the drug's target as a near-complete knockout, and has no adaptation (korkut reads at
   24 h).
2. **Long routes are wrong.** Past about 20 hops the calls are mostly UP and worse than chance: de-repression and loop
   artefacts accumulate along long paths through the large networks (the rail class of specs/108).
3. **Coverage.** 64% of korkut cases and 63% of decryptM cases have no route from a target; 17-25% of errors are
   changed readouts there.

Candidate changes (each would be a new spec, pre-registered, fitted only on FIT data and judged on unseen data): a
per-hop attenuation or a distance limit; a partial-inhibition dose instead of a knockout; restricting networks to
proteins expressed in the cell line (needs an external expression dataset, so Adam's approval).
