# specs/176 research (2026-11-26)

Registration `2f8027fd`. Code: `bench/analysis/176/` (`budget_176.py` T1/T2/T4/T5/T6, `routes_176.py` T3,
`followup_176.py`, `axis_176.py`, `ledger_176.py`; tests `test_176.py`). Predictions: the specs/175 file (abundance
rule). Outputs outside the repo (`~/deltasignal-catalogs/analysis/176/`). Lun half B was not read.

## The answer in one paragraph

Lun 2019 mostly does not measure what the model predicts. The model asks "is this readout downstream of this protein
by a curated route?". The measurements are dominated by two other things: **whether the over-expressed protein is
active at all** in these cells, and **which readouts respond to almost anything**. Knowing those two (from the
protein's other readouts and the readout's behaviour under other proteins) predicts a change with AUROC 0.83. A curated
route scores 0.56, and adds nothing once the two are known. The largest single pattern is chemical, not a pathway:
**an over-expressed tyrosine kinase raises tyrosine phosphosites across the panel**, with or without a route.

## Registered tests

**T1. The measured call is reliable for strong effects and noisy at the margin.**
- a. In the same construct's EGF well, 71% of basal change calls repeat in the same direction. EGF itself moves
  p-ERK, p-MEK and p-p90RSK in control wells; leaving those out (post hoc), 76% repeat: 60% of marginal calls
  (1–1.5 δ), 77% at 1.5–2 δ, 87% at 2–3 δ, 99% above 3 δ. 20% of unchanged pairs are called changed in the EGF well.
- b. Strong effects sit on routes more often: 23% of marginal changes have a route, 41% of those above 3 δ (21% of
  unchanged pairs). With the call at 2 δ or 3 δ, route enrichment is 2.2 and 2.3 (1.55 at 1 δ). The model stays level
  with reachability (balanced accuracy 0.400 vs 0.399 at 2 δ).
- c. Narrower bins (1–3 vs 6–8) keep 64% of change calls: a smaller dose contrast, as expected.

**T2. Protein and readout explain the data; the route does not.** AUROC for "changed" on the 3,127 scored pairs:

| predictor | AUROC |
|---|---|
| protein's change rate on its other readouts × readout's change rate under other proteins (leave-one-out) | **0.82** |
| protein's rate alone | 0.77 |
| readout's rate alone | 0.66 |
| a curated route exists | 0.56 |
| the model calls a change | 0.55 |

Logistic regression, protein-grouped cross-validation (post hoc): protein + readout 0.830; adding the route 0.831;
adding the model's call 0.831. At 2 δ and 3 δ: 0.882 → 0.884 and 0.886 → 0.885. The protein's rate uses its other
measured readouts, which the model does not have; it stands for "is this protein active here".

**T3. The missing routes are not our construction.** Of 505 measured changes with no route, 501 have no route from
any form of the protein, pinned or not, in any of the 28 networks; 4 do. (Chains that cross from one top-level network
into another were not tested.)

**T4. p-AKT falls broadly, and the S6K feedback does not explain it.** p-AKT S473 falls under 80 of 288 proteins and
rises under 31. Proteins that lower it raise p-S6 / p-p70S6K / p-mTOR somewhat more often (odds ratio 1.75, Fisher
p = 0.058). They are active proteins in general: a median of 7 other readouts changed, against 3.

**T5. Withholding the known false-call sources costs more than it gains.**

| calls withheld when the route… | fixed | broken | recall | balanced accuracy |
|---|---|---|---|---|
| (none) | | | 0.194 | 0.384 |
| passes through PIP3 | 101 | 35 | 0.147 | 0.373 |
| is longer than 12 steps | 168 | 54 | 0.121 | 0.370 |
| either | 183 | 61 | 0.112 | 0.368 |

Each removes about three wrong calls per right one, but recall is already low, so none would pass a sensitivity
guard. Neither candidate from specs/174 is carried forward.

**T6. The model is right when the protein is active, and cannot tell when it is.**

| other readouts the protein changes | proteins | model change calls | right |
|---|---|---|---|
| 0 | 23 | 32 | 0% |
| 1–3 | 95 | 116 | 10% |
| 4–9 | 69 | 192 | 21% |
| 10 or more | 41 | 171 | 53% |

Of 288 measured proteins, 17 change nothing and 53 change ten or more readouts.

## Post hoc: the dominant pattern

- **One direction holds 54% of the variance** of the response matrix (288 proteins × 30 readouts, in units of δ). Its
  loadings are p-STAT3, p-STAT5, p-ERK1/2, p-FAK, p-SRC, p-PLCγ2 (positive) and p-AKT (negative). The next holds 12%.
- **Those are the tyrosine sites.** Residues from the canonical sequences; protein class from the networks (a curated
  catalyst of a reaction that adds a phosphotyrosine):

| | pairs | UP | DOWN |
|---|---|---|---|
| tyrosine kinase × tyrosine-only site | 240 | **56%** | 1% |
| tyrosine kinase × tyrosine + S/T site (ERK, p38, JNK) | 120 | **58%** | 3% |
| tyrosine kinase × S/T site | 800 | 14% | 10% |
| other protein × tyrosine-only site | 1,488 | 13% | 6% |
| other protein × S/T site | 4,960 | 10% | 5% |

- The class from the networks is crude: it takes every member of a catalyst complex (40 proteins). Known tyrosine
  kinases raise 6–9 of their 9 tyrosine readouts (FLT3, ALK, FGFR1 9 of 9; SRC, LCK, BTK, RET, NTRK1 and eight other
  kinases 8 of 9). Most members it sweeps in do not (BRAF, CDK2, MAP2K2, MAP2K7 0 of 9; the adaptor CRKL is the
  exception, 8 of 9), and nor do the kinase-impaired receptors ERBB3 (0) and EPHB6 (3). 24 of the 40 raise at least
  half.
- **Sizing on half A** (crude class): the model's calls plus "tyrosine kinase × tyrosine-containing site → UP" add
  154 calls, 88 right. Recall 0.194 → 0.312, precision 0.28 → 0.35, balanced accuracy 0.384 → **0.430**: +0.046
  [+0.028, +0.065] over the model and +0.050 [+0.026, +0.076] over reachability (protein bootstrap). It is the first
  thing on this dataset to beat reachability. It was found on half A, so it proves nothing until half B.
- **Caveat on what it is.** Promiscuous phosphorylation by a concentrated, self-activated kinase, and phosphotyrosine
  antibodies binding the over-expressed kinase's own homologous sites (p-SRC Y419 under any SRC-family kinase, p-BTK
  under TEC), cannot be told apart in these data. Either way it is what an over-expression experiment shows.
- A curated "catalyses its own phosphorylation" class (41 proteins) is active more often (mean 36% of readouts changed
  against 16%) but separates poorly (AUROC 0.66 for "changes ten or more"); withholding calls outside it is 166 fixed /
  46 broken with recall 0.132. Not carried forward.

## Ledger (each case once, first cause that applies)

**744 measured changes: 144 caught, 600 missed.**

| missed because | cases |
|---|---|
| marginal effect (under 1.5 δ; about 40% of these do not repeat) | 235 |
| tyrosine kinase raising a tyrosine site | 83 |
| no route; the protein changes ten or more other readouts | 152 |
| no route; other | 85 |
| route exists, model reads no change | 26 |
| route exists, wrong direction | 19 |

**331 change calls where nothing moved.**

| because | cases |
|---|---|
| the protein changes at most 3 other readouts (not active here) | 131 |
| route through PIP3 (the switch) | 67 |
| route longer than 12 steps | 48 |
| other | 85 |

## What is fixable

1. **Tyrosine-kinase promiscuity: yes, as a registered class rule.** Needs a clean definition fixed in advance (the
   protein is the active unit of a curated tyrosine-kinase activity, from Reactome's catalyst annotations), then one
   scoring on half B.
2. **Marginal effects: by reporting.** A secondary truth at 2 δ beside the registered call, as korkut is reported by
   engagement beside its raw score.
3. **Whether the protein is active: not from topology.** It is the largest factor and the model has no input for it.
   With a measurement of the protein's own activity it could be supplied (the engagement idea, specs/138 and 161).
4. **The switch and long routes: not by withholding calls** (T5).
5. **The other no-route changes of active proteins (152): not from Reactome.** No form reaches them in any network.

## Correction to specs/174

`READOUT_ACC` in `bench/analysis/174/trace_174.py` had no entry for STAT3, so p-STAT3 routes were not traced there
(specs/174 reported CREB and p70S6K as unmatched, not STAT3). Added here. With it, 575 of 732 pairs with a route are
traced. The specs/174 counts and conclusions do not change class; its route-length table was computed without STAT3.
