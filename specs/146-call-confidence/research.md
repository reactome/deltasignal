# specs/146 research: the confidence table, measured

## Verdict

**H1 fails as registered** — and the failure is one stratum, which redefines
what should be served.

| validation split | n | accuracy | mean claimed | ECE | AUC [95%] | gates |
|---|---|---|---|---|---|---|
| curator held-out | 18,481 | 0.9079 | 0.9213 | **0.0499** | **0.6075 [0.5280, 0.6786]** | both pass |
| experimental | 845 | 0.7609 | 0.8798 | **0.1491** | **0.3900 [0.2783, 0.5466]** | both fail |

Fit: 5,030 curator tuning cases; every case joined a basis (0 fallbacks).

## The trace: one stratum carries the whole failure

Per-stratum claimed vs experimental accuracy — every stratum roughly
calibrated except one:

- **`not / untouched` (fold exactly 1.0, called "no change"): claimed 0.961,
  experimental accuracy 0.297 (n = 138).** Against curators, an untouched
  readout matching "no change" is right 96% of the time (the wide curator set
  is mostly true negatives). Against experiment it is mostly WRONG — the
  specs/080 class: in the experimental benchmark an untouched readout usually
  means the network lost the route, not that nothing happens.
- Everything else: `down/railed` claimed 0.913 vs 0.959 observed; `up/gt3`
  0.913 vs 0.859; `down/le3` 0.907 vs 0.850; the low-claim strata observe low.

Excluding the untouched stratum, the same table is calibrated on **both**
axes: experimental ECE 0.0487, AUC 0.671 [0.482, 0.781] (n = 707; ten
pathways cannot tighten the interval); held-out ECE 0.0769, AUC 0.676
[0.613, 0.738].

## Amendment 1 (2026-11-07, recorded): what gets served

The registered serving semantics would have presented "no change, 96%
confident" for a readout the perturbation never reached — exactly the
misleading-by-precision pattern the honesty programme exists to kill. Amended
serving:

1. **An untouched readout (fold exactly 1.0) is served with `confidence:
   null` and `confidence_stratum: "untouched"`** — the client presents "the
   perturbation does not reach this readout in this network", not a
   confident no-change. The table still carries both measured rates for the
   stratum (curator 0.96, experimental 0.30) so the scope of the refusal is
   itself documented.
2. Every other valued call is served with the curator-record confidence and
   its resolved stratum, presented as what it is — agreement with the
   validation record — with both-axes per-stratum observed rates shipped in
   the table for anyone who looks.

This amendment changes presentation, not measurement; the gate verdict above
stands as the registered FAIL. No fresh axis exists to re-validate the
amended semantics against, and none is claimed: the served numbers are the
measured per-stratum rates themselves.

## Provenance

Fit and validation: frozen record `20261026-1349_8b6527e / 05dc867`
(curator_cases.tsv, experimental_cases.tsv); basis per (pathway, key_output)
from `bench/analysis/141/bench_readouts.json`; `railed` imported from
`score_124.py`; TUNING_PATHWAYS imported from `holdout_report.py`. Script:
`bench/analysis/146/fit_confidence.py`, seed 146. Registered at 3bb468c
before measurement.
