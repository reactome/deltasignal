# specs/124: how much do the results depend on the 100× ceiling?

**Status:** registered 2026-11-02, before any code and before any arm. Adam, 2026-11-02: the 100 "has been mentioned to
be arbitrary. I feel like if it is possible to address this concern we should", and, on the plan to measure first,
"feel free to make it happen". A sensitivity analysis, not a rule: nothing here is adopted.

**Why.** The solver works on an internal 0–1 scale and every node rests at 0.01 (the design spec's x₀, 1 on the UI
scale). The top of the scale is 1, so no node can rise more than 1 / 0.01 = 100-fold. Read biologically, the model
assumes every protein rests at 1% of the most it could reach. That number was never measured or fitted. It matters
directly where a value reaches the top (a positive loop that rails stops at the ceiling, specs/054, 082, 108), and for
magnitudes; class calls use the 0.85 / 1.15 cut-offs, far from it. Before the ceiling becomes a learned parameter, we
need to know how much of what has been concluded depends on it.

## Change (code)

- **`DS_BASELINE`** (server, default `0.01`, must lie in (0, 0.5]): the internal resting level the TSV parser gives
  every node (today a literal 0.01 in two places), and the fallback wherever the solver reads a baseline (five
  `get(..., 0.01)`). The ceiling is then 1 / `DS_BASELINE` fold. The UI scale stays internal × 100, so at the default
  nothing changes; the solve and the parse report `resting_baseline`.
- Nothing else is rescaled. The de-repression limit (`DS_HILL_SAT_H_MAX`, 10×) and the depletion floor (0.1×) are
  fold limits of their own and stay as they are; they are separately arbitrary and not measured here.
- **Bench** (`benchmark_vs_mpbiopath.py`, both axes): reads `DS_BASELINE` (default 0.01) and refuses unless every solve
  reports the same `resting_baseline`. A readout's fold is its internal value / `DS_BASELINE` (today × 100, identical
  at the default). The over-expression pin is fold min(80, 0.8 · C) for ceiling C: 80× wherever the scale allows it,
  and 0.8 of the top otherwise (today's pin is 80× = 0.8 of the top), so C = 10 pins 8× and C = 30 pins 24×. Knockouts
  pin 0. Cut-offs stay 0.85 / 1.15.

## Arms

On the frozen canonical build `20261026-1349_8b6527e`, every other setting the default, both axes (curator and
experimental), each with `run_arm.sh --server DS_BASELINE=b --bench DS_BASELINE=b`:

| arm | `DS_BASELINE` | ceiling C | over-expression pin |
|---|---|---|---|
| `c100` (control, new code) | 0.01 | 100× | 80× |
| `c10` | 0.1 | 10× | 8× |
| `c30` | 1/30 (0.0333333333) | 30× | 24× |
| `c1000` | 0.001 | 1000× | 80× |

**D0 (must hold before anything is read):** `c100` reproduces the frozen canonical predictions case for case (status,
value to print precision, call). If not, the code change is wrong and the arms are void.

## Measures (all reported; no gate, no adoption)

For each arm against `c100`:
- accuracy on curator held-out, curator all and experimental, with RAF set aside and with (specs/054 D1); cases that
  change class, net, McNemar p; per-pathway net and how many pathways move;
- **rails**: share of scored readouts at the top (fold ≥ 0.99 C) or at the bottom (fold ≤ 0.01);
- how many changed cases are rail-dependent by specs/108 (`rail_audit.py` on the canonical build);
- magnitudes: Spearman correlation of log fold over the readouts scored in both arms, and the share of readouts whose
  fold moves by more than 10%;
- convergence counts.

The c10 arm changes two things at once (ceiling and over-expression dose, since 80× cannot exist under a 10× ceiling);
its over-expression cases are reported separately from its knockouts. c1000 changes only the ceiling.

## Predictions (written before any run)

- **P1.** c30 and c1000 each change the class of fewer than 2% of scored curator held-out cases and of scored
  experimental cases.
- **P2.** c10 changes more, between 2% and 8% on each axis, mostly over-expression cases.
- **P3.** At c30 and c1000, more than half of the changed cases are rail-dependent by specs/108.
- **P4.** The share of readouts at a rail moves by less than 5 points between c30 and c1000: a rail moves with the
  ceiling; it does not go away.
- **P5.** Magnitudes do depend on the ceiling: under c1000, more than 10% of scored readouts move by more than 10% in
  fold against c100.

## Reading, fixed in advance

- If P1 holds, the up/down results do not rest on the 100, and the ceiling can enter the learning plan as one more
  parameter (a resting level per node class, specs/125).
- If P1 fails, the headline numbers depend on an arbitrary constant; that is reported beside every headline before any
  further work, and the ceiling is resolved first.
- Whatever P5 shows applies to every correlation-type score (CPPA, decryptM magnitudes): those are reported with the
  ceiling named.
- No ceiling is chosen from these scores: the curator and experimental sets are SELECT data (specs/051).

## Amendment 1 (2026-11-02, after D0 and the first reading, before these arms run)

c30 and c10 change two things at once, the ceiling and the over-expression dose (24× and 8×). The registration flagged
this only for c10. A trace shows that at least one c30 change is the dose alone: MIB1 over-expression in NOTCH1 reads
2.075 at 80× and exactly 1.000 at 24× and at 8× under the default ceiling. Two dose-only arms at the default ceiling
separate the two effects:

| arm | `DS_BASELINE` | ceiling | over-expression pin |
|---|---|---|---|
| `d24` | 0.01 | 100× | 24× (`--bench DS_PERTURB_UI_UP=24`) |
| `d8` | 0.01 | 100× | 8× (`--bench DS_PERTURB_UI_UP=8`) |

Read: a c30 (c10) change that d24 (d8) also shows is a dose effect. A change in c30 (c10) but not in d24 (d8) is a
ceiling effect at that dose. The ceiling-only contrasts are c30 vs d24 and c10 vs d8. P1 and P2 are still read as
registered, against c100; the decomposition is reported beside them. c1000 needs no dose arm, because its pin stays
at 80×.
