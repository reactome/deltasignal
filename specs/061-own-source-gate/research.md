# specs/061 research

## R1. Arm `rtwin58g`: the gate fails; the MDM2 prediction failed

`scripts/run_arm.sh rtwin58g --server DS_RECYCLED_TWIN=on --server
DS_DEPLETION_OWN_PRODUCT=suppress_only_source --bench
DS_PIN_ROOT_IGNORE=depletion` at DS `cb53f3f`, results
`builds/20261012-1401_01206ca/results/cb53f3f/rtwin58g`. Scored with
`holdout_report.py --pair-on-valid --set-aside RAF_MAP_kinase_cascade`
(RAF moves nothing here; both readings identical).

| against | held-out | tuning | experimental |
|---|---|---|---|
| control `64d23cb` | **+24** (47 / 23, p 0.0056) | −26 (21 / 47) | **−9** (3 / 12, p 0.035) |
| `rtwin58s` | +3 (26 / 23, p 0.78) | −138 | −16 (2 / 18, p 0.0004) |
| `rtwin58e` | −7 (31 / 38) | +7 | +8 (8 / 0, p 0.0078) |

Held-out in release 89.46% (control 89.33%). Convergence: curator 1,655 of
1,725 (control 1,652), experimental 218 of 244 (control 215).

**Gate: fails** (experimental loss significant, p 0.035; TP53 −42 and DSB
−14 lose more than 10).

Against the control, by perturbation (curator, net):
- TGF-beta SMAD2/SMAD3: **0** — the registered fix held (`rtwin58s` −12);
- MET OE +10, CBL KO/OE +11, EPS15 OE +5, SRC +4, PTCH +8 (the specs/058
  pins, as in `rtwin58`);
- AKT1/AKT2 KO +5 each (`rtwin58s` +38 each), MDM4 +3;
- **MDM2 KO (TP53) −42**, experimental −9: the prediction failed;
- RNF168 OE (DSB) −16: the specs/057 non-converged coin; SH3GL1 −5.

### The MDM2 knockout, traced (`bench/analysis/061/mdm2_trace.jl`)

The TP53 solve under the MDM2 knockout does not converge under `full`,
`suppress_only` or `suppress_only_source` (500 iterations each). Under
`suppress_only_source`, TP53 copy 8c4aa191 (R-HSA-6793669) reads 8e-8 fold;
it is depleted by three entity-own depleters (the MDM2:TP53 pool 80e7a4cf and
two others), and in all three **its only source copy is itself**. Once the
iterate dips below baseline, the rule forbids de-repression, which keeps it
low, which keeps the cap: a self-confirming low branch. The readouts flip
from UP (canonical 2–58×) to DOWN (0.18×), not to the 1.0 the static entity
test gave.

The rule as registered is circular whenever the target is its own source
(202 of 315 pairs; 191 have no other source): the target's own value
includes the depletion being decided, so it cannot say which partner fell.
The SMAD2 recycled copies (7632c373 and two others) have a separate source,
the root copy c854b46a, which is why TGF-beta was fixed.

`ownsrc` (the rule alone) is reported in R2.

## R2. Arm `ownsrc` (the R1 rule alone): fails, same defect

`scripts/run_arm.sh ownsrc --server DS_DEPLETION_OWN_PRODUCT=suppress_only_source`
at `cb53f3f` (the pre-amendment rule), results `results/cb53f3f/ownsrc`.

| against | held-out | tuning | experimental |
|---|---|---|---|
| control `64d23cb` | −9 (14 / 23, p 0.19) | −34 (13 / 47) | **−11** (2 / 13, p 0.0074) |
| `ownsup` (`cc9ab80`) | −9 | −80 | −15 (2 / 17, p 0.0007) |

Convergence: curator 1,653, experimental 219 (control 1,652 / 215).
**Gate: fails.** The losses against `ownsup` are the self-sourced cap
(R1's MDM2 mechanism) undoing `suppress_only`'s TP53 gains. Not re-run after
amendment 1; the amended rule is measured only in the registered
combination (`rtwin58g2`).

## R3. Amendment 1, arm `rtwin58g2`: the registered gate passes

`scripts/run_arm.sh rtwin58g2 --server DS_RECYCLED_TWIN=on --server
DS_DEPLETION_OWN_PRODUCT=suppress_only_source --bench
DS_PIN_ROOT_IGNORE=depletion` at DS `9a0b839` (rule code `5992c91`), results
`builds/20261012-1401_01206ca/results/9a0b839/rtwin58g2`. RAF moves nothing;
both readings identical.

| against | held-out | tuning | experimental |
|---|---|---|---|
| control `64d23cb` | **+32** (33 / 1, p < 1e-4) | **+112** (113 / 1) | **+7** (7 / 0, p 0.016) |
| `rtwin58s` | +11 (12 / 1, p 0.0034; TGF-beta only) | 0 | 0 |
| `rtwin58g` | +8 (22 / 14, p 0.24) | +138 | +16 (18 / 2, p 0.0004) |

Headline, in release: curator held-out **89.51%** (control 89.33%), tuning
80.92% (78.70%), all 87.67%; experimental 7 more correct (598 / 849, from
591). Convergence: curator 1,645 of 1,725 (control 1,652), experimental 214
of 244 (215); the bench does not record which solves, and no case moved
outside the perturbations below.

**Gate (registered, unchanged): passes.** Held-out net > 0 at p < 1e-4;
no experimental loss (7 / 0); worst pathway −1 (TGF-beta, SMAD7 KO, one
case; CHEK2 KO in TP53, one case); no experimental loss to trace.

Against the control, by perturbation (curator net): AKT1 KO +38, AKT2 KO
+38, MDM4 KO +32, BRCA1 KO +3 (TP53, tuning); MET OE +10, CBL KO/OE +11,
PTCH KO/OE +8, BRCA1/BARD1 +4 (DSB), MDM2 OE +2 (Cell Cycle Checkpoints);
SMAD2/SMAD3 0. Experimental: AKT1 +2, MDM4 +2, AKT2 +1, BRCA1 +1, MDM2 OE
+1.

Predictions: MDM2 KO returned to canonical (as predicted), TGF-beta held
at 0 (as predicted), AKT1/AKT2 recovered to +38 each (as predicted).
**Wrong:** I predicted DSB would fail the gate on RNF168; RNF168 did not
move in this arm (the coin landed as in canonical). That is luck, not a fix:
the RNF168 solve is still non-converged (specs/057 R1).

### What this means

- `rtwin58g2` is `rtwin58s` (specs/059 R2) plus the 12 SMAD2/3 cases, with
  nothing else moving: the three rules together (specs/058 pins, specs/059
  twin edge, this amended own-product rule) do what each step was built for.
- It is post hoc (four amendments' worth of looking at these sets: 058 →
  059 → 060 → 061 → amendment 1). The held-out gain is concentrated: 33
  fixed cases in 5 pathways, 22 readouts, 6 genes. The spec's own condition
  applies: **I recommend an independent check before adopting it as the
  default** (Adam's call): the korkut2015 TEST set scored once, or a fresh
  Reactome release.
- Adoption would mean three changes: two solver defaults
  (`DS_RECYCLED_TWIN=on`, `DS_DEPLETION_OWN_PRODUCT=suppress_only_source`)
  and one benchmark protocol change (`DS_PIN_ROOT_IGNORE=depletion`: a form
  fed only by depletion counts as a root, as in MP-BioPath's networks), which
  is a protocol-of-record decision (Adam, 2026-09-25) and not only a rule.

## R4 (2026-10-15): amendment 2, arm `rtwin58g3` — identical to `rtwin58g2`

Code `639787b` (the two review fixes, amendment 2), canonical build, same
overrides. Results `639787b/rtwin58g3`.

| against the control (RAF set aside where it applies) | net | fixed / broke | p |
|---|---|---|---|
| curator held-out | +32 | 33 / 1 | < 1e-4 |
| curator tuning | +112 | 113 / 1 | < 1e-4 |
| experimental | +7 | 7 / 0 | 0.016 |

- Against `rtwin58g2`: **0** predicted classes differ on either axis (24,100
  curator, 849 experimental cases), validity identical; the largest
  `pred_ui` change is 1.9e-5. Convergence identical (curator 1,645 / 1,725,
  experimental 214 / 244).
- So both defects are real in the code (each has a test that fails on the
  old code) but latent on this catalog: no carrier free form without a root
  twin sat next to a recycled copy in a scored solve, and no modifier-only
  copy changed a gate that mattered.
- The registered gate passes as in R3. The candidate is the amended code.

## R5 (2026-10-15, after adoption): most tuning and experimental gains are iteration-cap artifacts

The week's adversarial review (finding: "the source gate puts TGF-beta into a
limit cycle") showed `_own_source_down` is an on/off switch read on the live
iterate, so where a source copy sits in its target's cyclic component the
cap toggles between sweeps and the solve never settles; the scored value is
wherever the 500-iteration cap stops. Census: 11 of 97 gated pairs, in 3
pathways (TGF-beta 4, R-HSA-8878166 6, R-HSA-5358351 1).

Check (`bench/analysis/061/capcheck.jl`, results `capcheck_9e93c11.tsv`):
every case whose class differs between the control (`64d23cb`) and the
adopted defaults (`9e93c11`), 176 in all, re-solved under the defaults at
iteration caps 500, 1000, 2000 and 4000. At 500 every case reproduces the
benchmark exactly.

| | fixed / broke | fixes converged and stable at every cap |
|---|---|---|
| curator held-out | 33 / 1 | **33 of 33**; the 1 break (SMAD7 KO, TGF-beta) is unconverged and flips with the cap |
| curator tuning | 113 / 1 | **2 of 113**: 111 are AKT1, AKT2, MDM4 and BRCA1 knockouts in Transcriptional Regulation by TP53 whose solves never converge; their class changes with the cap |
| experimental | 7 / 0 | **1 of 7**: 6 are the same unconverged TP53 solves |

So the held-out +32 that decided the gate is real (all converged, all stable),
and the reported tuning +112 and experimental +7 are not: they are the
iteration cap's choice in an oscillating TP53 solve. The experimental gate
condition (no significant loss) still holds; the experimental gain does not.
The adoption decision rested on "passes more tests"; on stable cases it still
does (held-out +33 stable fixes against 0 stable breaks), but it introduced
non-convergence in TP53 and TGF-beta, which amendment 3 addresses.
