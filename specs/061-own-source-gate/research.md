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
