# specs/083 research: arm sa

Solver `eb7cdd0`, `--server DS_SELF_ACTIVATOR_WEIGHT=0.1`, on the canonical rebuilt for
specs/081, `20261020-2038_23ef3e1`, against `results/8190d83/dflt`, its current-defaults
run. `dflt` reproduces specs/081's `rcp` on every prediction: 24,100 curator and 849
experimental cases. The logs report 49,040 curator and 5,984 experimental self-activator
slot-solves.

## Results

| split | net (fixed / broken) |
|---|---|
| curator held-out | **0** (0 / 0) |
| curator tuning | 0 (0 / 0) |
| experimental (RAF set aside) | **0** (0 / 0) |

**No case changes class on either axis.** Predicted values do move:
- **curator:** 464 change, 322 of them by more than 5%. They are in pluripotent stem cell
  transcription 96, RUNX1 65, TP53 42, RUNX2 35, HOX 32 and NOTCH1 12. For example, CBFB
  knockout at a RUNX2 readout goes from 31.5× to 11.2×.
- **experimental:** 60 change, 40 of them by more than 5%, mostly TP53 knockouts from 0.000×
  to 0.001×.

The benchmark perturbs at 0 or 80× and scores three classes, so counting a gene once
or twice lands on the same side of the 0.85/1.15 cutoffs.

**The gate passes trivially** (held-out 0, experimental 0, no pathway moves).

## Predictions

1. PTEN overexpression in PIP3, at least 5 of 9 read DOWN: **fails**, as the disclosed
   check expected. The two-copy PTEN mRNA and the MKRN1 loop still rail it.
2. PTEN knockout stays UP: **holds**.
3. No held-out prediction: held-out is 0.

## Reading

The rule fixes a counting error in 251 reactions, and the magnitudes it changes are real,
but nothing the 3-class benchmark scores moves. Like specs/048's homodimer fix, it is a
correctness change with no measured accuracy effect. It matters for magnitude claims
(specs/021) and for any later fit to graded data. Adopting it is Adam's call, on
faithfulness.

## Adoption (2026-10-20)

Adam: "ok". `DS_SELF_ACTIVATOR_WEIGHT=0.1` is the default, as measured in `sa`; `off`
restores the old product. Because `sa` changed no class, the scored predictions of main
under the new default equal `results/8190d83/dflt`.
