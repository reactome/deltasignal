# 027 — How a readout with several node copies is read

## Why (the canonical baseline's propagator failures)

- **Where it shows.** 262 held-out cases are missed changes where both our
  network and Reactome have a correctly signed route (`case_oracle.py`). At
  least three quarters of them predict exactly 1.0: the signal vanishes rather
  than fading. 184 are knockdowns expected DOWN.
- **What the solves show.** In a random sample of 40, re-solved on the
  canonical build, **20 have a readout copy that does move in the expected
  direction**, and the benchmark's `DS_KO_AGG=max` hides it: PALB2 KD reads
  [0, 1, 1, 0] across the readout's copies, and MAP2K1 KD reads [0, 1]. The
  copies still at 1 are not downstream of the pin, for example variants built
  from a different set member.
- **Why `max` is the default.** It was chosen under the old broad pins, which
  set many copies directly. It has not been measured under root pinning.

## Pre-registration (committed before any arm runs)

- **Arms**, bench-side, on the canonical build, paired against
  `184dfd5/baseline`: `--bench DS_KO_AGG=mean`, `min` and `extreme` (the copy
  that deviates most from baseline).
- **Decides:** curator held-out in-release net, McNemar p, pathways moved,
  perturbations. The experimental axis is reported.
- **Adopt** only if held-out net > +15 with p < 0.05, the gain is not one
  pathway, and experimental is no worse than −15. If more than one mode
  passes, adopt the one with the larger held-out net.
- **The trade to watch.** `min` and `extreme` let one moving copy decide the
  readout. They can create false changes where a copy moves for an unrelated
  reason. False change is reported separately.

## Result (build `20260925-2326_d4f4f64`, arms at `e3685f5`, against `184dfd5/baseline`)

| `DS_KO_AGG` | held-out net | p | pathways up / down | perturbations | tuning | experimental |
|---|---|---|---|---|---|---|
| **mean** | **+65** (132 / 67) | 5e-6 | 17 / 6 | 132 | −8 (ns) | 0 (15 / 15) |
| min | −33 (254 / 287) | 0.17 | 3 / 5 | 332 | −6 | −9 |
| extreme | +29 (153 / 124) | 0.09 | 18 / 8 | 161 | −2 | +3 |

- **`mean` is adopted as the benchmark default.** It meets every
  pre-registered condition; the largest single pathway is Pre-NOTCH at +16.
- **The trade.** On the held-out, in-release cases it fixes 132 misses and
  creates 64 false changes. `extreme` fixes 153 misses but creates 124 false
  changes. `min` breaks far more than it fixes.
