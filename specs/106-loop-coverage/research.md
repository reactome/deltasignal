# specs/106 results

Catalogs: benchmark `20261024-2211_7f4b28a` (92 pathways); sample `20261025-1318_7f4b28a_sample101`
(specs/097's 101 other pathways, same generator commit). Script `loop_coverage.jl`; per-loop and
per-pathway tables `bench92_*.tsv`, `sample101_*.tsv`. The older sample build
(`20261023-1648_c15e27d_sample101`, before specs/100–102) gave almost the same figures (358 loops,
286 dissolved, 55 untouched).

## Coverage

| | benchmark (92) | sample (101) |
|---|---|---|
| pathways with a loop | 71 | 53 |
| loops (cyclic components, network as generated) | 744 | 361 |
| touched by at least one rule | 690 (93%) | 306 (85%) |
| dissolved (nothing cyclic left) | 657 (88%) | 290 (80%) |
| partly handled | 33 | 16 |
| untouched | 54 | 55 |
| loop nodes in dissolved loops | 46% | 51% |
| loop nodes left in remnants | 37% | 42% |

By loop count the rules handle most loops in both catalogs. By node count about 40% of loop nodes
remain in unhandled remnants in both, concentrated in a few large components (below).

## Loop types

| type | benchmark: loops / pathways | sample: loops / pathways | handled by | why it is handled this way |
|---|---|---|---|---|
| carrier binding/release (E → E:X → E) | 564 / 52 | 194 / 30 | `DS_CARRIER_LOOP_MODE=closed_twin`, `DS_CARRIER_LOOP_OPEN=fed_single_ub`, `DS_CARRIER_COMPONENT_TWIN` (specs/049, 077, 084, 088) | a protein bound and released unchanged is conserved; its free form follows its own gene, not its own release |
| interconversion / modification pool | 122 / 19 | 107 / 14 | `DS_CYCLE_MODE=balance` (specs/039, 086, 091) | a protein's modified states are one conserved pool, solved as π = πP |
| self-contained regulator (TF:gene, sequestering inhibitor) | 82 / 17 | 11 / 11 | `DS_SELF_INHIBITOR_*`, `DS_SELF_ACTIVATOR_WEIGHT` (specs/022, 040, 083, 095) | the input the regulator is built from is counted once |
| recycle (p130/p107 back to E2F4 complexes) | 3 / 3 | 0 / 0 | `LNG_RECYCLE_POOLS` (specs/081) | uncurated return step written as a pool |
| **remnants (left after the rules)** | | | | |
| self-depleting ring (Goodwin) | 36 / 21 | 5 / 4 | numerically: LM fallback (specs/066, 073, 076) | no structural rule; the fixed point is found where relaxation fails |
| positive, closed by a catalyst | 20 / 17 | 37 / 13 | **none** | gain-1 knife-edge (specs/014) |
| positive, closed by material steps | 21 / 7 | 25 / **17** | **none** | gain-1 knife-edge |
| positive, closed through free ubiquitin | 13 / 10 | 3 / 2 | **none** (`DS_UBIQUITIN_MODE=inert` refuted, specs/098) | gain-1 knife-edge |
| other negative feedback | 1 / 1 | 3 / 2 | default fixed-point iteration | negative feedback has a stable fixed point |

Largest remnants: benchmark: Class I MHC antigen processing (890 nodes, ubiquitin-closed), DSB
Repair / DSB Response / DNA Repair (454 each, self-depleting), PIP3 (413, self-depleting). Sample:
Macroautophagy (337 + 254, material), G alpha (12/13) (314, catalyst), Hedgehog 'off' state (192,
material), neurotransmitter receptors (119, material).

## Pre-registered expectations

- **P1 holds:** touched 85% in the sample against 93% in the benchmark (8 points).
- **P2 holds:** the carrier rules touch 564 / 194 loops; the next mechanism (pools) 122 / 107.
- **P3 holds:** material-closed positive remnants occur in 17 sample pathways, more than any other
  unhandled class (catalyst 13, self-depleting 4).

## Reading

- The adopted loop rules are structural and fire outside the benchmark at a similar rate: carrier
  rules in 30 of 101 sample pathways, pools in 14, self-contained regulators in 11. This is
  coverage, not correctness: no ground truth exists for the sample.
- The benchmark over-represents self-depleting rings and ubiquitin-closed loops (signalling and DNA
  repair) and under-represents material-closed positive loops (17 sample pathways against 7).
- Unhandled by any structural rule: positive loops closed by catalysts or material steps (30 sample
  pathways together) and through free ubiquitin; self-depleting rings are handled only numerically.
  These are the targets for step 2, in that order by prevalence outside the benchmark: material,
  catalyst, self-depleting, ubiquitin. Each new rule is to be judged by a stated rationale and
  behaviour checks on networks without ground truth (convergence, no rail to 100× or 0 under a 2×
  input, label invariance, dose monotonicity), not by experimental cases.

## Amendment 1 (2026-10-25): held nodes

Tracing the sample's material-closed loops (specs/107) showed loops closed through H+, H2O, CoA-SH and
NAD+. These are on each bundle's cofactor list, and under the default `DS_COFACTOR_MODE=inert` the
solver holds them at baseline, so no signal goes round such a loop. The first run did not count that
mechanism. `loop_coverage.jl` now has a `held` mechanism: an edge into a cofactor, drug or pathogen
node (the default holds, specs/007, 055, 048) is handled. Re-run on the same builds:

| | benchmark (92) | sample (101), first run | sample (101), amended |
|---|---|---|---|
| touched | 690 (93%) | 306 (85%) | 319 (88%) |
| dissolved | 657 (88%) | 290 (80%) | 299 (83%) |
| untouched | 54 | 55 | 42 |
| loop nodes in dissolved loops | 46% | 51% | 53% |
| loop nodes left in remnants | 37% | 42% | 39% |
| held fires in | 1 loop | – | 15 loops, 10 pathways |

Sample remnants, amended: positive closed by catalyst 35 loops / 13 pathways, by material steps
20 / 15, through ubiquitin 3 / 2; self-depleting 5 / 4; other negative 1 / 1. The benchmark's
remnants do not change. P1–P3 still hold (P1: 88% against 93%; P3: material 15 pathways against
catalyst 13). The tables in this directory are the amended run; the figures above the amendment are
the first run's.
