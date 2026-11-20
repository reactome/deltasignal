# specs/166 research: literature edges (FullTextMining composite-h12b)

## Verdict: FAIL

The literature edges make predictions worse on every axis.

| | |
|---|---|
| net(with) on the primary pool (korkut + experimental) | **−461** (condition 1 needs > 0) |
| finished shuffles | 4 of 19: −917, −965, −1,063, −1,103 |
| beats every finished seed | yes |
| rank p with 4 seeds | 0.20 (registered bound ≤ 0.05, needs 19 seeds) |

Real edges do less damage than random ones with the same degrees and signs. So they carry
information, but at full strength they harm a curated network.

## Provenance

| | |
|---|---|
| Arms | `builds/20261118-lit-h12b-*`, built by `bench/analysis/166/build_arms_166.py` from the composite-h12b bundles (`external/literature-h12b/`, manifest verified) |
| Solver and bench | commit `17347537` (the pre-scoring subset-catalog fix, `bench/test_subset_catalog.py`) |
| `with`, `without`, slices | `results/17347537/lit166/` (curator, experimental, korkut) |
| shuffles | `results/17347537/lit166p/` (experimental, korkut; amendment 1) |
| Scorer | `bench/analysis/166/score_166.py`; output `analysis/166/score_166.json` |

## The three numbers (`with` vs `without`, cases scored in every primary arm)

| axis | cases | unreached in `without` | gained | lost | new change calls | their precision | new readouts' truth | previously covered fixed / broken | accuracy without → with |
|---|---|---|---|---|---|---|---|---|---|
| korkut | 2,811 | 2,146 | 427 | 65 | 388 | **0.26%** | 424 not / 3 changed | 58 / 90 | 0.785 → 0.637 |
| experimental | 493 | 80 | 15 | 12 | 15 | 67% | 15 changed | 20 / 74 | 0.726 → 0.637 |
| curator held-out (SELECT) | 5,692 | 4,141 | 560 | 45 | 409 | 2.2% | 547 not / 13 changed | 26 / 75 | 0.952 → 0.875 |
| curator tuning | 2,878 | 1,150 | 214 | 28 | 186 | 3.2% | 200 not / 14 changed | 120 / 233 | 0.845 → 0.746 |

Gains on the two lists handed over on 2026-11-18:
- **korkut:** 3 of the 108 scored pairs gained a route.
- **curator/experimental:** 789 of 5,376 gained a route.

### What happens on korkut

The new routes reach readouts that were measured as unchanged (424 of 427) and call them changed.
- 180 of the 427 are railed at the 0.01× or 100× bound.
- Previously covered calls are broken more often than fixed.

### What happens on experimental

The 15 newly covered cases are all measured changes, and 10 are called right. But on cases
already covered, the edges break 74 and fix 20.

## Cases lost: edges into curated roots (not in the nets above)

| axis | scored in `without` | scored in `with` |
|---|---|---|
| korkut | 4,428 | 2,836 |
| experimental | 545 | 508 |
| curator held-out | 7,622 | 5,692 |

The bench pins a gene at its root forms. A literature in-edge raises a root's in-degree, so the
root stops being one: 201 of about 4,900 curated roots across the 10 pathways. That lands on
exactly the genes that are perturbed.
- **korkut:** 1,245 (pathway, drug set) cases have no root left to pin. Another 372 have no
  readout copy left, because the abundance rule (specs/134) reads an edge into a starting copy as
  making it.
- **curator:** the losses concentrate in Apoptotic execution (592), RUNX1 (588) and RUNX2 (440).

## Shuffles (degree- and sign-preserving)

| seed | primary net | korkut net | coverage gained | new change-call precision | previously covered fixed / broken |
|---|---|---|---|---|---|
| shuf08 | −965 | −857 | 1,088 | 1.1% | 107 / 126 |
| shuf10 | −1,103 | −999 | 1,124 | 1.4% | 67 / 130 |
| shuf12 | −917 | −820 | 1,072 | 1.2% | 64 / 81 |
| shuf13 | −1,063 | −959 | 1,202 | 2.8% | 70 / 101 |

Random edges cover about 2.5 times as many readouts as the real ones and break more.

Unfinished seeds:
- **shuf09:** reached the 6 h cap.
- **shuf11:** one korkut solve took more than 3,600 s.
- **shuf14:** stopped when the run was ended.
- **shuf01–07:** stopped under earlier rules, with partial curator logs only, in `lit166/`.
- **shuf15–19:** not run.

The run was ended at 21:14 on 2026-11-19 by Adam's choice (option 2): the verdict does not depend
on the remaining seeds, and the rest would have taken about a day. Shuffled networks solve slowly
(DSB Repair 5,393 s in one seed). Partial curator per-pathway accuracies from the stopped seeds also
put the shuffles below `with`, e.g. IL-4/13 at 51–72% against 75.8%. They are descriptive only.

## Ablations: `only-<slice>` (descriptive, net vs `without` on cases scored in both)

| slice | primary | korkut | experimental | curator held-out | curator tuning |
|---|---|---|---|---|---|
| verified | −224 | −166 | −58 | −101 | −143 |
| extracted | −133 | −118 | −15 | −309 | −99 |
| sign stated | **+40** | +42 | −2 | −252 | −59 |
| sign inferred | −316 | −309 | −7 | −172 | −209 |
| sign na | −77 | −8 | −69 | −7 | −168 |
| binding/other | −15 | −20 | +5 | −65 | −18 |
| regulation | −383 | −395 | +12 | −326 | −102 |
| expression | −9 | −9 | 0 | −12 | −18 |
| prior < 0.4 | −11 | −20 | +9 | −33 | +21 |
| prior 0.4–0.6 | −284 | −218 | −66 | −222 | −259 |
| prior ≥ 0.6 | −164 | −122 | −42 | −111 | −109 |

Three readings, none of them a finding:
1. **Where the damage comes from.** The "regulation" group (regulator edges into molecule nodes)
   and inferred signs carry most of it.
2. **Verification does not help.** The verified layer and the highest prior_confidence bin do no
   better than the rest.
3. **"Sign stated" is the one positive primary net.** It is one of 11 slices, so treat it as a
   lead.

## The 105-pair question (FullTextMining 002; corrected by them to 80 in 012)

**Method.** `bench/analysis/166/trace_105.py`, run on `without` without solving. It resolves
targets and readouts exactly as `korkut_bench.py` does. Output:
`analysis/166/trace_105_structural.tsv`.

**Result.** Of the 284 unreached korkut pairs in the 10 pathways, none has a curated directed path
from a target node to the nodes the model scores for that readout.
- **56, all total-protein readouts,** reach some form of the readout protein, but only bound,
  modified or released forms. The abundance rule reads starting copies, and binding or modifying a
  protein does not change how much of it there is. Moving these readouts needs synthesis or
  degradation, which the networks mostly do not encode (specs/135).
- **228 (143 total-protein, 85 phospho)** have no path to any form of the readout protein.

FullTextMining's corrected count:
- 80 pairs, of which 10 are length 0 (one complex contains both the target and the readout);
- the other 70 are all total-protein readouts.

This is consistent with the 56. Their check reads any copy, including literature nodes as
endpoints.

## What this means

1. **Coverage is not the bottleneck on drug data; precision is.**
   - Real edges reach 427 more korkut readouts, and the measurements say 424 of them do not move.
   - specs/165 found that change-call precision on korkut was already about 11%. Literature routes
     at full strength lower it further.
2. **Edges into curated roots cost cases before any prediction is made.** If literature edges are
   tried again, roots must be decided without them (a bench rule; Adam's decision), or the edges
   must not target roots.
3. **The edges beat random ones, so they carry information.** Using it would need attenuation
   (prior_confidence or a learned per-class strength, the specs/125 direction) or a stricter
   selection, such as stated signs. Each would be a new registration, judged on data not used to
   choose it.

Nothing is adopted. The served catalog and the solver defaults are unchanged.
