# specs/168 research: the measurement (2026-11-20)

Arms run with `bench/analysis/168/sched_168.sh` at `91189f0b`, run name `lit168`, two at a time; log
`~/deltasignal-catalogs/analysis/168/sched.log`. Scorer `bench/analysis/168/score_168.py` (committed before any arm);
output `~/deltasignal-catalogs/analysis/168/score_168.json`. Before running, each variant's curated rows were checked
equal to its parent's `without` (no difference in any pathway).

**Not finished:** `pilot-with` and `pilot-single-with`. In each, the first curator solve did not return within the
registered 3,600 s (38,055 literature edges in PIP3 alone). They are recorded as not finished, not re-run.

## The gate: FAIL

v3-single against v3-without, primary pool: **net −374** (korkut −305, experimental −69). The net is ≤ 0, so per the
registration no shuffle was run.

## Per pair (fixed / broken on the cases scored in both arms; "prev" = readouts the base reached, "new" = newly covered)

| pair | axis | net | prev fixed / broken | newly covered | new change-call precision |
|---|---|---|---|---|---|
| **gate**: v3-single vs without | korkut | −305 | 72 / 51 | 401 | 0.009 |
| | experimental | −69 | 11 / 90 | 15 | 0.667 |
| | curator held-out (reported) | −355 | 37 / 119 | 461 | 0.022 |
| v3 vs without | korkut | −375 | 75 / 65 | 427 | 0.003 |
| | experimental | −11 | 8 / 29 | 15 | 0.667 |
| | curator held-out (reported) | −415 | 39 / 62 | 563 | 0.022 |
| **copy**: v3-single vs v3 | korkut | +99 | 118 / 16 | 3 | — |
| | experimental | −58 | 13 / 71 | 0 | — |
| | curator held-out (reported) | +113 | 158 / 45 | 3 | — |
| v3-strict vs without (in-sample) | korkut / experimental | −11 / −1 | 0 / 11, 0 / 1 | 0 / 0 | — |
| pilot-strict vs pilot-without (3 pathways, all tuning) | korkut | −81 | 51 / 105 | 36 | 0.000 |
| | experimental | −50 | 13 / 56 | 25 | 0.105 |

## The reading

1. **Attaching once roughly halves the damage on korkut and on curator readouts already reached, but does not
   make the edges useful.** v3-single against v3: korkut 118 fixed / 16 broken, korkut change calls 937 → 847.
   Against the curated network it is still −305 on korkut. Newly covered korkut readouts stay almost all wrong
   (precision 0.9%).
2. **On experimental, attaching once is worse** (−69 against v3's −11). 72 of the 90 already-reached cases it breaks
   are PIP3 over-expressions (ERBB2, PIK3CA, AKT1, KIT, …), which now read DOWN to about 0 where the curated network
   reads 50–100×. Traced on ERBB2 over-expression (trace servers at `91189f0b`): the two networks agree on every node
   at rest. The literature edge PHF8 → ERBB2 (into a curated root, present in both v3
   and v3-single) unroots ERBB2's free starting form, so it is no longer pinned (6 pins → 5). v3 still carries the
   signal through its other copies; v3-single does not, and 1,491 nodes diverge more than 10×. Across the 90 breaks,
   49 have the perturbed gene's starting form unrooted by a literature edge (v3: 19 of 29). The other 41 come through
   a non-local route in PIP3 not traced here.
   **Corrected by specs/170** (2026-11-23): with every curated pin restored the PIP3 breaks remain. They come from four
   literature edges PI3K ⊣ MTOR tipping a bistable component into its low state; the unrooting count above was a
   correlation, not the mechanism (specs/170 research.md).
3. **The strict selections move almost nothing on the 10 pathways** (v3-strict: 12 broken, 0 fixed). The pilot-strict
   bundle does reach new readouts on the three tuning pathways, but breaks more than it fixes on every axis.
4. **The dominant mechanism is still the specs/166 one.** A literature in-edge into a curated starting form removes the
   perturbation's pin, and edges attached to the readout region add change calls that are almost never measured
   changes on korkut.

## For the next bundle (sent to the literature source)

- No edge into a curated starting form, even when it is the only copy. Or a bench rule that ignores literature
  in-edges when deciding roots, which would be a protocol change needing its own registration and Adam's approval.
- Attaching once is right for AND slots, but a single picked copy can miss the perturbation's route. A pool node that
  averages the copies (as the generator's `variant_pool`) is the untried alternative.
