# specs/178: Lun 2019 on Release 98

**Status:** registered 2026-11-26, before the Release 98 build or any Release 98 prediction was made. Run 2026-11-27: every specs/177 verdict holds (`research.md`).
**Requested by:** Adam, 2026-11-26 ("you should be using Release98 anyways").

## Why

Every Lun prediction so far (specs/127, 173–177) came from the whole-Reactome build
`20261027-1346_8b6527e_top114`, which is Release 97: it was built on 2026-10-27, before Release 98 was served
(specs/150–152), and no Release 98 build of the 28 top-level pathways exists. The tyrosine-kinase class table was read
from Release 97 as well. The served catalog is Release 98, so the Lun numbers should be too.

## What is run

1. **Build.** The 28 top-level pathways (`specs/114-whole-reactome/top_level_pathways.tsv`, all present in Release 98)
   with the generator at `d6d026b` (the commit behind the served Release 98 catalog) against Release 98 (checksum
   525378806530), as a variant build. `current` and the served catalog do not move.
2. **Predictions.** `run_arm.sh --axis lun_predict` on that build at main's HEAD: the specs/127 protocol with the
   abundance rule (specs/175). Values-blind, as before.
3. **Class table.** `classes_177.py` against Release 98, written to `bench/analysis/178/`. The readout table does not
   depend on the release.
4. **Scoring.** `score_177.py`, unchanged apart from taking the class table's directory, on both halves with the
   specs/173 measured call.

## Reading

Nothing is chosen here: the rule, the measured call and the gates are those of specs/177. **Half B was seen once in
specs/177 and is FIT**, so this is a re-measurement on a new release, not a second held-out test. The pre-registered
test of record stays specs/177.

- Reported for each half, Release 98 beside Release 97: the specs/177 gates (candidate − reachability, candidate −
  model, added calls right vs wrong), the arms, route enrichment and lift, and the class table.
- Reported: how the networks differ for these pairs (pairs scored, calls that changed between releases, and which
  readouts they are), and the class table's differences between releases.
- **Expectation:** the specs/177 verdicts hold on Release 98. If any verdict differs, both are reported and the moved
  pairs are traced before anything is concluded.
- If they hold, the Release 98 figures are the Lun figures quoted from here on, with the release named.
