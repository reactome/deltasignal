# specs/079: what is left to win, and where

**Status:** analysis only, 2026-10-20. No rule is proposed or measured here.

Adam, 2026-10-20: "do you feel like we are getting somewhere with the approach we are
taking". This spec answers that with a census of every error under the current
defaults: main `3e96e7c`, canonical `20261017-1246_5e663d9`, `results/3e96e7c/dflt`.
`census.py` (this directory) reproduces every figure.

## Where the errors are

| | valid | wrong | no route in the network (`no_path`) | routed but wrong |
|---|---|---|---|---|
| curator | 23,511 | 2,763 (88.2% right) | 1,316 | 1,447 (790 false change, 657 missed) |
| experimental | 845 | 244 (71.1% right) | 124 | 120 (36 false change, 84 missed) |

About half of all errors have no route from the perturbed gene to the readout **in the
network of the pathway being scored**, the one the benchmark solves. That is not a
claim about Reactome as a whole (below). No solver rule can fix those: the most any propagator can score on
these networks is **94.4% curator** and **85.3% experimental**. The solver work since
specs/039 has been spent on the other half.

## The experimental set and the curators disagree

All 845 valid experimental cases are also curated. The curators' answer matches the
experiment in **693 (82.0%)**.

- Where they agree, we are right on **592 of 693 (85.4%)**.
- Where they disagree (152), we match the curators 134 times, the experiment 9 times,
  and neither 9 times.

So **134 of the 244 experimental errors (55%) are cases where we give the curators'
answer**. One example is RBL1 and RBL2 knockouts in Mitotic G1 (20 cases): E2F targets
read UP, as the curators say. The experiment says NORMAL, plausibly because the pocket
proteins back each other up. A rule that moves the model towards the curated network cannot win those cases. It
is also why curator and experimental gains have traded against each other more than once
(specs/046, 049, 077).

What is left where the curators and the experiment agree is **101 cases**: 53 missed,
43 with no route, and 5 false changes. The largest groups:

| perturbation | cases |
|---|---|
| PIP3: PTEN overexpression | 9 |
| Mitotic G1: CCND1 overexpression | 6 |
| RAF: NF1 knockout (RAF is set aside, specs/054) | 5 |
| Mitotic G1: RB1 knockout | 4 |
| WNT: APC knockout | 4 |
| WNT: APC overexpression | 3 |
| Mitotic G1: AKT1 knockout, SRC overexpression | 3 each |
| Cell Cycle Checkpoints: MDM2 knockout | 3 |

## What this says about the approach

- **It has worked where it can work.** Since specs/049 (2026-10-12) curator held-out in
  release has gone from 89.29% to 90.26%, each step traced and gated. Every recent
  adoption has removed an artefact found by tracing cases.
- **The returns are small and shrinking.** Each spec moves 20–65 held-out cases out of
  18,573, and the routed errors that remain are spread over 427 curator perturbations.
- **The curator held-out split has been used to choose rules for six weeks** (specs/051:
  SELECT, not TEST). Several adoptions were post hoc (061, 077). The curator gains are
  real on these cases, but they are not a clean estimate of accuracy on new ones.
- **The experimental axis is where to judge progress, and it has not moved:** 591/849 at
  specs/048, 606 before specs/077, 601 now. Most of what it counts as error is a
  disagreement with the curated pathway itself.

## Recommendations

1. **Report the experimental axis split by curator agreement**: 85.4% where the curators
   agree with the experiment, plus the 152 disputed cases on their own. It is a fairer
   headline, and it separates model error from curation-versus-experiment disagreement.
2. **Aim rule work at the 101 agreed experimental errors**, starting with PTEN
   overexpression in PIP3. They are the only errors that every reference says are ours.
3. **The `no_path` half is mostly not missing biology.** A 2026-07-16 check, on an older
   catalog and not re-run here, found the mechanism in Reactome for nearly all such cases:
   - Of 229 experimental failures, 99% were reachable in the full human reaction graph
     when complex and set membership are followed both ways. Following reactions alone,
     62% were.
   - Most curator `no_path` cases were co-descendants of one hub (gene and readout both
     downstream of it), not gene → readout.
   - Importing neighbouring pathways' networks (specs/056) cost held-out −164.

   So this half needs either a wider network that does not over-couple, or a
   non-directed (co-regulation) reading. Neither is a curation request.
4. **Hold an untouched test.** Release 98, when published, scored once with frozen rules
   (specs/051, 062), is the only clean estimate of whether these six weeks generalise.
