# specs/181: can kinase–substrate knowledge cover Lun's missing routes, and what does the model get wrong where a route exists?

**Status:** registered 2026-11-27, before any kinase–substrate table was downloaded or read.
**Requested by:** Adam, 2026-11-27 ("we need to figure out why it is running so poorly and find solutions"; "ok to your
plan").

## Why

On Release 98 (specs/178 predictions, both halves pooled, model alone), 1,001 of 1,437 measured changes have no route
from the protein to the readout in any network, so recall is capped near 30% whatever the solver does. 212 more have a
route and the model is silent (125) or calls the wrong direction (87). These are separate problems with separate
parts below. Both halves of Lun 2019 are FIT: **everything here is description on seen data; nothing is adopted from
it, and anything worth adopting needs a dataset not yet seen.**

## Part K: enzyme–substrate edges as reachability

**Table.** OmniPath's enzyme–substrate interactions, human, as served on the download date (file and sha256 recorded;
the file stays outside the repo). No licence terms are accepted on Adam's behalf; OmniPath is used because it needs
none. Only phosphorylation and dephosphorylation rows are used.

**Arms** (each a set of protein → readout calls, added where curated reachability has none):
- **K1, exact site:** the over-expressed protein phosphorylates a registered site of the readout
  (`bench/analysis/127/readouts.tsv`, gene and position) → UP; dephosphorylates one → DOWN. Both on one pair → no call.
- **K1p, any site:** as K1 for any site of a readout's proteins (reported only; the antibody does not read those sites).
- **K2, one hop into the curated network:** the protein phosphorylates any site of a Lun protein Z, and Z has curated
  reachability to the readout → UP (direction is not known from a site, so UP is the only call; reported with that
  caveat). Z is limited to Lun's proteins because curated reach is only computed for them.

Self pairs and specs/179's direct pairs are handled as everywhere else (self pairs out; direct pairs kept, since
overlap with the curated first hop is part of the question and is reported).

**Measured, pooled and per half:**
1. Coverage: how many of the no-route real changes each arm reaches.
2. Precision of the added pairs: share measured changed, and share matching the arm's direction, against the base
   rate among all no-route pairs.
3. A relabelling null, as specs/113: enzymes' substrate lists permuted among enzymes of the same type (kinase or
   phosphatase) that are Lun proteins, 1,000 draws, seed 181. An arm is informative only if its added-pair precision is
   above the null's 95th percentile.
4. Balanced accuracy of (model + tyrosine-kinase rule + arm) against (model + tyrosine-kinase rule), protein bootstrap,
   seed 181, and the same without the tyrosine-kinase rule.

**Registered reading.** K1 "covers the gap" if it reaches at least 20% of the no-route real changes with added-pair
precision above the null's 95th percentile. Below that, kinase–substrate knowledge is recorded as not the missing
piece on this dataset, whatever balanced accuracy says. K1p and K2 are reported, not judged.

## Part M: the 212 route-exists errors, traced

Each of the 125 silent and 87 wrong-direction cases is classified from the network (specs/174's tracer on the Release
98 build) and the prediction file: route length, whether the route passes a pool, a switch component, an inhibitory
step, a total-readout rule, or a set readout where the mean dilutes one member; and the predicted fold against the
0.85 / 1.15 cutoffs. Classes are counted and one case per class is traced by hand. Analysis only; no rule is
registered here.
