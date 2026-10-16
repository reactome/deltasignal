# specs/065: per-pathway signature cache (generator bug fix)

**Status:** pre-registered 2026-10-16, before the build. From the week
review (docs/reviews/2026-10-15-week-review.md, finding 12).

## The defect

`src/variant_emission.py` caches each set member's reference signature in
the module-level `_sig_cache`. `prime_entity_caches` clears every other
per-pathway cache between pathways, but not this one. A catalog build runs
all 92 pathways in one process, so a member that appears in two pathways
keeps the signature computed in the first, from that pathway's
prefetch-limited entity caches (a getter returns None or "no components"
for anything outside the prefetch). The signature then depends on build
order, not on the pathway being built.

## Fix

Clear `_sig_cache` in `prime_entity_caches` (LNG). Test: the cache is empty
after priming.

## Measurement and prediction

One variant build from LNG main `da08370` + the fix, benchmarked with solver
main `41f3a5a` against canonical `20261016-0006_da08370/results/0406a4b`.
Since specs/013 the solve is label-free, so identical network content gives
identical predictions.

- **Prediction:** most likely nothing moves. A stale signature needs a member
  whose leaves resolved differently in an earlier pathway's prefetch.
- **If nothing moves** (every prediction identical, per-pathway node and edge
  counts identical): merge, no new canonical build needed.
- **If something moves:** trace every moved case; the fix is a correctness
  fix (build order must not decide a network), so it is adopted unless a
  traced loss shows the per-pathway signature is the wrong one. Report
  held-out, tuning and experimental net with McNemar p, both with and without
  RAF.
