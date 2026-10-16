# specs/065: per-pathway signature cache (generator bug fix)

**Status:** pre-registered 2026-10-16, before the build; measured the same day (research.md): no network and no predicted class changes. Merge LNG #115, no new canonical build. Amendment 1 (five more caches, LNG #116) measured the same day: again identical networks and no class change. From the week
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

## Amendment 1 (2026-10-16, before the build): every derived per-pathway cache

Reading `neo4j_connector.py` directly (week review gap) found the same defect
in five more caches. Each memoises values computed from
`get_complex_components` / `get_set_members` / `get_labels`, which return
empty for anything outside the current pathway's prefetch, and none is
cleared between pathways:

- `logic_network_generator._variant_leafsets_cache`
- `logic_network_generator._matching_leaves_cache`
- `logic_network_generator._handoff_leaf_cache`
- `logic_network_generator._sets_chosen_cache`
- `variant_keys` memo (`reset_caches()`: `_memo`, `_variants_cache`, `_struct`)

Caches filled by their own direct Neo4j query (`_containing_cache`, drug,
pathogen, form-profile, cofactor) hold database truth and are left alone.
Every catch-all in the connector logs and re-raises; none falls back silently.

Fix: clear all five in `prime_entity_caches`. Measurement, prediction and
decision rule as above (most likely nothing moves; if something moves, trace
every case; adopted unless a traced loss shows the per-pathway value is the
wrong one), against canonical `20261016-0006_da08370/results/0406a4b`.
