# specs/043: rejoin copies of one entity that the generator split without cause

## Verified defect (specs/041 step 3b; specs/042 research)

The generator gives each occurrence of an entity its own node (uuid). Phase 2
then joins a producer's output to a consumer's input **only along curated
`precedingEvent` pairs** (`_get_or_create_entity_uuid`). That was deliberate.
- Merging every occurrence by stId welded loops Reactome does not contain:
  specs/018, and MP-BioPath's hand-cut loops (specs/034 §11).
- Per-position copies resolve about 86% of those.

**But where a producing reaction has no `precedingEvent` link to a consumer of
the same entity, the copies are never joined**, even when Reactome curates
them as one entity flowing from producer to consumer.

| case (verified in Neo4j) | Reactome | our build |
|---|---|---|
| WNT, CTNNB1 [cytosol] R-HSA-448839 | one entity: bound by the destruction complex (195304), released by it (201685), imported to the nucleus (201669) | 3 copies: a root that only binds, a release product that only goes to the nucleus, a dissociation dead end. The curated cycle is severed. |
| WNT, AMER1 [cytosol] | the same pattern | 3 copies |
| TP53, TP53 Tetramer R-HSA-3209194 | one entity: 3 producers, 21 consumers | 4 copies. The p14ARF (6804996) and USP7 (3215310) products have **no following events**, so they dead-end in dissociations. |

These copies account for 18 traced experimental errors, and an estimated
499 root-split stIds in 50 pathways (census, specs/042 derivation-fable §0).

## The fix, stated narrowly

**Join a produced copy P and a consumed copy C of the same stId, in the same
pathway, only when all of these hold:**
1. **P's producing reaction has no `precedingEvent` link to any consumer of
   that stId.** The curation is silent, not contradictory. Where a link exists
   to some other consumer, Phase 2 already decided, and it is left alone.
2. **C is a root copy** (nothing produces it), or a copy whose only producers
   are derived edges (dissociation, assembly). That is, C would otherwise be
   an unexplained source.
3. **The join does not create a strongly connected component that Reactome's
   own reaction graph lacks.** Every new cycle must exist among the curated
   entities and reactions (the specs/041 `loops.py` check). This is the guard
   against re-welding.

Behind a generator flag (`LNG_REJOIN_ORPHANS=1`, default off) until measured.

## What it does NOT do

It adds no solver rule. A rejoined binding/release cycle with no supply (WNT)
may still read wrong until a conservation or turnover rule handles it (specs/039
pools; specs/042 R1). This spec measures the rejoin alone first.

## Census before any arm

Count:
- joins per pathway;
- new strongly connected components, and their sizes;
- rejected joins (by condition 1, 2 or 3);
- cases whose route changes.

List every join in WNT and TP53 by name, for review.

## Gates

The same as specs/040, both axes:
- curator held-out > +15, p < 0.05;
- experimental ≥ 0, not significantly negative;
- no pathway loses > 10.

Plus:
- **zero new strongly connected components absent from Reactome;**
- relabel churn no worse than the control.

The rejoin is a faithfulness fix: it makes the network match Reactome. So
**neutral results on both axes are acceptable for adoption, with Adam's
agreement**, as variant sharing was in specs/020.
