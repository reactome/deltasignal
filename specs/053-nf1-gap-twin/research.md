# specs/053 research: why an NF1 knockout cannot reach RAS

**Status:** trace only (2026-10-13). No rule is registered. This was meant to
register "an enzyme knockout raises its substrate" (specs/049 open item,
specs/052 R4); the trace shows the solver already does that once the enzyme
moves, and that the NF1 enzyme never moves.

## T1. NF1 does nothing under the default, in either direction

Canonical results `1f58a1b` on build `20261012-1401_01206ca`, RAF/MAP kinase
cascade (R-HSA-5673001), six readouts:
- NF1 knockout: curators and experiment say UP; predicted NORMAL, `pred_ui`
  exactly 1.0 (6 curator, 5 experimental).
- NF1 overexpression: curators say DOWN; predicted NORMAL, exactly 1.0 (6).
- Under specs/052 `hold`, the knockout reads UP (41-100x, one DOWN), the
  overexpression is unchanged.

## T2. The pin reaches 19 nodes and no readout

`root_cycle` pins free NF1 (R-HSA-947770, uuid e0eb0581). Everything it
reaches is SPRED dimer:NF1 → ubiquitination (R-HSA-5658424) → proteasome
(R-HSA-5658430). None of the six readouts is reachable.

## T3. The RAS GAP is an identity twin the pin cannot reach

Curated in Neo4j (release 97):
- R-HSA-5658438 SPRED dimer binds NF1 → SPRED dimer:NF1 (R-HSA-5658232);
- R-HSA-5658435 RAS GAPs (set R-HSA-5658216, which has SPRED dimer:NF1 as a
  member) bind p21 RAS:GTP (R-HSA-109783);
- R-HSA-5658231 RAS GAPs stimulate RAS GTPase activity → RAS GAPs +
  p21 RAS:GDP.

In the network, SPRED dimer:NF1 exists twice for each SPRED variant, with
`nodes.csv` rows equal in every column but the uuid:

| copy | made by | used by |
|---|---|---|
| 22bcec19 / 303d14d0 / a4b32081 | R-HSA-5658438 (from NF1) | ubiquitination only |
| 3152e339 / 920c10b3 / 598a10a8 | R-HSA-5658231 release only | the GAP binding RAS:GTP; a cycle-balance pool carrier |

The GAP copies have no producer except their own release step, so they
stay at baseline whatever happens to NF1. Under `hold` they are held at 0
through `member_leaves`, which is why `hold` got the knockout right: it
bypassed the split, it did not model the enzyme.

## T4. How narrow the shape is

- **Twins in general are not the lever.** `bench/analysis/twin_silo_census.py`:
  18,915 of 23,511 curator cases (787 of 845 experimental) reach some node
  whose identity twin they do not reach; the build has 104,572 twin nodes.
  That is the variant-copy structure the UUID-silo work closed at 1:3.9
  against.
- **The NF1 shape is unique.** A cycle-balance pool carrier whose only
  producers are its release steps, with an identity twin that has other
  producers: 3 carriers in the whole canonical build, all SPRED dimer:NF1
  in RAF. The cases whose pins reach the twin but not the carrier are the
  17 NF1 cases (12 curator, 5 experimental), all wrong today.

## T5. Options (none chosen)

1. Solver: a closed pool carrier follows its produced identity twin, as
   `closed_twin` (specs/049) does for carrier loops from root copies.
2. Generator: join the release copy to the produced copy where the
   identities are equal and the release copy has no other producer.
3. Leave it. Either fix is designed on, and could only move, these 17
   cases in one tuning pathway (RAF is one of the paper's ten), so its
   own score is not evidence for it; it rests on faithfulness (the GAP copy
   *is* the SPRED dimer:NF1 that R-HSA-5658438 makes).
