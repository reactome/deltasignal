# specs/053 research: why an NF1 knockout cannot reach RAS

**Status:** trace (2026-10-13). The rule it led to is registered in spec.md (generator join, Adam's choice). This was meant to
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

## R1. Result: not adopted (2026-10-13)

**Runs.** Both by `scripts/run_arm.sh` at **a807851** (solver identical to
main 89e8d28).
- `cj0`: canonical build `20261012-1401_01206ca`. Reproduces `1f58a1b` on
  every prediction and every `pred_ui` (24,100 curator, 849 experimental).
- `cj`: build `20261013-1452_cf3fced_ctjoin` (generator cf3fced, LNG PR
  #112, `LNG_CLOSED_TWIN_JOIN=1`).
- Results: `builds/<build>/results/a807851/{cj0,cj}`.

**Structural (J-P1): held.** 461 joins, with every pathway's count equal to
the census. Nodes 213,171 → 212,710 (−461); edges unchanged, because each
closed copy and its twin had disjoint edges.

**Accuracy: fails the gate on every line.**

| | fixed | broke | net | p |
|---|---|---|---|---|
| curator held-out | 40 | 533 | **−493** | < 1e-4 |
| curator tuning | 13 | 120 | −107 | < 1e-4 |
| experimental | | | +2 | |

By pathway (curator): DSB repair −488 (18/506), HDR −56, TP53 −47,
RUNX2 −30, TGF-β −8, IL-2 family −11. Three pathways lose more than 10.

**The losses are far outside the registered "touched" set (382 curator
cases).** That census counted cases whose pins reach the fed twin but not
the closed copy. It missed two effects of a merge, which is two-way:
- the fed twin now also reads the cycle's release, so any perturbation
  that reaches the cycle moves the twin and everything downstream of it;
- `root_cycle` pin selection changed: 176 DSB cases, 36 TP53 and 14 VEGF
  pinned a different number of nodes, because a merged node is no longer
  a closed cycle form.

The DSB loss is not traced case by case; both mechanisms are present
there.

**NF1 (J-P2, J-P3): both missed, and the trace shows why.** On the joined
build the GAP does follow NF1: solved directly, SPRED2 dimer:NF1 reads 1.0
(control), 6e-5 (NF1 knockout) and 28× (NF1 overexpression).
- **Knockout, 0 of 11 UP.** RAS:GTP (R-HSA-109783) reads 1.000007. NF1 is
  one of ten members of the RAS GAPs set (R-HSA-5658216), and losing one
  barely moves it. Not traced into the pool's drive.
- **Overexpression, 1 of 6 DOWN.** RAS:GTP itself reads 0.19, DOWN, which
  is correct. The five p-MAPK and RAS:RAF readouts read UP at 41–100×: the
  RAF-loop inversion recorded in specs/049 (P7/P8). So the join fixed the
  edge it targeted, and the readouts are decided further downstream.

**Conclusions.**
1. A uuid merge is two-way, and here that is the damage. It is the weld
   family of specs/018 again. A one-way form (the closed copy reads its
   twin; the twin does not read the cycle) is what `closed_twin` does for
   root copies (specs/049), and it would not have moved the twin's
   downstream. It is not registered.
2. The NF1 cases would still fail: the knockout because NF1 is one of ten
   GAPs, and the overexpression because of the RAF inversion. The NF1
   "enzyme knockout raises substrate" item is therefore not one defect, and
   the specs/052 gain on it was an artefact of holding every NF1 node.
3. LNG PR #112 is left unmerged (draft); the flag would expire under the
   specs/020 policy, and the branch is the record.
