# specs/080: do the `no_path` cases have a route in the database?

**Status:** analysis, 2026-10-20. No rule is measured here.

Adam, 2026-10-20: "I think that way more of the no_path cases should have a path and we
are just doing something wrong… the best check is with the database itself", and the
publication's list is to be checked against but not treated as the gold standard.

## Method

`db_reach.py` (this directory) runs every valid benchmark case of
`results/3e96e7c/dflt` (canonical `20261017-1246_5e663d9`) three ways:

- **Database (Neo4j, Release 97), scoped to the pathway's own reactions** (`hasEvent*`):
  - start from every participant that contains the gene's product;
  - step from an entity to a reaction where it is an input, catalyst or regulator (either
    sign), and from a reaction to its outputs;
  - the readout is reached when a reached entity is, or contains, the readout entity.

  Graded readings:

  | reading | rule |
  |---|---|
  | v1 | reactions only |
  | v2 | v1 + set membership both ways |
  | v3 | v2 + a component reaches a complex containing it when no reaction in the pathway makes that complex (a root complex) |
  | v4 | v2 + every complex containing it (loose) |
  | **v5** | v1 + member → set one way + the v3 root-complex step (**the strict reading used below**) |

- **Our network:** forward over every edge, from the gene's network copies to the
  readout's forms. This is the benchmark's own `no_path` test, run on every case.
- **MP-BioPath's publication** (`analysis_results/input_to_keyoutput_distance.tsv`,
  10 pathways): whether their network had a distance for the pair.

## Results (`summarise.py reach_3e96e7c.tsv`)

| | valid | our network reaches | database v1 | v5 | v3 | v4 |
|---|---|---|---|---|---|---|
| curator | 23,511 | 9,207 | 8,992 | 10,312 | 11,590 | 11,658 |
| experimental | 845 | 679 | 629 | 713 | 725 | 727 |

Under v5:
- **2,254 curator cases have a route in the database and none in our network**; 1,149 are the reverse.
- Of the 1,316 curator `no_path` errors, **392 (30%) have a database route**. The other
  924 have none within the pathway, even under v4 (438 reachable).
- Experimental: 60 database-yes/ours-no; **36 of the 124 `no_path` errors** have a route.
- MP-BioPath's networks had a path for 708 of the 845 experimental cases (ours 679,
  database v5 713). For 63 cases they had one and we do not.

So the hypothesis holds for about a third of the `no_path` errors: the database has the
route and our network loses it. For the other two thirds the pathway's own reactions
contain no directed route; that half of specs/079's ceiling stands.

## Where our network loses the route (`breaks.py`, `summarise_breaks.py`)

For each database-yes/ours-no case: the shortest v5 route, and the first node on it that
our network has but does not reach. Nodes our network does not have, such as expanded
sets and merged reactions, are passed through.

| step where it breaks | cases | wrong |
|---|---|---|
| a reaction edge (the node exists, but not from the copy we reach) | 1,342 | 81 |
| the readout form (the readout entity is reached, but not the forms scored) | 542 | 103 |
| a component into a root complex that contains it | 430 | 244 |

Most of the 1,342 reaction-edge breaks are cases the curators call NORMAL, and we are
right. The database has a route the curators did not expect to matter. Restoring routes
wholesale would add false changes, as every connectivity addition has so far
(specs/016, 019, 056).

**The largest single break:** Mitotic G1, free RBL2 (p130, R-HSA-389102) does not reach
RBL2:E2F4,E2F5:TFDP1,TFDP2 (R-HSA-1226089). That is **166 cases, 160 wrong**. The same
break for RBL1 (p107, R-HSA-389104 → R-HSA-1226088) adds 21 wrong.
- PP2A dephosphorylates p-RBL2 back to RBL2. In the database no reaction makes
  RBL2:E2F4/5:DP, so it is a root complex, and the generator decomposes it into its own
  root copy of RBL2. The RBL2 that PP2A produces never feeds it. The return of p130 to
  E2F4, which closes the pocket-protein cycle, is the step Reactome does not curate as a
  reaction.
- The route then runs through the DREAM complex to E2F and DREAM targets (TOP2A, cyclin A,
  CDC6 …).
- The curators expect it in every case: CCND1 and CDK4 knockout DOWN; CDKN2A, CDKN1B and
  PPP2R1A knockout UP; the reverse for overexpression. We predict NORMAL in all 187.
  13 of them are experimental cases.

## Next

A rule for the component → root complex break: a produced copy of a protein also
supplies the root complex that contains it. That means the generator's boundary leaf,
or a one-way supply edge like `DS_RECYCLED_TWIN`. It covers 430 cases (244 wrong), 181
of the wrong ones in Mitotic G1. It closes cycles, as the p130 one does, so it has to be
pre-registered with the usual gate. It is not measured here.
