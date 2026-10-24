# specs/099: what the unhandled positive loops are

**Status:** analysis, 2026-10-23. No rule is proposed here.

Adam, 2026-10-23: "ok" to tracing the unhandled positive loops (specs/097), starting with the
catalyst-closed ones.

## 1. Catalyst-closed loops are small, and self-catalysis dominates (`closing_catalysts.jl`)

Every catalyst or positive-regulator edge inside an unhandled positive component that has no
cycle without such edges, classified by how the catalyst E relates to the reaction R:

| | benchmark (92) | sample (101) |
|---|---|---|
| catalyst, **self** (E is also an input of R) | **400 edges, 137 reactions, 5 pathways** | 15 edges, 7, 3 pathways |
| catalyst, regenerated (E among R's outputs) | 29 edges, 2 reactions | 1 |
| catalyst, feedback (a different entity) | 17 edges, 13 | 2 |
| regulator, feedback | 22 edges, 3 | 7 edges, 7 |

Tracing them exposed an audit error (specs/097 correction): pairs carrying both an input and a
catalyst edge were labelled catalyst-closed. Corrected, catalyst-closed loops are 12 components
in 9 benchmark pathways and 8 components (47 nodes) in 7 sample pathways. Self-catalysis squares a
complex's own fold; specs/012 measured de-duplicating it as harmful (−15 held-out, on September's
networks).

## 2. Material-closed loops (`material_cycles.jl`, `classify_cycles.py`, `listed.py`)

The shortest material cycle of each material-closed loop, classified by what moves around it:

| class | benchmark | sample |
|---|---|---|
| **interconversion** (a protein leaf in every entity on the cycle) | **23 cycles, 9 pathways** | 32 cycles, 4 pathways |
| exchange (material between different entities) | 9 cycles, 8 pathways | 10 cycles, 9 pathways |
| **metabolic** (small molecules only) | 0 | **15 cycles, 11 pathways** |
| self-loop (a reaction outputs its own input) | 0 | 2 cycles, 2 pathways |

The interconversion cycles split again:
- **listed open carrier loops the solver keeps** (a displacer drives the release, specs/077's
  `fed_single`): benchmark 15 cycles in 4 pathways (intrinsic apoptosis R-HSA-109606, RIG-I/MDA5
  in R-HSA-168928, R-HSA-512988, Cell Cycle Checkpoints); sample 1. Kept on purpose: dropping
  them (`fed`) cost intrinsic apoptosis 26 cases;
- **listed in no table** (neither the pool nor the carrier-loop finder took them): benchmark 8
  cycles in 5 pathways (R-HSA-194138, R-HSA-5693571, R-HSA-6783310, R-HSA-74752, R-HSA-8848021);
  sample 31 cycles in 3 pathways (R-HSA-166058, R-HSA-75893, R-HSA-975155).

## Reading

The unhandled positive loops are not one class:
- **catalyst-closed:** small, mostly self-catalysis, a measured negative (specs/012) to undo;
- **interconversion, finder gaps:** a generator matter in a handful of pathways (5 benchmark, 3
  sample, the sample's concentrated in TNF and MyD88 signalling): the most promising structural
  target, a pool-finder extension like specs/081 and 086;
- **interconversion, kept open releases:** deliberate (specs/077), not a gap;
- **metabolic cycles:** common outside the benchmark (11 sample pathways) and absent from it, so
  no benchmark case can test a rule for them. Recorded as an unvalidated class for the external
  test (specs/096).
