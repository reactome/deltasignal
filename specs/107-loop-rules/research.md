# specs/107 results

Code: `src/core/material_loops.jl` (rules 1–3), `leaf_refs` on `ReactionNetwork` and
`parse_leaf_refs` (`src/io/tsv_parser.jl`), applied first in `solve_steady_state` and in
`solved_carrier_network`. Tests: `test/test_material_loops.jl` (44). Generator: `leaf_refs.csv`
(logic-network-generator PR #128, branch `leaf-refs-107`).

Catalogs: the reference table needs the generator change, so both catalogs were copied as overlays
with `leaf_refs.csv` added and nothing else changed: `20261024-2211_7f4b28a_refs107` (benchmark) and
`20261025-1318_7f4b28a_sample101_refs107` (sample). Each overlay's `BUILD.json` records its source build
and how the file was written.

## Implementation notes found while building

- A protein node is a `simple_entity` too; rule 1's "small molecule" therefore means every leaf
  is an `R-ALL-` entity (caught by the first test fixture).
- The specs/106 audit missed cofactor holds; that is specs/106 amendment 1, applied here.

## B1 coverage (specs/106 audit, `loop_coverage.jl` with the rules applied)

Sample (101 pathways, 361 loops):

| arm | dissolved | partly handled | untouched | loop nodes in dissolved loops | material-closed remnants |
|---|---|---|---|---|---|
| rules off | 299 | 20 | 42 | 53% | 20 loops / 15 pathways |
| rule 1 | 303 | 20 | 38 | 53% | 16 / 13 |
| rules 1 + 2 (`strict`) | 304 | 20 | 37 | 53% | 15 / 12 |
| rules 1 + 3 (`wide`) | 334 | 14 | 13 | 83% | 11 / 9 |

Benchmark (92 pathways, 744 loops):

| arm | dissolved | partly handled | untouched | material-closed remnants |
|---|---|---|---|---|
| rules off | 657 | 33 | 54 | 21 / 7 |
| rules 1 + 2 | 660 | 33 | 51 | 21 / 7 |
| rules 1 + 3 | 670 | 27 | 47 | 15 / 4 |

## Predictions

- **P1 fails:** rule 1 dissolves 4 of the 11 sample metabolite cycles (predicted at least 8). The
  rest are two-substrate reactions, which a one-in/one-out conversion does not cover: transaminations
  (GOT2, OAT), transketolase, GPX2 (GSH with H2O2), pyruvate carboxylase (pyruvate with HCO3-), an
  N-acetyltransferase with acetyl-CoA. Cytochrome c oxidised ⇄ reduced is a protein redox pair the
  generator's pool finder does not list; rule 3 reaches it through the reference table.
- **P2 fails:** `strict` changes no benchmark material-closed remnant. Traced
  (`why_not.jl`, per entity node on each remaining cycle):
  - the carriers of PAK (Signaling by Rho GTPases) and L1/microtubule pass through set-variant
    nodes, which `strict` does not walk and `wide` does;
  - tBID:BCL-2 is listed by the generator as an open loop and kept on purpose (BAD displaces tBID,
    specs/077);
  - the DDX58/IFIH1 TRAF3:MAVS cycle is listed as an open scaffold loop and kept on purpose
    (specs/077);
  - cyclin B1:CDK1 between nucleus and cytoplasm changes CDK1's modification (different leaves),
    which only the reference table joins.
- **P3 fails narrowly:** rules 1–3 leave material-closed remnants in 9 of the 15 sample pathways
  (predicted at most 8).
