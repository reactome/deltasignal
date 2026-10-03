# specs/044: the boundary decomposition is a derived layer, not the curated network

**Decision (Adam, 2026-10-03):**
- The logic network ships what curators curated.
- The generator's decomposition of root-input and terminal-output complexes is
  written to its own file, `boundary_edges.csv`.
- DeltaSignal loads both files.
- Set semantics inside that layer (finding F7) are fixed there, and measured.

## Why

The generator decomposes every root-input complex into its components
(`assembly` edges, one `hasComponent` level at a time since specs/030) and
every terminal-output complex into readout sinks (`dissociation`). This is how
a single protein can be perturbed or read when a complex is what the curators
drew. But it is our inference, not curation, and it was mixed into
`logic_network.csv`, which causes two problems:

1. **Other users of the network get a different network from the one curated.**
   A root complex stops being a root: its curated inputs are replaced by
   protein leaves. Filtering on `edge_type` cannot undo this, because since
   specs/030 the same `assembly` edges also join complexes to nodes the
   pathway computes.
2. **F7 (LNG code review 2026-10-02).** The decomposition flattens a set
   inside a root complex into its members, and makes every member a required
   AND input.
   - Example: "VAV1 Rho/Rac effectors:GDP" (R-HSA-114543) is GDP plus a
     CandidateSet of RHOA, RHOG, RAC1, RAC2 and CDC42. In GPVI, a knockout of
     any ONE of the five drives the complex to 0, the same as knocking out all
     five (solved through the API, 2026-10-03).
   - Catalog-wide (build 20260928-1110_06ccb63): 3,828 of 10,165 assembly
     edges reach their complex only through a set. 557 complex nodes, in 81
     pathways, need two or more alternatives at once.
   - Counted from Reactome's own structure by
     `~/deltasignal-catalogs/analysis/044/f7_scope.py`.
   - Benchmark exposure: 3,920 curator cases (100 genes) and 206 experimental
     cases (20 genes) perturb such an alternative. They currently score
     **above** average (88.2% vs 85.8%; 78.6% vs 70.2%). So the fix is NOT
     expected to raise accuracy by default, and it must be measured, not
     assumed (`f7_cases.py`).

`containment.csv` cannot replace the decomposition as it stands. It flattens
to leaves, so it cannot tell a required subunit from one alternative of a set.

Why keep the decomposition at all: without it, a gene knockout reaches a
complex only through curated reactions. The specs/030 hierarchy (held-out
+197) works by joining a root complex's components to nodes the pathway
computes.

## Contract

The generator, logic-network-generator branch `feat/boundary-edges-file`:

| File | Contents |
|---|---|
| `logic_network.csv` | Edges from curated reactions only. |
| `boundary_edges.csv` | Same columns. The edges the boundary pass adds: `assembly`, `dissociation`, and `composition` when that switch is on. They are tagged where they are emitted, not recognised by `edge_type`: specs/036's capped-set pools also use `assembly`. |
| `containment_structure.csv` | Reactome's structure, one level at a time: `parent_stable_id`, `parent_type`, `child_stable_id`, `relation` (`component` / `member` / `candidate`), `stoichiometry`, `reactome_release`. |
| `stid_to_uuid_mapping.csv` | Covers the nodes of both edge files. |

DeltaSignal:
- `parse_logic_network` reads `boundary_edges.csv` beside `logic_network.csv`,
  appended after the curated rows. With the default generator switches this
  is the original row order (see Row order below).
- The columns must match, or the read fails.
- Every Python reader goes through `bench/network_files.py:open_network`.
- `catalog.sh` counts both files in `BUILD.json`, and treats a pathway whose
  boundary file has different columns as not built.
- An upload or a renamed copy cannot use the "beside logic_network.csv" rule.
  `/api/parse` therefore takes a `boundary_edges` part, the CLI `parse` takes
  `--boundary`, and `parse_logic_network` / `parse_complete_network` take
  `boundary_path` (`:sibling`, a path, or `nothing`).
- Row order (review of #89): in the served build every boundary edge is
  already at the end of `logic_network.csv` (checked: GPVI, NODAL, RAF). The
  diagram set-member and handoff passes, both off, are the only ones emitted
  after the boundary pass. So the split keeps the original row order. Two
  regenerations of one commit can still differ in row order (RAF, from row
  3,243). That is the pre-existing ordering nondeterminism (code review F14),
  present without the split as well.
- A bundle from before the split has no boundary file and reads as before.

Still inferred, and still in `logic_network.csv`: `depletion` edges and set-pool
nodes (specs/033). They are generator inferences too. Moving them is a separate
decision.

## Verification (before either side merges)

1. **Generator:** GPVI, NODAL and RAF regenerated from the branch. The union of
   the two files equals build 20260928-1110_06ccb63's `logic_network.csv`, edge
   for edge by stable id. `logic_network.csv` has no assembly or dissociation
   edges. Done 2026-10-03.
2. **DeltaSignal:** split the served build's `logic_network.csv` by
   `edge_type` into the two files (`builds/20260928-1110_06ccb63_split044`). For that build the boundary pass is exactly
   `assembly` + `dissociation`: composition and capped pools are off, and
   every assembly edge there has no `edge_reaction_id`. Score it with this
   branch, then compare against the unsplit build scored by the same code.
   **Required: every prediction identical on both axes.**

   **Result (2026-10-03, solver 8030591):** identical. All 24,100 curator cases
   and all 849 experimental cases match in every field, the numeric
   prediction `pred_ui` included. Both arms score 20,482 of 24,100 (84.99%),
   the canonical figure. The split moved 19,903 edges (10,165 assembly and
   9,738 dissociation) and left 194,294 curated. Results:
   `builds/20260928-1110_06ccb63/results/8030591/ctrl044` and
   `builds/20260928-1110_06ccb63_split044/results/8030591/split044`.
3. Merge DeltaSignal first, then the generator. A DeltaSignal that does not
   read `boundary_edges.csv` would silently lose those edges from a new
   catalog.

## Next

Fix F7 inside `boundary_edges.csv` only: the members of a set component feed
one OR node, which feeds the complex, as specs/033 does for set-valued
catalysts. Predictions are registered here, and the arm runs on both axes,
before anything is adopted.
