# 025 — Member-specific regulators, and whether our derived edges earn their place

Decisions and rules from Adam, 2026-09-25:

1. **Recycled roots become the protocol** (specs/024): *"mimicking this while
   keeping the loops for the test sounds like the most faithful way of running
   the tests."* `DS_PIN_SCOPE=root_cycle` is now the benchmark default.
2. **The IL-2 family entry is fixed** to R-HSA-451927. MP-BioPath's list had
   447115, which is Interleukin-12 family signaling; the published network
   holds 58 of 451927's 59 reactions and 0 of 447115's (specs/024). Every other
   entry was checked the same way: median 97% of each catalog pathway's v97
   reactions are in its published network. The only other zero is Rho GTPase
   cycle, whose published reactions no longer exist in v97 (a Reactome
   rewrite, not a typo).
3. **A regulator naming a specific set member** (CREBBP:NS1) *"should do it for
   every reaction it impacts"* if the member is known exactly. The generator
   attaches a Reactome reaction's regulators to every variant.
   `DS_SIBLING_REGULATORS=member_only` drops a regulator edge from a variant
   when the regulator names a member that only a sibling variant uses and none
   that this variant uses. That is 113 of 4,634 regulator edges, in 9 pathways.
   Traced example: CREBBP:NS1 is dropped from *EP300 + IRF3* and kept on
   *CREBBP + IRF3*.
4. **Edges our generator adds that are not in Neo4j.** *"They are not in the
   pathways in neo4j and will be hard to visualize on the pathway diagrams.
   And they are ruining some things. Unless you think there is good reason to
   keep them."*

   | type | edges | in MP-BioPath's published networks? |
   |---|---|---|
   | `assembly` (entity → complex) | 9,642 | yes: entity → complex is 18.9% of their edges |
   | `dissociation` (complex → entity) | 9,600 | yes, fewer: 3.2% |
   | `depletion` (catalyst −→ substrate) | 1,199 | **no counterpart** |

   `assembly` is the "A and B → AB complex" decomposition that makes a root
   entity pinnable at all, so it is kept. `dissociation` and `depletion` are
   measured by removal.

## Pre-registration (committed before any arm runs)

- **Build:** a new catalog from generator HEAD with the corrected pathway
  list.
- **Baseline:** code defaults with `DS_PIN_SCOPE=root_cycle`. It becomes the
  new canonical in docs/RESULTS.md, with the old build's numbers alongside.
- **Arms**, each paired against that baseline through `scripts/run_arm.sh`:
  - `S`: `--bench DS_SIBLING_REGULATORS=member_only`
  - `D1`: `--bench DS_SKIP_EDGE_TYPES=depletion`
  - `D2`: `--bench DS_SKIP_EDGE_TYPES=dissociation`
- **Reading rules:**
  - **S** (a model correction): adopt if held-out net > +15 with p < 0.05, the
    gain is not one pathway, and experimental is no worse than −15.
  - **D1 and D2** (removals Adam prefers for faithfulness and visualisation):
    adopt removal if held-out is **no worse than −15** (the noise floor) and
    experimental is no worse than −15. A larger loss means the edge carries
    information that has to be put back another way; it is traced before any
    decision.
  - All results are reported with pathways moved, perturbations, and both axes.
- **Known interaction:** the pin rules read root-ness from the unmodified
  network. Under `root_cycle`, dissociation and depletion in-edges already
  count as non-disqualifying, so removing them does not change who is pinned.
