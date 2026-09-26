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

## Results (build `20260925-2326_d4f4f64`, Reactome 97)

**New baseline** (`d68f6d9/baseline`): root_cycle pins, the corrected list,
code defaults.
- The first run (`8afaab8/baseline`) silently dropped IL-2 family. The
  benchmark took pathway ids from MP-BioPath's `pathway_list.tsv`, which still
  said 447115, found no such directory in the new catalog, and skipped it.
  `apply_catalog_ids` now takes ids from `bench/catalog_pathways.tsv`, the list
  the catalog is built from.
- The two runs are byte-identical outside IL-2 family.

| | cases | accuracy | macro-F1 |
|---|---|---|---|
| curator held-out (71 pathways), every case | 19,000 | **85.12%** | **0.8030** |
| curator, all 82 scored pathways | 24,100 | 82.56% | 0.7805 |
| experimental (10 pathways) | 849 | 65.61% | 0.5736 |

- **IL-2 family** now scores **247 / 260 (95%)**. On the old build all 260
  were unscorable, because the IL-12 network contains none of the genes.
- Experimental is 0.8pp below the same protocol on the previous build
  (66.43%), which is within its regeneration noise floor (15 cases, 1.8pp).

**Arm S, member-specific regulators** (against `8afaab8/baseline`; IL-2 is
absent from both sides of this comparison):
- held-out **+1** (1 / 0); tuning +15 (WNT, 2 perturbations); experimental +2.
- It does not clear the pre-registered +15 held-out bar, so it is **not
  adopted as a measured improvement**.
- It broke 2 cases in 23,176. The decision to adopt it on faithfulness
  (Adam's rule) is his.

**Arm D1, depletion removed: kept.**
- held-out −8 (115 / 123); tuning **−88** (TP53 −90, PIP3 −42, MET −40,
  IFN-γ −34); experimental **−45** (36 / 81). It fails the experimental
  condition.
- Traced: 425 of the 436 breaks are not at the readout. They are lost
  propagation, 258 DOWN→NORMAL and 138 UP→NORMAL.
- The edge encodes "an enzyme uses up its substrate": an E3 ligase knocked out
  leaves its substrate UP. Reactome draws only the substrate as an input to
  the ligase's reaction, so without the edge nothing upstream of the substrate
  can see the enzyme.

**Arm D2, dissociation removed: kept.**
- held-out **−72** (8 / 80); tuning −64; experimental −15.
- Traced: **all 152 breaks** are readouts fed by a dissociation edge, and all
  go DOWN→NORMAL.
- Curators name free subunits as readouts while Reactome outputs them inside
  complexes (the specs/015 released-subunit handles). The dissociation edge is
  their only route. MP-BioPath's published networks had the same complex →
  entity edges (3.2%).

**On visualising them.** All three derived types sit between entities that
are already on the diagram:
- `assembly` and `dissociation` connect a complex to its own components;
- `depletion` connects an enzyme to a substrate of its own reaction.

They can be drawn as annotations on existing glyphs rather than as new arrows.

## Corrections after review of PR #73, and the canonical baseline (`184dfd5/baseline`)

**Defect fixed: dissociation sinks were being pinned.**
- The recycled-root rule pinned dissociation sinks: out-degree-0 nodes fed only
  by dissociation. For 7 perturbations (ATR, ERCC2, ERCC3, RAC1, RHOB, SH3GL1,
  LMNA) those were the only pin, so 66 curator cases counted as scored while
  nothing moved.
- A node with no out-edges is no longer a root. "Recovers 75 of 76" was
  effectively 68 of 76. The newly-scorable accuracies quoted in specs/024
  included those inert cases; the every-case figures did not.
- Wording: the rule is "every producer is downstream (same strongly connected
  component) or arrives by a derived edge". That covers catalytic cycles, but
  not only them: 12 of the recycled pins are in no cycle, qualifying only
  through the derived-edge exemption. "0 of the other 780" is 0 of 801.
- specs/024 cites `geneToRootNodes`, which is dead code in mp-biopath. The live
  path is `IdMap.getIDmap` → `Evidence.getGenomic` → ROOT-only in nlmodel.jl.
  It agrees with the spec's mapping on all 877 perturbations.

**Gene names corrected** (`GENE_NAME_CORRECTIONS`), only where the intended
gene is unambiguous and in the pathway: PRKDC1→PRKDC, FOXM→FOXM1,
FBX7→FBXW7, CAP9→CASP9. Against the previous baseline: +44 curator cases (44
fixed, 0 broken) and +5 experimental.

**Release-skew exclusions** (Adam: tests whose entities are not in this
release cannot measure the generator or the solver). Each excluded case is
tagged in the dump:
- `gene_not_in_pathway`: 274
- `readout_not_in_release`: 124 (+3 experimental)
- `readout_not_in_pathway`: 54
- `gene_name_unresolved` (CACNAD1, IQGAP): 23

Cases that fail for our reasons, or Reactome's, are kept:
- MIR675: no root form; its root is the host gene H19.
- IFNA1: Reactome annotates the IFNA1 gene entity with IFNA13's Ensembl id.

| | cases | accuracy | macro-F1 |
|---|---|---|---|
| **curator held-out, in release** (70 pathways) | 18,573 | **85.82%** | **0.8133** |
| curator held-out, every case (71) | 19,000 | 85.28% | 0.8053 |
| curator, all pathways, in release (81) | 23,625 | 83.20% | 0.7896 |
| experimental, in release | 846 | 66.43% | 0.5803 |

**Other corrections to the text above:**
- Counts on build 2326: 4,587 regulator edges; 9,680 assembly, 9,738
  dissociation, 1,199 depletion.
- D1: the 425 non-readout breaks include 34 UP→DOWN sign flips, not only lost
  propagation. The E3-ligase explanation is a mechanism, not a traced break.
- D1 and D2 were measured against a baseline without IL-2 family, so IL-2 is
  unmeasured in those arms.
- IL-2 alone accounts for +81 cases of the change from the old build; the rest
  comes from root_cycle.
