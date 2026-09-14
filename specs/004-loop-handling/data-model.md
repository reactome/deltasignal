# Data Model: Handle loops properly, and compare honestly

Phase 1. Entities are analytical rather than persisted; the "storage" is the
existing catalog files and benchmark outputs.

## Component (strongly connected component)

A maximal set of nodes mutually reachable in the directed network as solved.
Computed per pathway network, not across the catalog.

| field | type | rule |
|---|---|---|
| `pathway` | string | catalog directory name |
| `nodes` | set of uuid | size > 1 to be a component of interest; size-1 components are not reported |
| `intra_edges` | int | edges whose source and target are both in `nodes` |
| `negative_edges` | int | of `intra_edges`, those with `pos_neg == "neg"` |
| `reaction_stids` | set of string | distinct `edge_reaction_id` over `intra_edges`; the empty string is a distinct value and is counted (it marks bridges with no owning reaction) |
| `nodes_per_reaction` | float | `len(nodes) / len(reaction_stids)`; undefined and reported as such if `reaction_stids` is empty |
| `edge_type_counts` | map | `intra_edges` grouped by `edge_type` |

Current catalog: 34 components larger than one node, holding 4,561 nodes.

## Component classification

| class | rule | current count |
|---|---|---|
| `recycling_artifact` | `nodes_per_reaction >= 15` | 12 components holding 4,205 of the 4,561 cycle-resident nodes (92%) |
| `candidate_feedback` | `nodes_per_reaction < 15` | 22 components holding 356 nodes (8%) |

**The threshold is a reported parameter, not a constant to hide.** Measured
sensitivity: thresholds 10 and 15 give an identical partition (12 / 4,205);
at 20 two components move across, giving 10 / 4,078. So the partition is
stable across 10–15 and mildly sensitive above that — it is **not**
threshold-independent, and any result depending on it must state the value
and this sensitivity. Every result that depends on it must state the value used and
show the sensitivity across that range. Polarity is recorded but is NOT part
of the rule: 21 of 34 components have no negative edge, yet ERBB2's 465-node
artifact has 121, so polarity alone misclassifies.

## Case cyclicity

Attached to each scored benchmark case, for the acyclic/cyclic split
required by FR-003.

| field | type | rule |
|---|---|---|
| `readout_in_scc` | bool | any uuid the readout resolves to is cycle-resident |
| `perturbation_in_scc` | bool | any perturbed uuid is cycle-resident |
| `largest_scc_on_path` | int | largest component containing any node on a shortest signed path from perturbation to readout; 0 if none or no path |

A case is **cyclic** if `readout_in_scc` or `largest_scc_on_path > 0`.
Reporting uses this, not pathway-level cyclicity, because a pathway can be
36% cyclic while a given case never touches a cycle.

## Intervention arm

| field | type | rule |
|---|---|---|
| `name` | string | e.g. `drop_diagram_bridge` |
| `catalog_dir` | path | a derived catalog; the baseline catalog is never mutated |
| `edges_removed` | int | reported always |
| `cycle_nodes_before/after` | int | reported always |
| `destroys_curated_causality` | bool | true if it removes `catalyst`, `regulator`, `depletion`, `input` or `output` edges; such an arm is diagnostic-only under constitution principle I |

## Arm result

Every arm reports all of these; a result missing any field is not reportable.

| field | rule |
|---|---|
| `correct`, `total`, `accuracy`, `macro_f1` | on the scored subset |
| `per_class_f1` | DOWN / NO_CHANGE / UP |
| `per_pathway_net` | gained minus lost, per pathway |
| `changed_both_converged` | changed predictions where both arms converged — the noise guard |
| `coverage_delta` | cases scoreable in baseline but not in this arm, and the reverse. **Reported as coverage loss, never absorbed into accuracy** (FR-007) |
| `mcnemar_p` | paired test against the baseline arm |

## Control arm (MP-BioPath network)

| field | rule |
|---|---|
| `source` | `mp-biopath-pathways/pathways/<name>.tsv` |
| `columns` | parent dbid, child dbid, polarity (`1`/`-1`), conjunction (`0` = AND, `1` = OR) |
| `node identity` | Reactome database identifier — the same key the benchmark cases already use, so no translation layer is required |
| `common_case_set` | cases mappable in **all three** arms; a case unmappable anywhere is excluded everywhere |

Only 1.8% of MP-BioPath's edges are negative (449 of 25,196), against a much
higher share in ours. This is a second structural asymmetry and must be
reported next to the acyclicity one, not silently.
