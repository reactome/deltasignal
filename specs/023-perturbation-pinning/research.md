# 023 — Pin where the gene enters the network, and let the rest propagate

## What the benchmark does today (measured 2026-09-25, build `20260925-1039_d4f4f64`)

A perturbed gene resolves to the entities that reference it (Neo4j
`referenceEntity`, `neo4j_gene_to_stids`). **Since 2026-07-14 (`0bd4565`)** each
of those stids also matches every node whose `member_leaves` include it, so the
benchmark pins **every complex that contains the gene, wherever it sits**:

- **Scale:** 856 perturbations pin 9,378 nodes. Only **1,023 (11%) are roots**
  (in-degree 0).
- **Mid-pathway pins:** 8,355, made up of 6,371 complexes, 1,333 dissociation
  sinks and 651 entity forms produced by a reaction (e.g. a glycosylated
  secreted form).
- **Coverage:** 805 of 856 perturbations pin at least one mid-pathway node. 77
  have no root occurrence at all.
- **Example:** WNT5A in Signaling by WNT pins 1 root entity and 36 downstream
  nodes. These include generic set complexes (`WLS:WNT`, `WIF1:WNT`), and
  pinning them at 0 knocks out *every* WNT ligand. That is the specs/021
  reversal.

A hard pin overrides the node's own producing reaction. A downstream complex
is then not *computed* from the perturbed gene; it is *set*, including for
partners the knockdown should not affect.

## Intended protocol (Adam, 2026-09-25)

*"the idea is to pin the root input that is either a particular entity or a
complex containing that entity … breaking apart the complexes into root
entities where the root would look like A entity and B entity -> AB complex ->
reaction. In this way both the entity that is not part of the complex and the
entity inside the complex can be perturbed."*

## The rule (`DS_PIN_SCOPE=entry`; default `all` = today)

Let G be the nodes the gene resolves to today. Pin only the **entry
occurrences**: the members of G that cannot be reached from any other member
of G. Roots in G are always entries.
- A gene first produced mid-pathway (e.g. by transcription) gets that first
  occurrence pinned.
- When every member of G is reachable from another (G lies inside one cycle),
  all of G is pinned, as today.
- Everything else propagates.

## Pre-registration (committed before the arm runs)

- **Arms:** control = production scoring (`ae84de9`, byte-identical to the
  `ffb3aa1` control in specs/022); arm = `DS_PIN_SCOPE=entry`, code defaults
  otherwise.
- **Decides:** curator held-out net and macro-F1, McNemar p, pathways moved,
  distinct genes. **Reported:** tuning, experimental axis, per-pathway net,
  pinned-node counts, and the WNT5A dose ladder.
- **This is a protocol correction, not a model change.** The primary question
  is whether the benchmark measures what it claims to, so it is adopted if the
  held-out result is not worse than the noise floor (−15). Adam's statement of
  intent is the justification, not the accuracy. A gain is reported, not
  required. A loss beyond the noise floor is traced before any decision.
- **Prediction:** large churn, because 94% of perturbations change their pin
  set. Direction unknown. The false-change rate should fall, because
  downstream complexes stop being set directly.
