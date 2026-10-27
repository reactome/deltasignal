# specs/113: Reactome interaction edges (census before any rule)

**Status:** analysis only, 2026-10-27. No rule is registered and nothing has been measured.

Adam, 2026-10-26: use the database's interaction records to add connections. Reactome 97 holds 78,863
`UndirectedInteraction` nodes (IntAct), each joining two or more `ReferenceGeneProduct`s through `interactor`. They
are unsigned and undirected, and they attach to proteins, not to the curated forms (EWAS, complexes) that the networks
are built from. HPN-DREAM is FIT after specs/112, so any rule would need a fresh test. The candidate test is LINCS P100,
which is TEST and still unseen. This census reads only P100's analyte table (probe → UniProt, site) and the graph.

## LINCS P100 readouts and what interaction edges could reach

| | count |
|---|---|
| probes / distinct probe proteins | 96 / 86 |
| probes whose phospho-site is a phosphorylated EWAS in Reactome 97 (specs/111) | 7 |
| probes whose site occurs in the 92 benchmark networks (specs/111) | 4 |
| probe proteins present in some form in the networks (`leaf_refs.csv`, frozen build) | 38 |
| probe proteins absent from the networks | 48 |
| … of which a `ReferenceGeneProduct` in Reactome | 43 |
| … of which have any EWAS in Reactome | 17 |
| … of which have at least one IntAct partner in the networks | 39 |

Interactions touching the networks:
- 47,942 of the 78,863 interactions touch at least one of the 3,956 network proteins.
- 14,812 join two network proteins.

## Reading

- **The readout is the bottleneck, not connectivity.** A P100 probe is one phospho-site. Interaction edges can make a
  protein reachable, but they cannot create a curated phospho-form of it. Only 7 of the 96 sites exist in Reactome at
  all, whatever is added.
- **What edges could reach.** They could reach most of the 48 missing proteins (39 have a partner in the networks).
  Scoring those proteins would need a readout at the protein level ("any form of the protein moves"), not at the
  site. That is a change of protocol, decided before any value is read.
- **No direction or sign.** The edges are unsigned and undirected. They could support a reachability-style score
  (does a perturbation reach the protein), not a signed prediction. specs/112 found reachability does most of the
  work anyway, so this is not useless. But it would measure the edges, not the solver.
- **Scale.** Adding 14,812 protein–protein links among network proteins (47,942 if partners outside are imported)
  to networks of 3,956 proteins would couple most of each pathway to most of the rest. Added edges have cost the curator
  axis every time they were tried (specs/016, 019, 053, 056, 085), and over-coupling is the main error class.
