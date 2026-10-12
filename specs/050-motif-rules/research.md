# specs/050 research

## R1. Substrate accumulation and competition (2026-10-12): census, nothing built

Started from the open item in specs/049: an enzyme or binder knockout should
raise the species its reaction consumes (SMAD7 → TGFBR1, PTPN11 → SPRY2).

**What is curated (Neo4j, Release97).** SMAD7 does not act through
dephosphorylation in a cycle. "I-SMAD competes with SMAD2/3 for type I receptor
(TGFBR1)" (R-HSA-173512) consumes the active receptor
TGFB1:TGFBR2:p-TGFBR1 (170841) with SMAD7. "PP1 dephosphorylates TGFBR1"
(R-HSA-178178) releases a receptor (2167873) that no reaction uses, so no
cycle closes and no pool forms (TGF-β has no `pools.csv` rows). In our
network, 170841 is an `input` to five reactions and has no depletion in-edge.
The curated mechanism is **competition for a shared species**.

**Why not accumulation.** Two earlier attempts added consumption as edges
from co-inputs. A naive version lost 14pp; general consumption cut
experimental macro-F1 from 0.66 to 0.53. Both over-coupled. On the canonical
build, 23,054 entity nodes have exactly one consuming reaction. Any
supply/demand rule on a node's own level sends such a node to the cap when
its single partner is knocked out.

**The census** (`bench/analysis/competition_census.py`; build
`20261012-1401_01206ca`, results `1f58a1b`):
- **Motif:** a species S is an input to ≥ 2 distinct reactions, and is not a
  cofactor, carrier free form or pool node. The gene is a partner (co-input
  or catalyst) of one consumer Ri. The readout is reachable from a sibling
  consumer but not from Ri.
- **Prediction:** competition says the readout moves opposite to the gene.
- **The count is an upper bound.** It assumes the effect crosses the call
  threshold (1.15×).

| motif restriction | curator: wrong, would fix | curator: right, would break | experimental fix / break |
|---|---|---|---|
| any shared species | 186 | 2,275 | 20 / 27 |
| gene's branch is an inert sink (nothing downstream acts) | 170 | 716 | 5 / 18 |
| inert sink, gene has no other route to the readout | 156 | 504 (all expected no change) | — |
| Ri annotated GO "negative regulation of …" (reaction or parent pathway) | 39 | 30 (26 WNT) | 0 / 0 |
| GO negative regulation + inert sink | **29** | **4** | 0 / 0 |

The last two rows count a closed carrier free form as its gene, as
`closed_twin` makes it. Without that, the SMAD7 cases in TGF-β are missed.

**Reading.**
- As a structural rule, competition is refuted by the census. Curators
  applied sequestration logic only to known antagonists: SMAD7, PMEPA1,
  follistatin, BCL2, the IL13RA2 decoy and PTCH. Topology cannot tell those
  from ordinary binding partners.
- Curation can tell them apart, partly. GO "negative regulation of …" marks
  323 human reactions (74 on the reaction itself), including the I-SMAD
  competitions and BCL2 sequestering tBID. It misses FST, IL13RA2 and PTCH.
- With that restriction the motif is clean (29 / 4) but small, all curator,
  and concentrated in two pathways (TGF-β 17, BMP 12) and two genes.

**Decision: not built.** It needs a generator table of GO-annotated
reactions and a solver term that reads sibling consumers' partners. Those
are dependencies that do not follow network edges and can close new cycles.
That is too much for ≤ 29 curator cases in two pathways, with nothing on the
experimental axis. Kept as a measured candidate. It is worth revisiting if
the motif census (step 2 of the plan) finds GO-annotated competition in
pathways with experimental data.

**PTPN11 → SPRY2** is the other shape: the substrate's only consumer is the
knocked-out enzyme's reaction (n = 1). That is the accumulation regime
refuted above, so it stays open.
