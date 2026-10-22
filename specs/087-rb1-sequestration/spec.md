# specs/087: RB1 represses E2F only by sequestration, which Reactome does not encode (analysis)

**Status:** analysis, 2026-10-22. No rule proposed.

Adam, 2026-10-21: "ok" to tracing RB1 knockout in Mitotic G1 after specs/086.

## The cases

In Mitotic G1 (R-HSA-453279) on canonical `20261021-2218_f99c4ab`, every RB1 knockout and
overexpression case reads exactly 1.000:
- **Curator:** E2F targets (CCNA1, CCNE1, CDC6, CDT1, DHFR, TK1, …) are expected UP under
  knockout and DOWN under overexpression; 16 and 13 cases are wrong.
- **Experimental:** 4 knockout cases expected UP, all wrong.

## The trace (`rb1.jl`)

1. RB1 [cytosol] (R-HSA-9660621), the pinned starting copy, translocates to the nucleus.
2. Nuclear RB1 is an input of "RB1 binds and inhibits E2F1/2/3:DP complexes"
   (R-HSA-9018017, 360 variant copies), whose product is RB1:E2F1,E2F2,E2F3:TFDP1,TFDP2
   (R-HSA-68644).
3. **RB1:E2F leads nowhere.** In our network it has only dissociation sinks. In Release 97,
   no reaction in any human pathway consumes, is regulated by or is catalysed by a complex
   of RB1 and E2F (Neo4j, checked).
4. The E2F target genes are driven by E2F1:TFDP1,TFDP2 and its gene complexes. Free E2F:DP
   reads 1.0 whatever RB1 does: binding to RB1 uses it as an input, and the model does not
   lower a reaction's inputs.

So RB1's repression of E2F is curated only as sequestration: RB1 takes E2F:DP out of
circulation by binding it. The model cannot represent that except through depletion edges.

## Why no rule

specs/063 measured exactly this class, a dead-end binding complex depleting its free
partner, in three forms, and all three were refuted (held-out −1,180, −1,686, −1,121). No structural
signature separates a sequestering complex (RB1:E2F, follistatin, USP18) from a terminal
product (RUNX1's co-activator complexes). The reaction's name here says "binds and
inhibits", but names are not data.

**What would fix it is curation, not modelling:** a curated statement that RB1:E2F
negatively regulates the E2F target-gene reactions, or a curated flag marking a binding
reaction as sequestering. This is a suggestion for Adam as a Reactome curator; no one
outside the project has been contacted.

The same limit is behind part of specs/063's 348 held-out misses.
