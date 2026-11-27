# specs/183: why does p-AKT fall under kinase over-expression?

**Status:** registered 2026-11-27, before the questions below were measured.
**Requested by:** Adam, 2026-11-27 ("find solutions"); the open lead from specs/181.

## Why

p-AKT S473 is the readout that falls most often in Lun 2019 (106 of 435 measured falls on Release 98's scored pairs),
and 36 of the model's 87 wrong-direction calls are p-AKT called UP and measured DOWN. Nothing curated explains it.
Both halves are FIT; this is a trace, and nothing is adopted from it.

## Questions, in order

1. **Who lowers it?** The proteins behind the falls by library (kinase or phosphatase), tyrosine-kinase class
   (specs/177), and how many other readouts each changes (specs/176's protein activity). Is the fall a property of
   active proteins in general (share of proteins lowering p-AKT by activity band)?
2. **What falls with it?** For proteins that lower p-AKT against those that do not: the rate of UP and DOWN at every
   other readout. The registered expectation distinguishes three explanations:
   - *feedback from a raised growth pathway* (ERK or mTORC1 inhibiting upstream of AKT): p-ERK1/2 or p-S6 / p-p70S6K
     rise in the same proteins;
   - *the PI3K branch as a whole falling*: p-GSK3β, p-S6, p-p70S6K, p-4EBP1 and p-PDK1 fall together with p-AKT;
   - *a general cost of over-expression*: the fall does not depend on what the protein is, and scales with nothing but
     activity.
3. **Is it real?** Repeat rate of the p-AKT falls in the EGF well (specs/176 T1a's check), and their size against δ.
4. **Does Reactome hold a cause?** For the proteins with a p-AKT fall and a curated route to p-AKT: the sign of the
   shortest route, and whether any negative route exists (the tracer of specs/174 on the Release 98 build).

## Reading

Each explanation is recorded as supported or not by question 2's table. If one is supported and Reactome holds the
edges (question 4), the specific loop is named as a candidate for a later, separately registered arm. If none is, the
fall is recorded as outside the curated network.
