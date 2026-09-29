# specs/041: validate against established biology, with the route classified

**Why (Adam, 2026-09-29).** Both benchmarks are frozen 2019 answer sets.
- The curator axis is the curators' *expectations*, and MP-BioPath's hand
  edits were made to match those same expectations (specs/034 §12d). So
  matching them is partly circular.
- The goal is a model useful for analysing biological data, so it should be
  tested against what is established in biology, and scope must be accounted
  for.
- An expected outcome can be real, yet arise through a mechanism outside the
  pathway (Adam).

## Route classes (applied before any scoring)

For each (perturbed gene, readout, expected direction), the class comes from
the Reactome graph (Neo4j, release 97):

| class | meaning | counts as |
|---|---|---|
| **a. in-pathway** | a route from the gene to the readout exists within the pathway's own events | a test of the model |
| **b. cross-pathway** | routes exist only through other Reactome pathways (named) | a test of cross-pathway stitching, not of the single-pathway model |
| **c. not in Reactome** | no route anywhere in Reactome | a curation gap (a list for Reactome), not a model error |

- A route is a directed chain of reactions (input/output, catalyst, regulator)
  connecting an entity containing the gene to an entity containing the readout.
- The search is bounded in depth and records the reactions used.
- The exact definition is fixed by the script before any scores are read.

## Step 1: classify the existing experimental axis (849 cases)

- Route class per case.
- Our canonical prediction per case (build `20260928-1110_06ccb63`, solver
  `8d6d27c`): right or wrong, by class.
- Result: how much of the 5.9pp experimental gap to MP-BioPath is reachable
  from the pathways as curated.

## Step 2: pilot on RAF/MAPK and PI3K/AKT

- Expected outcomes are drafted blind twice (Opus, Fable), **without reading
  the curator or experimental ground truth**. For each key gene (KO and
  overexpression), and each main readout:
  - the direction;
  - its basis (textbook consensus, or a named result; citations to be
    verified);
  - a confidence level.
- The drafts are reconciled, then **reviewed by Adam**.
- Each outcome is route-classified, then the model is run on exactly those
  perturbations.
- The output is a per-pathway table: expected, class, predicted, match. Every
  class (a) mismatch is traced (network, maths or curation).
- Only then is it compared with the curator answers, to see where curators and
  established biology disagree.

## Use

- A class (a) mismatch pattern that recurs across pathways becomes a candidate
  rule. It is judged on other pathways and on the experimental axis, never on
  the pathway it came from.
- Class (c) outcomes go to Adam as curation gaps.
