# specs/113: does Reactome's interaction layer predict drug responses? (LINCS P100, scored once)

**Status:** pre-registered 2026-10-27, before any P100 value was read.

Adam, 2026-10-27: approved (on research.md's census) testing his idea, extra connections from the database's
interaction records, as a reachability layer scored once on LINCS P100. The frozen model (specs/111) is not changed and
not scored: only 4 of P100's 96 sites occur in its networks. This tests the wiring, with and without interactions,
the part that did most of the work on HPN-DREAM (specs/112). Scoring spends P100: it becomes FIT afterwards.

## Data

LINCS P100 (Litichevskiy et al. 2018, GEO GSE101406, CC BY 4.0), already downloaded and checksummed (specs/111):
89 compounds × 6 cell lines, 3 h, 96 phosphopeptide probes. Level 4 (`..._Level4_DIFF_n1684x96.gctx`) is the value
file. Only the analyte, instance and perturbagen tables have been read. `cell_info` is not read beyond `cell_id`.

## The drug → target map (fixed before any value is read)

Each compound is identified in ChEMBL 37 (the release specs/051 R4 pinned) by its standard InChIKey from
`pert_info`. If there is no exact match, the first InChIKey block (connectivity) is used, then the preferred name or a
synonym. Unmatched compounds are listed. Targets come from the first of these sources that gives any:

1. **Kinobeads** (specs/051 R4, unchanged): the molecule is in the Klaeger deposit, and the proteins have
   Kd_app ≤ 1,000 nM. The map is also built at 100 nM and 10 µM and reported, but not chosen on outcome.
2. **ChEMBL 37 `drug_mechanism`**: human targets of the molecule's curated mechanisms. A protein-complex or
   protein-family target is expanded to its listed components (curated mechanisms name families such as HDACs).
3. **ChEMBL 37 binding activities**: human single-protein targets with pChEMBL ≥ 6 (≤ 1 µM, R4's cutoff).
4. Otherwise the compound has no target and is not scored. It is listed.

Proteins are UniProt accessions without isoform suffix. The map and each compound's source are committed before any
value is read.

## Scores (computed and committed before any value is read)

A probe is scored by its protein (UniProt), not its site. Interactions attach to proteins, so this is the only level
at which they can act (research.md). A protein is *in the networks* if a node of the frozen catalog
`20261026-1349_8b6527e` holds it: a `member_leaves` leaf whose `leaf_refs.csv` reference is that accession.

- **A, curated reach.** 1 if, in any one of the 92 networks, a node holding a target has a directed path (any edge of
  `logic_network.csv`, length ≥ 0) to a node holding the probe's protein, else 0. Networks are not stitched together.
- **B, curated reach plus one interaction at each end.** 2 if A is 1. Otherwise 1 if curated reach succeeds once the
  sources are widened to the targets plus their IntAct partners in Reactome 97, and the probe is widened to its own
  protein plus its IntAct partners. Otherwise 0. A partner is any other `ReferenceGeneProduct` sharing an
  `UndirectedInteraction` (Reactome 97, 78,863 interactions).
- **C, interactions alone (secondary).** Distance from a target to the probe's protein in the IntAct graph: 0 or 1
  hops gives 2, 2 hops gives 1, more gives 0.
- **Constant:** AUROC 0.5.

## Truth and scoring

- **Context:** a compound in a cell line, at the highest dose given there.
- **Value:** the median of that context's replicates, per probe.
- **Changed:** |value − median of the cell line's DMSO wells| ≥ 3 × 1.4826 × MAD of those DMSO wells, per probe and
  cell line. A probe is dropped in a cell line if it has fewer than three DMSO values or zero MAD. NaN values are
  dropped.
- **Score:** AUROC per context, ties counted half. A context is scored if it has at least one changed and one
  unchanged probe.
- **Headline:** the mean AUROC over scored contexts, for each arm.
- **Interval:** a bootstrap over compounds (cluster resampling; contexts of one compound share targets): 10,000
  resamples, numpy seed 113, 2.5–97.5 percentiles of each arm's mean AUROC and of B − A.
- **Random-interaction control:** IntAct partner sets are relabelled by permuting protein identities among proteins
  that have interactions, which keeps each protein's degree. Then B is recomputed. This is done 200 times, seed 113,
  and B's percentile among the draws is reported.

## Predictions

- **H1.** A's mean AUROC is above 0.5, and its interval excludes 0.5. Predicted value at most 0.56: most compounds
  (HDAC, BET, EZH2 and other chromatin drugs) have targets outside the signalling networks.
- **H2 (Adam's idea).** B − A > 0 with an interval excluding 0. **I predict it fails** (interval includes 0): an
  unsigned hop at each end reaches most proteins, which leaves little to rank.
- **H3.** B is above the 95th percentile of the random-interaction control. I predict it fails.
- **H4 (secondary).** C ≥ B.

The result is reported whatever it shows. Nothing is tuned after scoring: no change to cutoffs, hop counts, dose or
threshold.
