# External dataset shortlist (specs/096 section 4)

Drawn up 2026-10-23 from published descriptions only. Nothing has been downloaded or scored.
Adam approves the choice before any download (no account, terms or contact without approval).

## Criteria

- Human cells; perturbation of genes in the benchmark pathways (knockout, knockdown,
  overexpression, or a drug with a known target).
- Readouts that map to Reactome entities in those pathways: proteins or phosphoproteins (best: the
  benchmark readouts are mostly protein forms), or transcripts (only for transcriptional targets).
- Public, research use allowed. Not used before: decryptM 1 and korkut2015 are already FIT
  (specs/051, 062).

## Candidates

| # | dataset | perturbations | readouts | access | fit | concerns |
|---|---|---|---|---|---|---|
| 1 | **CPPA**, Cancer Perturbed Proteomics Atlas (Zhao et al. 2020, *Cancer Cell*; tcpaportal.org/cppa) | ~170 compounds across ~320 cell lines | RPPA, ~200+ proteins and phosphoproteins | public portal; "not required for academic use" | **high**: phospho-protein readouts in AKT/MAPK/mTOR/cell-cycle signalling, many perturbations | drugs, not genes: needs a pre-registered target map (as specs/051 R4 used kinobeads); same assay family as korkut2015 |
| 2 | **HPN-DREAM** breast cancer network inference (Hill et al. 2016, *Nat Methods*) | kinase inhibitors (and DMSO) under 8 ligand stimuli, 4 breast cancer lines, time courses | RPPA, ~45 phosphoproteins | Synapse (needs an account and acceptance of terms: Adam's call) | high for signalling (EGFR/ERBB/PIP3/RAF) | inhibitors only; time courses (a steady-state call needs a rule fixed in advance); small protein panel |
| 3 | **LINCS L1000 CRISPR knockouts** (CMap; GEO GSE92742 / GSE70138) | ~7,500 single-gene CRISPR knockouts, many cell lines | 978 measured transcripts (+ inferred) | public GEO | genetic perturbations (no target map) | mRNA readouts: only transcriptional targets in the pathways (TP53, MYC, E2F targets); cell-line context varies |
| 4 | **Genome-wide Perturb-seq** K562 CRISPRi (Replogle et al. 2022, *Cell*) | knockdown of nearly every expressed gene | single-cell transcriptomes | public (scPerturb) | genetic, very large | mRNA readouts as #3; one cell line (K562) |
| 5 | **LINCS P100** (Litichevskiy et al. 2018, *Cell Systems*) | 90 drugs × 6 cell lines | ~96 phosphopeptides | public (clue.io, PanoramaPublic); CC BY-NC-ND | phosphosite readouts | drugs only; small fixed peptide panel, many sites outside the pathways |

## Recommendation

- **Primary: CPPA (#1).** Phospho-protein readouts are what the model predicts best, the
  perturbations cover the benchmark's signalling pathways, and access needs no account. The target
  map for its compounds must be fixed and committed before scoring.
- **Secondary, genetic: LINCS L1000 CRISPR knockouts (#3)**, restricted to readouts that are
  transcriptional targets in the network, as a check that does not depend on drug targets.
- HPN-DREAM (#2) is a good fit but needs a Synapse account and terms, so it is Adam's decision.

## Sources (descriptions only)

- CPPA portal and paper: https://bioinformatics.mdanderson.org/public-software/cppa/,
  https://odin.mdacc.tmc.edu/~hliang1/publication/CancerCell_PerturbedRPPA.pdf
- HPN-DREAM: https://dreamchallenges.org/dream-8-hpn-dream-breast-cancer-network-inference-challenge/,
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4854847
- LINCS L1000: https://omicsdi.org/dataset/geo/GSE92742, https://pubmed.ncbi.nlm.nih.gov/40308216/
- Perturb-seq: https://pertpy.readthedocs.io/en/stable/usage/data/pertpy.data.replogle_2022_k562_gwps.html
- LINCS P100: https://dspace.mit.edu/handle/1721.1/125709, https://link.springer.com/10.1038/s41597-021-01008-4
