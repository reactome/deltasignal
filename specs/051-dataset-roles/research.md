# specs/051 research

## R1. Survey of public perturbation datasets (2026-10-12, web check by a subagent)

The criteria: cell lines only, open download, a licence recorded, and
perturbations of genes or of drugs with known targets, with protein or
phosphosite readouts. "Verified" means the landing page was opened.

| dataset | accession | licence / access | perturbations | readout | proposed role |
|---|---|---|---|---|---|
| LINCS PCCSE P100 (Litichevskiy 2018; Dele-Oni 2021) | Panorama "LINCS PCCSE P100 2020" (verified), LDS-41234 | CC-BY 4.0, open | 119 compounds, 7 cell lines, mostly 1 µM, 3 h | targeted MS, 96 phosphopeptides | **FIT** |
| CPPA (Zhao 2020 Cancer Cell) | tcpaportal.org/cppa | **no licence stated**, open | ~170 drugs, 319 cell lines, ~15k samples | RPPA, ~210 proteins / phospho | FIT, if MDACC confirms the licence |
| Hijazi 2020 Nat Biotech | PRIDE PXD015943 (verified) | CC0, open | kinase-inhibitor panel (~60, unverified) | phospho-MS, >6,000 sites | **TEST**, cell-line samples only; the set also has primary AML and breast tumour samples, which are excluded |
| Korkut 2015 eLife (CellBox data) | github sanderlab/CellBox `data/` (verified) | eLife CC-BY; repo MIT | 12 drugs and pairs, 89 conditions, SK-Mel-133 | RPPA, 98 readouts | **TEST**, registered below |
| HPN-DREAM8 (Hill 2016) | Synapse syn1720047 (verified) | Synapse login and challenge terms | 4 inhibitors × 8 stimuli × 4 lines, 7 time points | RPPA, ~45 phospho | TEST, if the terms allow |
| MCF10A LINCS (Gross 2022) | Synapse LINCS_MCF10A | Synapse account | 6 ligands, time course | RPPA ~300 | TEST, ligand direction only |
| Mitchell 2023 Nat Biotech | MassIVE MSV000089282 | public | 875 compounds, HCT116 | total proteome only | low: abundance, not signalling |
| Francavilla 2016; Bekker-Jensen 2020 | unverified | PRIDE | EGF / TGF-α time courses | phospho-MS | minor TEST |
| Klaeger 2017 kinobeads | PRIDE PXD005336 (verified) | EBI terms | 243 inhibitors | target affinity | **mapping**: drug → targets and potency |
| HMS LINCS KINOMEscan; Davis 2011 | LDS-1045 (verified) | HMS LINCS terms (acknowledge) | binding panel | Kd | mapping |

**Excluded:**
- CPTAC-derived sets (Mertins, Mundt) and Casado 2018 (PXD005978): patient samples.
- PhosphoSitePlus (CC BY-NC-SA, licence agreement) and SIGNOR (CC BY-NC):
  priors at most, not data.
- Perturb-seq and the CITE-seq variants: RNA, or a few surface proteins.

**Why this split.**
- The FIT sets are large, from many cell lines, and close to decryptM in
  kind: drugs, single dose, phospho.
- The TEST sets come from other labs and other technologies: RPPA, and a
  different MS pipeline.
- **The drug → target map should come from Klaeger (measured affinity), not
  from our hand-written primary-target table.** That table (in
  `decryptm_coverage.py`) is a judgement call, and a measured map makes it
  auditable.

## R2. TEST registered: Korkut 2015 (2026-10-12)

- Files: CellBox commit `6bed6b94`, `data/` (expr.csv, pert.csv,
  node_Index.csv, expr_index.txt, loo_label.csv). Kept in
  `~/deltasignal-catalogs/external/korkut2015/` with `SHA256SUMS`.
- **Read so far: the design only.**
  - The 98 readout names, e.g. MAPKpT202, AKTpS473, 4EBP1pS65, SRCpY527,
    STAT3pY705, p53, CyclinD1.
  - The 89 condition labels, as drug code|dose and pairs.
- `expr.csv` and `pert.csv` values have not been opened. They stay unread
  until a scoring protocol for this set is committed.
- Note for that protocol:
  - The readouts are antibodies, so each needs a Reactome entity mapping
    (protein, or protein plus site), made from the names alone.
  - The conditions include drug pairs, a test of non-additivity, which
    decryptM could not give.
  - It is one cell line, so the per-site results are not independent.

## Open decisions (Adam)

1. HPN-DREAM8 and MCF10A need a Synapse account and acceptance of their
   terms. That is a person's agreement, not something to accept on the
   project's behalf.
2. CPPA has no stated licence. A note to the TCPA team at MDACC would
   settle it. It is the largest FIT candidate by far.
