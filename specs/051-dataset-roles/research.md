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

## R3. Rewired-network null: the curated wiring carries signal (2026-10-13)

Registered in spec.md before any build (495e9e9). Builds and arms:

- Source build `20261012-1401_01206ca`. Null builds `..._rewire1` to
  `..._rewire5`, made by `bench/rewire_network.py` at 0159108.
  - 257,125 reaction-level edges per build. Between 247,364 and 247,450
    moved (96%), from about 2.5 M accepted swaps.
- All six arms ran at 0159108 with `DS_CYCLE_MODE=off`,
  `DS_CARRIER_LOOP_MODE=off` and `DS_PIN_NETWORK_ROOT` set to the source
  build, so the same nodes were pinned in every arm.
  - Control `wire0`: `results/0159108/wire0` on the source build.
  - Seeds: `results/0159108/rw<seed>` on each null build.
- Scored by `bench/analysis/rewire_null.py`, with the JSON at
  `results/0159108/rewire_null.json` on the source build.
- Every case was paired across the six arms: 24,100 curator and 849
  experimental cases, none lost.

| axis | n | genes | wire0 acc / mF1 | seeds acc (mean) | seeds mF1 | wire0 − seeds, correct, 95% by gene | rule |
|---|---|---|---|---|---|---|---|
| curator tuning (FIT) | 5,100 | 94 | 0.761 / 0.749 | 0.434 | 0.406–0.427 | [+0.264, +0.390] | **pass** |
| curator held-out (SELECT) | 19,000 | 463 | 0.876 / 0.838 | 0.588 | 0.483–0.517 | [+0.264, +0.311] | **pass** |
| experimental (FIT) | 849 | 81 | 0.702 / 0.621 | 0.417 | 0.322–0.423 | [+0.174, +0.375] | **pass** |

Canonical macro-F1 on the same cases is 0.765, 0.853 and 0.616. Turning
off the pool and carrier rules costs about 0.015 on the curator axes, and
nothing on the experimental axis.

**Verdict, per the registered rule:** on every axis, `wire0` beats all five
seeds on macro-F1, and the gene-clustered interval excludes 0. That
includes the experimental axis, which was the informative one.

**What the null is and is not.**
- The rewired networks do not fall back to "no change". They over-call UP
  (38% of tuning predictions against 27% expected), because rewiring
  scatters inhibitors and depletion edges.
- On the experimental axis they score below always-majority (accuracy 0.417
  against 0.468; majority macro-F1 0.212). So this null is a floor, and
  passing it shows the model is not running on degree and pinning alone.
  - It does not show the model beats simpler uses of the same wiring.
  - That is the reachability baseline's job, and on decryptM reachability
    won (specs/021).
- Rewired solves converge less often (curator 1,511–1,631 of 1,725 against
  1,646; experimental 175–214 of 244 against 218). The case dump does not
  record convergence per case, so the gap is not split by it.

**Traced: cases every seed got right and `wire0` missed.**
- Curator tuning 146, held-out 194, experimental 28. The reverse (`wire0`
  right, every seed wrong) is 905, 2,304 and 150.
- 93 of the tuning 146 are in Mitotic G1, and 85 of those 93 are
  `no_path`.
  - They are overexpressions (ABL1, SRC, PTK6, JAK2, CCND1, CDK4) whose
    readout the real network does not reach. The rewired network connects
    them at random, and its UP bias lands.
  - These are the severed routes already traced in specs/047, not evidence
    for the null.

## R4. Drug → target map from measured affinity: registration (2026-10-13, before any value is read)

Protocol point 7 retires the hand-written `DRUG_TARGETS` table
(`bench/analysis/decryptm_coverage.py`). It is replaced by a map built from
measured affinity.

**Source.**
- Klaeger 2017 is not open access. Science returns 403 to a scripted
  download, and Europe PMC reports the manuscript (PMC6542668) as not open
  access. The supplementary tables are therefore not used.
- The same measurements were deposited by the Kuster lab in **ChEMBL**
  (source 48, `TUM_PROTEOMIC_KUSTER`; document `CHEMBL3991601`). That is
  70,505 activities, all `Kd apparent` in nM, under ChEMBL's
  **CC BY-SA 3.0** licence.
  - Pinned to **ChEMBL 37** (released 2026-05-01).
  - The fetched records and their SHA256 are kept under
    `~/deltasignal-catalogs/external/klaeger2017_chembl37/`.
- What was read before this registration: the record structure only.
  - Measured binding is `relation '='` with comment Active.
  - No binding is `relation '>'` with comment Inactive.
  - No drug–target value has been read.

**Rule (primary).**
- Protein P is a target of drug D if the deposited set has a record for
  (D, P) with relation `=` and **Kd_app ≤ 1,000 nM**.
- Only single-protein ChEMBL targets count, mapped to UniProt, then to the
  HGNC symbol. Complex and family targets are skipped and counted.
- Non-kinase proteins that kinobeads co-capture are kept as recorded, and
  their number is reported. They are not filtered by hand.

**Drug names.**
- Each decryptM drug name maps to one ChEMBL molecule by preferred name or
  synonym. The mapping is written out and committed before targets are
  joined.

**Drugs the assay cannot see.**
- Kinobeads measure ATP-site competition in lysate. Trastuzumab and
  pertuzumab (antibodies) and SHP099 (an allosteric phosphatase inhibitor)
  cannot appear.
- Any drug with no record, or no target at the cutoff, keeps its hand
  target. Its `source` column says `hand`, and it is listed. Allosteric MEK
  inhibitors and tideglusib may land here.
- They are not dropped, because the hand table is the better definition for
  a drug the assay cannot measure.

**Sensitivity.**
- The map is also built at 100 nM and 10,000 nM. Both are reported, and
  neither is chosen on outcome.
- The cutoff is not changed after a score.

**Use.**
- The map is a data definition for every future drug dataset (TEST sets
  included). It is fixed by this registration, not by any score.
- decryptM 1 is FIT, so it may be re-scored under the new map. The number
  is reported as a FIT-set measurement, not as accuracy.

## Open decisions (Adam)

1. HPN-DREAM8 and MCF10A need a Synapse account and acceptance of their
   terms. That is a person's agreement, not something to accept on the
   project's behalf.
2. CPPA has no stated licence. A note to the TCPA team at MDACC would
   settle it. It is the largest FIT candidate by far.
