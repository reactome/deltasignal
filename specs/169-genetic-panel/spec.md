# Feature Specification: An unselected readout panel under genetic perturbation

**Feature Branch**: `specs/169-genetic-panel`

**Created**: 2026-11-21

**Status**: RESULT (2026-11-21, research.md): reached base rate 0.057, the "comparable to drug panels" reading. Registered BEFORE any expression value was read. Only the files' structure was looked at: shapes, column
names, perturbation ids and the control flag.

**Input**: specs/167 found that on every record precision is the share of reached readouts that changed plus at most
0.024. That share is 0.88–0.94 on the curator and MP-BioPath records and 0.14 on korkut. Every unselected panel
scored so far is a drug panel, so perturbation type and label source cannot be told apart. Adam approved
downloading the Replogle 2022 K562 genome-wide Perturb-seq data (2026-11-21) to separate them.

## Data

- **Source:** Replogle et al. 2022, *Cell* 185:2559. Processed figshare+ files, doi
  10.25452/figshare.plus.20029387.v1, licence CC BY 4.0.
- **Local copy:** `~/deltasignal-catalogs/external/replogle2022/`, md5 checked. It holds the pseudo-bulk files only,
  normalized and raw: 11,258 rows (514 core controls) × 8,248 measured genes.
- **Split** (`bench/datasets.tsv`): by perturbed gene symbol, fixed by a hash (`common_169.split_of`), not by data.
  - **A** (`replogle2022_a`, SELECT): scored here.
  - **B** (`replogle2022_b`, TEST): sealed, never read here. The truth reader refuses a split-B gene.

## Cases

- **Model:** the served model on the served Release 98 catalog (`20261113-1311_baf608b_r98fix`), on a dedicated
  server at this commit.
- **Perturbed genes:** for each pathway, every split-A gene whose forms `/api/pins` resolves to pins (basis root or
  root_cycle).
- **Knockdown:** pinned at 0, the knockout protocol of record. Partial knockdown is not modelled.
- **Readouts:** every measured target gene T ≠ G that is the product of a curated gene-expression-regulated event
  (Reactome 98) whose output entity is in the pathway. It is valued with `/api/readouts` forms mode.
- **Unit:** (pathway, G, T).
- **Reached:** the predicted fold is not 1 at the record's precision (`fold_rules`).

## Truth

- **Engagement:** G's row is a non-control row with knockdown `fold_expr` ≤ 0.5 and ≥ 25 cells; if several pass,
  the one with the most cells. A G with no such row is dropped.
- **Measured fold:** the raw pseudo-bulk X[row, T] divided by the mean over the core-control rows.
- **Primary label:** korkut's convention, |log2 fold| ≥ 0.585 (1.5×). This is the scale on which korkut's reached
  base rate was 0.142–0.225.
- **Secondary label:** the model's own cutoffs (< 0.85 / ≥ 1.15). These are more exposed to pseudo-bulk sampling
  noise.

## Measurement (`score_169.py`, the specs/167 reader)

Per label:
- the base rate among reached and among unreached readouts;
- the call rate;
- precision;
- lift = precision − reached base rate, with a 95% bootstrap over perturbed genes (10,000 resamples, seed 167);
- the change rate among declined readouts.

Descriptive: the model re-called at 1.5× (matched cutoff).

## How it is read (fixed now)

Primary label, base rate among reached readouts:
- **≥ 0.50:** genetic perturbations do carry through the curated routes. The gene axes' high precision is not only
  readout selection, so specs/165's reading stands in substance.
- **≤ 0.25:** comparable to the drug panels. The gene axes' precision comes from how their readouts were chosen.
- **Between:** intermediate; neither reading holds.

The lift is reported either way. An interval excluding 0 means the call adds to reachability on this panel.

## Caveats stated now

- **Readout type:** the readouts are transcripts. The model reads the curated product of the expression event.
- **Cell type:** K562 is one cell line, whose regulatory state the pathways do not encode.
- **Time:** CRISPRi readouts come days after knockdown, close to the model's steady state.
- **Pairs:** a (G, T) pair in several pathways counts once per pathway, as on the other axes.

## Requirements

- **FR-001:** this file, the scripts (`bench/analysis/169/`) and their tests (`bench/test_169.py`) are committed and
  pushed before `truth_169.py` or any read of X.
- **FR-002:** nothing is fitted. Split B stays sealed.
- **FR-003:** the result goes in `research.md` whatever it says. Raw data and per-case outputs stay outside the repo,
  under `~/deltasignal-catalogs/analysis/169/`.
