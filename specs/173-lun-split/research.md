# specs/173 research

## Format check (2026-11-26, identifiers and headers only)

Archive sha256 verified against specs/126. Extracted to `~/deltasignal-catalogs/analysis/173/work/` (outside the repo):
`complete_dat.csv` (screening), `bindat_perbin.csv`, `bindat_persample.csv`, `Data S1 description.docx`.

- **`complete_dat.csv`**: columns `experiment, marker, perturbation, row_col, timepoint`, then 31 channels (30 readouts
  and GFP). 1,450,108 rows; 3 experiments (KP01-02, KP03, KS01-02); 115 `marker` values (about 105 proteins plus
  FLAG-GFP and untransfected controls); `perturbation` is always `EGF`; `timepoint` 0 or 10. **The file is truncated in
  the deposit**: its size equals the archive entry (655,181,824 bytes), it ends in a comma, and its last row has 24 of
  36 fields.
- **`bindat_perbin.csv`**: `experiment, marker, perturbation, row_col, timepoint, origin, target, bin, counts, fit_var,
  median, overall_var, var_ratio`. 423,630 rows; `origin` is always GFP; 30 targets; bins 0–9; 14 plate pairs
  (KP01-02, KP03, KS01-02 … KS19-20, PS01-02, PS03-04), each at timepoints 0 and 10. 702 `marker` values: all 649 Table
  S1 markers (636 genes; split 292 constructs A / 357 B) plus FLAG-GFP1–27 (no 3) and untransfected1–27.
- **`bindat_persample.csv`**: per-sample statistics, including the authors' significance filters (`bin_dat_sigfil…`).
  Header and the identifier columns only were read. It is not used.
- **`Data S1 description.docx`**: Data S1 is result plots (biaxial and bar plots). Not opened.
- **Paper (Europe PMC full text)**: only the data-availability statement (all data at the one Mendeley record) and the
  methods paragraphs on binning (ten equal-width bins over the 2.5–97.5% range of each sample) were printed.

This led to amendment 1.

## Split and sealing (2026-11-26)

`lun_a_173.py split` on `bindat_perbin.csv`: 177,480 rows in A, 214,440 in B, 31,710 control rows in both, none
unmapped. Files under `~/deltasignal-catalogs/external/lun2019/split173/`:
- `lun2019_a_perbin.csv`: sha256 `66511fd8b4854bcfcd12bcbd292072001ada65d490e51d247bf3b7f06707c5ef`, read-only;
- `lun2019_b_perbin.csv`: sha256 `6ae9507d228df33f7a79d9b34be2368971201bbf2c6ecbb38303f150e46a1a17`, **mode 000**.

The per-cell extract, the per-sample statistics file and the Data S1 description were deleted. The archive is back at
mode 000.

## Half A, scored (code at `bad27ef0`; output `~/deltasignal-catalogs/analysis/173/score_a.json`)

Frozen predictions: specs/127's values-blind top114 file (sha256 `0c0ec033…`, unchanged).

**Positive control: passes.** 10 of 11 measured self pairs read UP (RPS6KB1 at p-p70S6K reads NOT).

**Measured:** 8,640 (gene, readout) pairs on half A. P ∩ A has 3,214 pairs, of which 3,127 are measured: 2,383 NOT,
538 UP, 206 DOWN.

**R1 (three-class balanced accuracy, protein bootstrap):**

| | balanced accuracy |
|---|---|
| frozen model | 0.379 |
| reachability | 0.372 |
| always-NOT | 0.333 |

- frozen − reachability: +0.007 [−0.008, +0.024], **unresolved**;
- frozen − always-NOT: +0.046 [+0.023, +0.070], **gain**.
- Kinases (2,805 pairs): frozen 0.385, reachability 0.375. Phosphatases (322): frozen 0.308, reachability 0.330.

**R2 (does a route matter?):**
- Reached pairs (792) change in 30.8%, unreached (2,335) in 21.4%. **Enrichment 1.44 [1.13, 1.82].**
- Frozen model change calls on reached pairs: 548. Precision 0.334 (any change), 0.261 (right direction).
  **Lift +0.026 [−0.000, +0.054].**
- Measured direction among reached pairs: 186 UP, 58 DOWN, 548 NOT.

**Reading against the registration.** The registered reading named two outcomes: an enrichment interval above 1.5
(fit for development) and one including 1 (not fit). The result is neither: the interval excludes 1, so routes carry
signal here (K562: 1.0), but its lower end is 1.13, short of the 1.5 bar. The registration did not cover this middle
case. It is recorded as **weaker route signal than registered, but present**. The lift interval just includes 0
(unresolved).

## Traced by readout (half A, P ∩ A)

- **Background change.** STAT5, ERK1/2, FAK, p38, STAT3, AKT, p70S6K and SRC change for 35–46% of over-expressed
  proteins, reached or not. That is why the unreached base rate is 21% (K562: 6%; korkut: 5%).
- **Route-specific readouts.** p-PLCg2 changes in 82% of reached pairs against 13% unreached: BTK, TEC, ITK, SYK, LCK,
  SRC and JAK3 raise it, and the model calls 7 of 9 right. p-p38 (52% / 30%), p-STAT3 (56% / 31%), p-MKK3/6
  (38% / 13%), p-MKK4 (29% / 13%) and p-MAPKAPK2 (21% / 5%) follow routes too. p-p38 has 23 right UP calls (SRC family,
  receptor tyrosine kinases, MAP2K6).
- **Where the model fails.**
  - **p-AKT: 4 of 48 calls right.** The model reads UP under PI3K subunits, receptor tyrosine kinases and SRC.
    Measured, 24 do not move and 16 fall, including EGFR, SRC, ALK, BTK, PIK3CG and the regulatory subunit PIK3R3.
    The network has no way to show an over-expressed regulatory subunit damping PI3K.
  - **Total β-catenin: 0 of 44.** It is called UP under 44 proteins and changes for 7% of all proteins.
  - Also GSK3β (3 of 40), CREB/ATF1 (5 of 40), MAPKAPK2 (9 of 45). p53, SMAD2/3 and c-RAF change for under 10% of
    proteins, and every model call there is wrong.

## What this means

- Lun 2019 half A is the first genetic panel where routes carry signal (enrichment 1.44, against 1.0 for K562). The
  signal is concentrated in specific readouts (PLCγ2, p38, STAT3, the MKKs).
- The frozen model is no better than reachability. Its failures are specific and traceable: AKT under PI3K and
  receptor over-expression, and total β-catenin. That makes half A usable for development: each failure is a
  candidate rule or curation question, to be judged on half B.
- Half B stays sealed.
