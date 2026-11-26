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
