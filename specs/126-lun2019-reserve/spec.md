# specs/126: Lun et al. 2019 as the reserve test

**Status:** registered 2026-11-02. Adam, 2026-11-02, on the dataset survey: "I don't see any reason why Lun 2019
shouldn't be used as the reserve. stick to best practices. make sure to be as thorough as possible. no rushing to
conclusions." `bench/datasets.tsv` lists it as TEST from today. Nothing has been downloaded, and no value has been read.

**Why a reserve.** Every dataset we hold has shaped the rules, or will shape a learned model (specs/051, 125). The CPPA
test portions are the one clean test left, and one scoring spends them. A reserve chosen now, before any learning, by
its design and coverage alone, keeps a second independent test. Lun 2019 also covers the area our data do not:
genetic perturbation. Everything we hold is drug data.

## The dataset (from the paper's methods and the data record; no results read)

- **Lun X, et al.** Analysis of the human kinome and phosphatome by mass cytometry reveals overexpression-induced
  effects on cancer-related signaling. *Mol Cell* 74(5):1086–1102 (2019), doi:10.1016/j.molcel.2019.04.021, PMC6561723.
  The article is CC BY-NC-ND 4.0. Data: Mendeley Data doi:10.17632/3kh7ypz232.1, version 1 (published 2019-05-16), CC BY 4.0.
- **Perturbations.** 649 proteins of interest (541 kinases from Addgene kit 1000000014; 108 phosphatases from Dharmacon
  OHS4941), each as an N-terminal GFP fusion (pcDNA5 FRT TO), transiently transfected into HEK293T. Measured 18 h
  after transfection.
- **Conditions.** Basal, and EGF 100 ng/mL for 10 min. A 1 h EGF time course for 39 proteins. Nine kinase–phosphatase
  pairs. A375 (BRAF V600E) with vemurafenib 1 µM and/or CI1040 10 µM, for 3 h or 2 days.
- **Controls.** 108 control samples (FLAG-GFP and untransfected), kinase-dead mutants (AKT3 K177M, AXL K567R, MAPK3
  K71R, PRKCE K437W, MAP2K1 K97M), and MEK1DD in A375.
- **Readout.** CyTOF, 36 markers: ~30 phosphosites, 5 non-signalling markers, and anti-GFP. Anti-GFP measures the
  over-expressed protein in each cell, so every sample is a dose series across cells.

## What has been seen, and what must not be

- **Seen.** The methods and Key Resources Table, through a reader told to exclude results; the data record's metadata
  (title, version, licence); and the antibody targets below. A survey agent read public abstracts and landing pages,
  and reports recording no results.
- **Not to be read before scoring.** Any data file. Supplementary Tables S3–S7, which are results (shape clusters, the
  signed relationship matrix, non-catalytic kinases, signalling relationships, ERK-activating proteins). The paper's
  results text and figures.
- **Design files** may be read. Table S1 lists the proteins of interest; Table S2 is the antibody panel.

## Coverage (graph only, `bench/analysis/coverage_114.py`; panel in `~/deltasignal-catalogs/analysis/126/`)

32 of the 36 markers are scorable readouts. Four are excluded: GAPDH, cleaved PARP, p-HH3, and GFP (GFP is the dose).
Antibody coordinates are translated to human numbering as in specs/062, from memory; they still need checking
against UniProt.

| build | protein present | site curated in Reactome 97 | site present in network |
|---|---|---|---|
| frozen 92 pathways (`20261026-1349_8b6527e`) | 29 / 32 | 28 / 30 | 23 / 30 |
| whole Reactome, 28 top-level (`20261027-1346_8b6527e_top114`) | 32 / 32 | 28 / 30 | 28 / 30 |

Two sites are not curated in Release 97: RELA S529 and MTOR S2448. Two readouts are total-protein only (cyclin B1 and
non-phospho β-catenin). The frozen build lacks p70S6K, MARCKS and 4E-BP1.

**Perturbation coverage is not yet known.** It needs Table S1, the list of 649 proteins. Table S1 is a design file but
a download, so it waits for Adam's approval.

## Open before anything else

1. **Do the data exist?** The Mendeley record's metadata reports a size of 0 and an empty file list through the public
   API. The web page shows only a "Download All" button. This may be an API quirk, or the files may be hosted
   elsewhere. Confirming it needs a download.
2. **Acquisition and sealing**, proposed and awaiting Adam's approval:
   - download the data record and Table S1 once, to `~/deltasignal-catalogs/external/lun2019/` (outside the repo);
   - record each file's sha256 and size in this spec;
   - make the data files unreadable (`chmod 000`) until the scoring spec is registered;
   - read only Table S1 and the file names.
3. **HEK293T resting levels.** HEK293T may not be in CPPA or DepMap. Its resting levels will come from class defaults
   or another source chosen before unsealing.

## Scoring: to be registered before the data are opened, not now

These are questions only. The scoring will be decided after the learned-model plan (specs/125) settles, and before any
file is opened:

- **Condition.** Basal, 18 h after transfection, is close to a steady state and fits the steady-state model. EGF 10 min
  needs the time model.
- **What is predicted.** For each (protein over-expressed, readout), the direction of the readout's change with
  increasing dose of that protein. This is computed by us from the cell-level data with a method fixed in advance, not
  taken from the authors' tables. Possibly also the shape.
- **Comparators.** Reachability, always-no-change, and the frozen model alongside any learned model.
- **What it can and cannot say.** It tests over-expression of a GFP fusion in one non-cancer cell line, at transient
  and supraphysiological doses. It cannot test knockouts or drugs (except the A375 subset). The tag may change function.
- **Multiplicity.** One primary measure. Everything else is reported as secondary.
