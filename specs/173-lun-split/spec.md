# specs/173: Lun 2019 split into a development half and a sealed test half

**Status:** registered 2026-11-26, before the archive is unsealed. No value of the dataset has been read.
**Decision (Adam, 2026-11-26):** add a genetic-perturbation dataset with protein-state readouts beside the drug panels,
using Lun 2019 split by over-expressed protein, half for development and half kept sealed ("okay go ahead"). It amends
specs/126 (Lun as the whole reserve) and specs/127 (scored once, at specs/125 stage 5). specs/127 left the timing of a
frozen-only scoring to Adam; this is that decision, for half A only.

## Why

- Every panel we develop on is a drug panel. On korkut, 94% of false alarms are in drug sets holding a drug that did
  not engage its target (specs/138, 161). Over-expression has no engagement question.
- The one genetic panel we scored, Replogle 2022 K562 (specs/169), reads mRNA. Reached and unreached readouts changed
  equally often (0.057 vs 0.057). Most Reactome routes are protein steps, so a transcript readout is a poor match.
- Lun 2019 is genetic, unselected (a whole kinase and phosphatase library), and reads phospho-signalling state, which
  is what our readouts are. Splitting it keeps an independent test.

## The split (fixed now)

- A protein of interest is in half **A** if `sha256("173:" + gene)[0]` is even, otherwise **B**
  (`bench/analysis/173/common_173.py`; `gene` is the Table S1 symbol in `bench/analysis/127/pois.tsv`). Constructs of
  one gene (P1/P2, I1/I2) share its half.
- Panel: 288 genes in A, 348 in B. Population P of specs/127 (6,953 pairs from the committed values-blind predictions,
  sha256 `0c0ec033…`, checked today): **A 3,214 pairs from 198 proteins; B 3,739 pairs from 230 proteins.**
- Control samples (FLAG-GFP, untransfected) belong to both halves: every effect and threshold is computed against
  them, and they carry no protein of interest.
- `bench/datasets.tsv`: `lun2019_a` FIT (development: it may choose and fit rules), `lun2019_b` TEST.

## Unsealing and sealing B (in this order)

1. This spec and the split code are committed and pushed.
2. The archive is unsealed (mode 400) and `complete_dat.csv` of the screening data is extracted to a working directory
   outside the repository. Nothing else is extracted (the time course, binned and pre-analysed files stay in the
   archive).
3. **Format check** (specs/127, unchanged): the header and the distinct values of the sample-identifying columns only.
   Committed to research.md, with the mapping from sample identifier to gene and to control type. If the file cannot
   support specs/127 steps 1–6, the steps are amended and committed before any marker value is read.
4. **Split.** A script committed before it runs streams the file once and writes two files: A = control rows plus rows
   of half-A proteins; B = control rows plus rows of half-B proteins. It prints only row counts per half and per
   control type. A sample that maps to no Table S1 gene and is not a control goes to B (sealed) and is counted.
5. The B file is set to mode 000, the full extract is deleted, and the archive is set back to mode 000. Checksums of
   both halves go in research.md. From then on only the A file is read. The B reader refuses any half-A request, and
   the A reader refuses any half-B gene (as specs/169).

## What is read on half A

The measured call is specs/127's, unchanged (steps 1–7, `bench/lun_score.py`), computed on A's proteins against the
controls. The loader is committed and tested on synthetic data before it reads a marker value.

**Gate before anything else: the positive control** (specs/127 step 7, on A's self pairs). If it fails, scoring stops
and the failure is recorded.

**R1 (primary, specs/127's measure on A):** three-class balanced accuracy on P ∩ A of the frozen model, reachability
and always-NOT; 10,000 bootstrap resamples of proteins; frozen − reachability and frozen − always-NOT as gain, loss or
unresolved. A development reading, not the test of record.

**R2 (the question that decides whether this panel is worth developing on):** does a route matter here?
- Reached = the reachability comparator predicts UP (a route from the pinned roots to a readout node) on P ∩ A.
- Base rate = the share of measured changes (UP or DOWN), among reached and among unreached pairs; enrichment =
  their ratio, with a protein bootstrap interval.
- Lift = precision of the frozen model's change calls minus the reached base rate (specs/167), interval over proteins.
- **Registered reading.** Enrichment's interval above 1.5: routes carry signal on this panel, unlike K562 (1.0), and it
  is fit for development. Interval including 1: it is not, and we say so before building anything on it. Lift's
  interval above 0: the model adds information beyond its routes.

**Also reported:** the whole A panel (unpredicted pairs as NOT), the frozen 92-pathway build, kinases against
phosphatases, and the direction of measured changes among reached pairs (over-expression should mostly raise
readouts downstream of kinases).

## Non-goals

- Half B is not read, summarised or plotted, and no count of its values is made.
- No rule, default or threshold changes here. Anything developed on A is judged later on B, under its own registration.
- The paper's results (Tables S3–S7, figures, results text) stay unread.

## Amendment 1 (2026-11-26, after the format check, before any marker value is read)

**What the format check found** (research.md has the details):
- The deposited per-cell file (`complete_dat.csv`) is **truncated in the deposit**. Its listed size equals the archive
  entry, its last row stops at 24 of 36 fields, and it holds 3 of the 14 plate pairs (about 105 proteins). It cannot
  support the screen.
- The authors' per-sample binned file (`bindat_perbin.csv`) covers all 14 plate pairs: 649 constructs (all Table S1
  markers map, 636 genes) plus 26 FLAG-GFP and 27 untransfected control samples. It holds basal (timepoint 0) and EGF
  10 min, and 30 readouts. Per sample and readout, it has ten GFP bins with the cell count and median per bin.
- From the paper's methods: the bins are ten equal-width bins over each sample's own GFP range (2.5% trimmed). No cells
  are gated as transfected. So bin k is not the same dose across samples, and specs/127 steps 2–4 (dose-matched
  comparison with FLAG-GFP) cannot be built.
- `Data S1` (pre-analysed) is result plots and is not opened. The per-sample statistics file (`bindat_persample.csv`)
  holds the authors' significance calls and is not read.

**The amended measured call** (`bench/analysis/173/lun_a_173.py`, tested on synthetic samples). It asks specs/126's
question directly: the direction of a readout's change as the protein's dose rises.
1. Basal samples only (timepoint 0). A bin counts if it holds at least 20 cells.
2. Dose slope S(s, r) = median of the bin medians in the top three bins − median in the bottom three bins. At least two
   usable bins are needed in each group, otherwise (s, r) is not measured.
3. Control correction: E(s, r) = S(s, r) − the median S of the FLAG-GFP samples on the same plate pair (all FLAG-GFP
   samples if the plate pair has none). This removes what GFP expression itself does to a readout.
4. A gene's E is the median over its samples (constructs and replicates).
5. Threshold δ_r: each FLAG-GFP sample's E against the other FLAG-GFP samples (leave one out), δ_r = the 95th
   percentile of |E|, or the maximum under 10 controls. The call is UP if E > δ_r, DOWN if E < −δ_r, otherwise NOT.
6. The positive control (specs/127 step 7) is unchanged and gates everything.

Readouts in the file but not predicted (p-NFkB, p-mTOR: sites not curated) are measured and not scored. Predicted
readouts the file lacks (cyclin B1, p-RB) are not scored.

**Splitting.** The binned file is split as registered (controls to both halves). The truncated per-cell extract and the
per-sample statistics file are deleted unread beyond the format check. The time-course files were never extracted.
