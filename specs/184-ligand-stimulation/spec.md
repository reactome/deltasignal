# specs/184: ligand stimulation, the perturbation Reactome's signalling pathways were curated from

**Status:** registered 2026-11-27, before any prediction was made and before any measured value was compared with one.
**Requested by:** Adam, 2026-11-27 ("ok", on the proposal to test ligand stimulation; step 1 of three).

## Why

Every external dataset scored so far perturbs the cell in a way the pathways were not curated from: a drug whose
engagement is uncertain (korkut, specs/138), a kinase over-expressed without evidence that it is active (Lun,
specs/181–183), a CRISPRi knockdown read as transcripts (specs/169). On all of them the model is about as good as
reachability. Reactome's signalling pathways are mostly curated from ligand experiments: add the growth factor and
follow the cascade. A ligand added to a responsive cell is active by construction. No ligand-stimulation axis has
been scored (specs/112 scored only the mTOR-inhibitor sub-challenge of this dataset).

Three earlier findings predict that this is where the model should do best: calls were right far more often when the
perturbation demonstrably engaged (specs/179), short routes carried signal (specs/181), and the one win over
reachability was on this dataset (specs/112).

## Data

HPN-DREAM training data (Hill et al. 2016), already on disk and **FIT** since specs/112. Reverse-phase protein arrays
in four breast cancer lines (BT20, BT549, MCF7, UACC812), serum-starved, then stimulated with EGF, FGF1, HGF, IGF1,
insulin, NRG1, serum or PBS (vehicle) and read at 5, 15, 30 and 60 minutes, 2 and 4 hours; 41–48 phospho readouts
per line (`CSV/<line>_main.csv`); 16–32 unstimulated replicates at time 0 per line. Only DMSO rows are used (no
inhibitor). Because the data are FIT, **nothing here is a held-out test**: this is step 1, a development measurement.
Steps 2 (Lun's EGF well) and 3 (a fresh ligand dataset scored once) are separate registrations.

Read before registering, to fix the design: the file headers, the counts of rows per (inhibitor, stimulus, time
point), and the overall range of values (linear, 0.04–48). No value was looked at by condition or readout.

## Predictions (`bench/hpn_stim_bench.py`, values-blind)

The Lun protocol (specs/127, specs/175) with ligands in place of over-expressed proteins:

- **Build:** Release 98, `20261126-2322_d6d026b_top98` (28 top-level pathways, each network solved on its own), the
  solver at this spec's implementation commit with default settings.
- **Conditions:** EGF, FGF1, HGF, IGF1, INS (insulin), NRG1, each pinned alone at its root forms
  (`pin_occurrences`, `DS_PIN_SCOPE=root_cycle`, `DS_PIN_ROOT_IGNORE=depletion`) at the bench's over-expression level
  (80×); Serum pins all six at once. One solve per condition per network that holds a root form.
- **Readouts:** specs/112's site map (`specs/112-hpn-dream/site_map.tsv`, amendment 4 corrections included), the 47
  readouts with a site. A readout is the mean over the unpinned nodes that are, or contain, an EWAS phosphorylated at
  one of its sites; called with the curator cutoffs (0.85 / 1.15).
- **Reachability:** UP when a readout node is reachable from the pinned roots, else no change.
- **Across networks:** `lun_score.combine_networks` (the call every change-predicting network agrees on; a conflict
  counts as no change).
- **Secondary contrast, reported only (withdrawal):** the same roots pinned at 0 instead; the stimulation call is the
  opposite of the class that solve gives (starved cells have no ligand, so "ligand absent → ligand at rest" is the
  other reading of the experiment).

The model has no cell-line context: the four lines get the same predictions.

## Measured call (`bench/analysis/184/stim_184.py`)

For a cell line, ligand L and readout r:

- d(t) = mean over replicates of log2 value(L, t) − mean over replicates of log2 value(PBS, t), the same time point.
- **E = the mean of d(t) over t in {5, 15, 30, 60 min}** (the early signalling window; at least two time points
  needed). The model is a steady state; the window is fixed here, not chosen from results.
- σ = the standard deviation of log2 value over that line's time-0 unstimulated replicates for r.
- **δ = max(2σ, log2 1.15).** UP if E > δ, DOWN if E < −δ, otherwise no change. log2 1.15 is the model's own cutoff.

Scored pairs: (line, ligand, readout) with a measured call and a prediction in at least one network. Serum is
reported separately and is not in the primary set (serum holds far more than six ligands).

## Registered reading

Three-class balanced accuracy on the six single ligands × four lines, as on Lun. Arms: model, reachability, never.

- **G1 (primary gate):** model − reachability, 95% bootstrap interval above 0. The bootstrap resamples the 24
  (line, ligand) conditions with replacement, 10,000 draws, seed 184. The four lines share one prediction per ligand,
  so the interval understates the uncertainty; a bootstrap over the six ligands alone is reported beside it.
- **Reported, no gate:**
  - route enrichment: the share of reachable pairs that change over the share of unreachable pairs that change, with
    the same interval; the reached base rate against the model's precision (specs/167's comparison);
  - direction: among measured changes the model calls a change, the share with the right direction;
  - recall and precision per arm, per ligand and per line;
  - by engagement, report-only (specs/179's idea): a condition is engaged when the ligand's own receptor readout is
    measured UP (EGF: any of EGFR pY1068, pY1173, pY992; HGF: c-Met pY1235; NRG1: HER3 pY1298 or HER2 pY1248);
    FGF1, IGF1 and insulin have no receptor readout on the array and are filed as engagement unknown; the receptor
    readouts themselves are left out of the split's pairs;
  - the same scores for the late window (2 h and 4 h), for δ doubled, for the withdrawal contrast, and for Serum;
  - the pairs with no prediction in any network, counted.

## Written before the run

- H1. Route enrichment is above 1.5, higher than Lun's (1.5) and below korkut's (3.4).
- H2. Direction accuracy is above 80%.
- H3. The model beats reachability (G1 passes). I hold this at about even odds: the model can only win where it calls
  DOWN correctly or stays silent on a reachable readout that does not move.
- If G1 fails and enrichment is no better than on Lun, the limit is not the perturbation's activity, and step 3 is
  not worth a fresh dataset. If enrichment is clearly higher but G1 fails, the route carries the signal and the
  solver adds none, on the fairest data available.

## Outputs

Predictions under the build's `results/<solver-sha>/`; scoring output outside the repo under `analysis/184/`
(the dataset stays outside the repo); `research.md` here with the figures, whatever they show.
