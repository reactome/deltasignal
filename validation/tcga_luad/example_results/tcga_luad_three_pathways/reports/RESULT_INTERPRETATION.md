# TCGA LUAD DeltaSignal Validation Results

## What Was Run

- Cohort: 502 TCGA lung adenocarcinoma tumor samples with RNA-seq and survival metadata.
- Pathways: 3 Reactome pathways selected after mRNA observability audit.
- Scoring: DeltaSignal steady-state solve for each sample/pathway pair.
- Sample summary metric: terminal_mean_activity, the average predicted activity of terminal/output nodes in each pathway network.
- Survival test: median split into low/high pathway activity groups, then log-rank test.
- Negative control: 1000 shuffled-label survival controls for terminal_mean_activity.

## Main Finding

Cell Cycle Checkpoints (R-HSA-69620) is the strongest signal in this validation run.

Biologically, this Reactome pathway represents checkpoint mechanisms controlling whether cells continue through the cell cycle: DNA damage checkpoints, G1/S transition control, G2/M control, and mitotic checkpoint behavior. In cancer, high activity in this kind of pathway can indicate highly proliferative tumors or tumors under replication/checkpoint stress, which can be associated with aggressive disease.

For this LUAD run, Cell Cycle Checkpoints had 6552 mRNA-backed root nodes out of 8526 root UUIDs, with 0 missing stable IDs after the LNG mapping fix.

The median split produced two equal groups:

- Low activity group: 251 patients, 73 deaths, KM median survival 1725 days.
- High activity group: 251 patients, 109 deaths, KM median survival 1229 days.

Plain English: patients whose tumors had higher predicted Cell Cycle Checkpoints pathway activity died more often and had shorter median survival in this cohort.

The log-rank p-value was 8.92e-05; BH FDR across the three terminal-mean pathway tests was 0.000268.
The shuffle control found 0/1000 shuffled controls as strong or stronger than the observed split.
Plain English: the observed survival association is stronger than the label-randomized nulls generated for this run.
A univariate Cox model using Cell Cycle Checkpoints terminal mean as a continuous score estimated HR=1.36 per 1 SD increase (95% CI 1.17-1.58, p=6e-05, FDR=0.00018).

## Other Pathways

### Cell Cycle Checkpoints (R-HSA-69620)

- log-rank p: 8.92e-05
- terminal-mean BH FDR: 0.000268
- low/high events: 73/109
- low/high KM median survival: 1725/1229 days

### Mitotic G1-G1 S phases (R-HSA-453279)

- log-rank p: 0.0555
- terminal-mean BH FDR: 0.0832
- low/high events: 81/101
- low/high KM median survival: 1622/1293 days

### PIP3 activates AKT signaling (R-HSA-1257604)

- log-rank p: 0.578
- terminal-mean BH FDR: 0.578
- low/high events: 99/83
- low/high KM median survival: 1492/1622 days

## Caveats

- This validates an association, not causality.
- The current expression-to-observation transform is cohort percentile activity; alternative transforms should be compared.
- Terminal outputs are graph-derived terminals, not manually curated biological readouts.
- Only three pathways were tested in this focused run. Reactome-wide expansion requires strict multiple-testing correction and more careful model selection.
- Clinical covariates beyond survival/event are not yet included in the main result tables.

## Generated Figures

- Kaplan-Meier curves for all pathways: `figures/km/km_all_pathways.png`
- R-HSA-69620 Kaplan-Meier curve: `figures/km/R-HSA-69620_km.png`
- R-HSA-453279 Kaplan-Meier curve: `figures/km/R-HSA-453279_km.png`
- R-HSA-1257604 Kaplan-Meier curve: `figures/km/R-HSA-1257604_km.png`
- Activity by survival event: `figures/activity_by_event/activity_by_survival_event_all_pathways.png`
- R-HSA-69620 activity by event: `figures/activity_by_event/R-HSA-69620_activity_by_survival_event.png`
- R-HSA-453279 activity by event: `figures/activity_by_event/R-HSA-453279_activity_by_survival_event.png`
- R-HSA-1257604 activity by event: `figures/activity_by_event/R-HSA-1257604_activity_by_survival_event.png`
- Terminal mean histograms: `figures/histograms/terminal_mean_histograms_all_pathways.png`
- R-HSA-69620 terminal mean histogram: `figures/histograms/R-HSA-69620_terminal_mean_histogram.png`
- R-HSA-453279 terminal mean histogram: `figures/histograms/R-HSA-453279_terminal_mean_histogram.png`
- R-HSA-1257604 terminal mean histogram: `figures/histograms/R-HSA-1257604_terminal_mean_histogram.png`
- Samples x pathways heatmap: `figures/heatmap_pathway_activity_zscore.png`
- Pair plot of pathway scores: `figures/pairplot_pathway_scores.png`
- Shuffle null plots: `figures/shuffle_nulls/shuffle_nulls_all_pathways.png`
- R-HSA-69620 shuffle null: `figures/shuffle_nulls/R-HSA-69620_shuffle_null.png`
- R-HSA-453279 shuffle null: `figures/shuffle_nulls/R-HSA-453279_shuffle_null.png`
- R-HSA-1257604 shuffle null: `figures/shuffle_nulls/R-HSA-1257604_shuffle_null.png`
- Univariate Cox forest plot: `figures/cox_univariate_terminal_mean.png`
