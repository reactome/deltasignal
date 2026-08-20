# TCGA LUAD clinical-utility analysis

## Scope

This is a secondary external-association analysis of the completed historical 502-tumor run. It is not perturbation ground truth and it does not represent the later current-stack LNG networks, whose pathway-level aggregation remains unresolved.

DeltaSignal terminal-mean scores were tested continuously. Clinical models adjust for age, sex, pathologic stage, and ever-smoking history. The expression baseline is the mean cohort percentile across unique mapped genes used as DeltaSignal inputs, so position-aware UUID duplication does not reweight genes.

## Results

### PIP3/AKT

- Complete-case cohort: 471 tumors, 168 deaths.
- DeltaSignal after clinical and expression adjustment: HR 0.886 (95% CI 0.713-1.100), p=0.274, three-pathway FDR=0.274 per 1 SD increase.
- Nested test for adding DeltaSignal beyond clinical plus expression: p=0.275.
- Repeated 5-fold CV median concordance: clinical 0.653; plus DeltaSignal 0.654; plus expression 0.650; plus both 0.648.
- Adding DeltaSignal to the expression model changed median held-out concordance by -0.002.
- DeltaSignal/expression-baseline Spearman rho: 0.687; two-predictor VIF: 2.0.

### Mitotic G1/G1-S

- Complete-case cohort: 471 tumors, 168 deaths.
- DeltaSignal after clinical and expression adjustment: HR 0.651 (95% CI 0.392-1.084), p=0.0991, three-pathway FDR=0.149 per 1 SD increase.
- Nested test for adding DeltaSignal beyond clinical plus expression: p=0.1.
- Repeated 5-fold CV median concordance: clinical 0.653; plus DeltaSignal 0.668; plus expression 0.673; plus both 0.671.
- Adding DeltaSignal to the expression model changed median held-out concordance by -0.001.
- DeltaSignal/expression-baseline Spearman rho: 0.960; two-predictor VIF: 12.4.

### Cell Cycle Checkpoints

- Complete-case cohort: 471 tumors, 168 deaths.
- DeltaSignal after clinical and expression adjustment: HR 1.620 (95% CI 1.014-2.589), p=0.0437, three-pathway FDR=0.131 per 1 SD increase.
- Nested test for adding DeltaSignal beyond clinical plus expression: p=0.0397.
- Repeated 5-fold CV median concordance: clinical 0.653; plus DeltaSignal 0.679; plus expression 0.678; plus both 0.676.
- Adding DeltaSignal to the expression model changed median held-out concordance by -0.002.
- DeltaSignal/expression-baseline Spearman rho: 0.949; two-predictor VIF: 9.0.

## Sensitivity and interpretation

Cell Cycle Checkpoints is the prespecified TCGA signal of interest. Its clinical-adjusted estimates across terminal mean, median, lower quartile, and upper quartile are retained in `tables/aggregation_sensitivity.tsv`; disagreement between these summaries is evidence that the pathway-level readout is aggregation-sensitive.

Cell Cycle Checkpoints retains a nominal coefficient after clinical and expression adjustment, but that coefficient does not survive the three-pathway FDR correction. The DeltaSignal and expression scores are strongly collinear, and adding both does not improve repeated-CV concordance over expression alone. The defensible conclusion is that the historical DeltaSignal score captures clinically relevant cell-cycle variation, but this analysis does not establish incremental predictive value from network propagation.

The analysis used 324 pathway-specific unique-gene entries across the three expression baselines. These counts overlap between pathways and must not be interpreted as that many distinct genes overall.

A survival association can demonstrate clinical relevance, but not causal perturbation accuracy. Any in-sample likelihood-ratio result is explanatory; the repeated cross-validation is the more conservative discrimination check. No result here repairs the current-stack terminal aggregation problem.

## Files

- `tables/cox_coefficients.tsv`: continuous hazard-ratio estimates.
- `tables/model_comparison.tsv`: nested tests and in-sample fit.
- `tables/repeated_cv_concordance.tsv`: all repeated-CV results.
- `tables/repeated_cv_summary.tsv`: compact discrimination summary.
- `tables/score_correlations.tsv`: score/baseline/stage correlations.
- `tables/aggregation_sensitivity.tsv`: alternate output summaries.
- `analysis_manifest.json`: revisions, parameters, and input hashes.
