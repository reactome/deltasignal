Diagnostics for specs/054 (RAF rail). Each reads `/s/raf_cases.tsv`, the RAF
rows of a case log (axis, gene, direction, readout, expected, predicted,
pred_ui, pinned uuids, readout uuids), made from cj0's `*_cases.tsv` by:

    awk -F'\t' -v ax=$ax 'NR>1 && $1=="RAF_MAP_kinase_cascade" && $7==1 {print ax"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5"\t"$11"\t"$12"\t"$13}'

Run in the dev container with the build at /catalog and that file's directory
at /s (see specs/054 research.md for the build and results used).

T5 (loop gain, the threshold): `raf_gain.jl`, `raf_sweep_cap.jl`,
`raf_thresh.jl`, `raf_sweep_fine.jl`. `raf_gain.jl` needs only the build;
the others read `/s/raf_cases.tsv` for the pinned and readout uuids.
