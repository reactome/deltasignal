# specs/155 data model

## Pair (fitting input; never shipped)

| field | meaning |
|---|---|
| axis | `korkut_rppa` (and `decryptm_phospho` only if R1a admits it) |
| unit | the LODO unit: korkut drug set / decryptM drug |
| readout | the antibody / site identity from the recorded cases file |
| pred_fold | the recorded `pred_ui` (UI fold) of the canonical scoring |
| meas_fold | the dataset's own numeric: 2^log2 per korkut condition; decryptM per R1a |
| stratum | the specs/146 level string, derived with the server's rules |

One pair per (axis, unit, readout, condition). Pairs exist only for scored
cases with a recorded predicted fold; nothing is imputed.

## Calibration table (`src/api/fold_ranges.json`)

Top level: `version`, `spec`, `fit_record` (results path + dataset identity),
`catalog_build`, `solver_commit`, `min_cell` (30), `nominal_shares`
([0.5, 0.9]), `gates` (measured verdicts for both arms, both axes, pass or
fail), `arm` ("a" or "b", per gate 2), `table`.

`table`: axis → stratum level → either
- arm (a): `{share → {lo, hi, n}}` — quantile interval of measured folds, or
- arm (b): `bands`: list of `{fold_lo, fold_hi, share → {lo, hi, n}}` with the
  level's (a) entry kept as `fallback` for folds outside every band.

Every `{lo, hi}` is in fold units (×), `n` is the pair count behind the cell,
and a cell with `n < min_cell` is not written at all.

## Served range object (per readout entry, additive)

`measured_fold_ranges`: list of `{axis, share, lo, hi, n}` — one element per
covered axis × nominal share; the field is ABSENT (not null, not empty) when
the readout is untouched, has no call, or no axis covers its class.

`protocol` additions: `fold_ranges` ("ok: <version>" | the refusal reason),
`fold_ranges_min_cell`.

## State rules

- Table mismatch with the serving record ⇒ no ranges anywhere; the refusal
  named in `protocol.fold_ranges` (specs/143/146 rule).
- Band miss under arm (b) ⇒ the level's fallback (a) entry, never
  extrapolation.
- No table file ⇒ exactly today's response plus the protocol reason.
