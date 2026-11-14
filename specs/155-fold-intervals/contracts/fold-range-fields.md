# specs/155 contract: the fold-range fields on POST /api/readouts

## Request

Unchanged. The feature is a server-side lookup; the client sends nothing new.

## Response, per result entry (additive)

A result entry that (1) has a `call`, (2) has a stratum other than
"untouched", and (3) matches at least one covered (axis, class) cell MAY carry:

```json
"measured_fold_ranges": [
  {"axis": "korkut_rppa", "share": 0.5, "lo": 0.31, "hi": 0.92, "n": 214},
  {"axis": "korkut_rppa", "share": 0.9, "lo": 0.12, "hi": 1.4,  "n": 214}
]
```

- `axis`: the record the pairs came from; never pooled across axes.
- `share`: the nominal share of recorded measured folds inside [lo, hi].
- `lo`, `hi`: folds (×), the recorded interval; `lo <= hi`.
- `n`: the recorded pair count behind the cell (always >= the shipped
  `min_cell`).
- The field is ABSENT otherwise — never null, never an empty list.

## Response, protocol object (additive)

```json
"fold_ranges": "ok: 2026-11-xx",          // or the refusal/absence reason
"fold_ranges_min_cell": 30
```

Refusal reasons mirror the specs/146 wording: absent file, unreadable file,
provenance mismatch (catalog build or solver commit), each stated as itself.

## Invariants

- Every existing field of every entry and of the protocol is byte-identical
  with the feature on (FR-001; proven by the specs/141/142 dump replays).
- The served numbers come from the shipped table by lookup only; the server
  computes no quantile, no bin edge, no interval.
- The wording contract for clients: a range is "the middle <share> of
  recorded cases", a share of the record — never a probability.
