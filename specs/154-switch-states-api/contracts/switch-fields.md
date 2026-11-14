# specs/154 contract: the additive switch fields

## POST /api/solve (additive)

Response additionally carries, when `DS_SWITCH_DETECT=on` (default) and the
solve converged:

```json
{
  "switch_components": [
    { "id": 3, "size": 41,
      "alt_activities": { "<uuid>": 0.95, "<uuid>": 48.7 } }
  ],
  "diagnostics": { "switch_detect": "on",
                   "switch_components_searched": 5,
                   "switch_components_found": 1,
                   "switch_sweeps_capped": 0 }
}
```

- `switch_components` is `[]` when nothing bistable was found; ABSENT when
  detection is off or the solve did not converge (specs/066: a cap-decided
  point anchors no second-state claim).
- Every existing field is byte-identical to today (FR-001).

## POST /api/readouts (additive)

Request additionally accepts:

```json
{ "switch_components": [ { "id": 3, "alt_activities": { "<uuid>": 48.7 } } ] }
```

- Optional; omitted → today's response byte for byte.
- Malformed (non-list, entry without `alt_activities` map, non-numeric or
  out-of-range values) → the endpoint's existing 400 error shape, naming the
  field.

Response entries valued with activities additionally carry:

```json
{ "switch": true, "state": "high", "other_state_fold": 0.01 }
```

with the presence rules of data-model.md, and the protocol block gains
`switch_rule`, `switch_floor`, `switch_max_sweeps`, `switch_components`.

## Compatibility

- A client that ignores the new fields sees today's behaviour exactly.
- A client that never passes `switch_components` back gets no switch fields —
  and the readouts response is unchanged byte for byte.
- docs/API.md gains both blocks with the presence/absence rules spelled out.
