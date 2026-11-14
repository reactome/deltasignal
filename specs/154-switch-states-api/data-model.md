# specs/154 data model

## SwitchComponent (solve-side diagnostic, served)

One per detected bistable component in a solve response's
`switch_components` list:

| field | type | meaning |
|---|---|---|
| `id` | int | the component's index in this solve (opaque, stable within the response) |
| `size` | int | member node count |
| `alt_activities` | map uuid → float (UI scale) | ONLY nodes differing from the canonical solve beyond the floor: the component's members at the second fixed point plus re-solved feed-forward downstream |

Validation: `alt_activities` is non-empty by construction (a component with
no differing node is not bistable); every value is in [0, 100]; every key is
a node of the solved network.

## Solve diagnostics (additive)

`switch_detect` ("on"/"off"), `switch_components_searched`,
`switch_components_found`, `switch_sweeps_capped` (components whose search
hit `DS_SWITCH_MAX_SWEEPS` — honest absence, FR-004).

## Readout switch fields (readouts-side, additive)

Present only on entries valued with activities (as `railed`):

| field | type | presence rule |
|---|---|---|
| `switch` | bool | always present when the request carried `switch_components` and the readout is valued; absent otherwise |
| `state` | "high" \| "low" | only when `switch: true`; by the readout's own fold at the two states |
| `other_state_fold` | float | only when `switch: true`; same scale/rounding as `fold` |

Readout rule: `switch: true` iff some component's overlay moves the
readout's fold across the served call cutoffs (different call in the two
states). With several bistable components, the first (lowest `id`) that
changes the call provides `state`/`other_state_fold`; the protocol block
says how many components the response carried.

## Protocol block additions (readouts response)

`switch_rule` (the detection description string), `switch_floor` (1e-3
internal), `switch_max_sweeps` (the served cap), `switch_components` (count
received in the request).

## State transitions

None — stateless API; the alternate rides the solve response and is passed
back verbatim.
