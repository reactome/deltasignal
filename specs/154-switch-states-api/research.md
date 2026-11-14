# specs/154 research: how the second fixed point is found and served

All questions resolved from the code, 2026-11-13.

## R1. The stateless-readouts constraint and the serving shape

**Decision**: `/api/solve` additionally returns `switch_components` — a list,
one entry per detected bistable component: `{id, alt_activities}` where
`alt_activities` is a SPARSE map (uuid → UI-scale activity) holding only the
nodes whose alternate value differs from the canonical solve beyond the noise
floor (the flipped component's members plus its re-solved feed-forward
downstream). `/api/readouts` accepts an optional `switch_components` field
(the same list passed back), resolves every readout once against the
canonical activities (exactly today's path) and once against
`canonical ⊕ alt_activities` per component, through the SAME specs/141
resolution (entity-mean, abundance rule, context filter — one resolution,
two value lookups), and serves `switch`/`state`/`other_state_fold`.

**Rationale**: `/api/readouts` is stateless (the client passes activities
back; `src/api/server.jl` readouts handler reads the request's `activities`).
The readout machinery separates resolution (which nodes) from valuation
(their activities), so valuing one resolution under two vectors is the
designed reuse — the same shape specs/141's entity-mean already uses for
groups. A sparse overlay keeps the payload proportional to the flipped
region, not the network (a full second vector would double every solve
response for the common monostable case).

**Alternatives considered**: a full `alt_activities` vector (payload, and
ill-defined with two bistable components); computing switches purely in
/api/solve per readout (readouts are not known at solve time; resolution
lives in /api/readouts); server-side state (violates the API's stateless
contract).

## R2. Where the second solve hooks, and what it reuses

**Decision**: a post-pass at the end of `solve_scc_ordered!`
(src/solvers/steady_state.jl), where every needed structure is in scope: the
final `x`, per-component reaction lists, `comp_nodes`, `comp_pools`,
`comp_shared`, `managed`, `obs_set`, and — critically — the MUTATED
`IndexedReaction`s (every adopted rule that edits flags or drops edges during
the canonical pass — self-fed holds, carrier drops, own-source removals —
has already been applied to them). For each cyclic component (size > 1, not
fully pinned): copy `x`; set the component's unpinned nodes to the high init
(the internal ceiling 1.0) and, in a second attempt, the low init (0.0; pins
keep their values); re-run the canonical per-sweep sequence confined to the
component — the damped gauss-seidel reaction sweep over the component's
reactions followed by `update_pools!(x, comp_pools[c], comp_shared[c], ...)`
(the exact call at steady_state.jl:2129) — until the component residual is
below the same tolerance or `DS_SWITCH_MAX_SWEEPS` is spent. A converged
point differing from the canonical one beyond the floor (R3) is the second
state. Downstream: re-run the acyclic/feed-forward updates for nodes
topologically after the component (the same single-pass updates the main
loop uses for acyclic nodes), yielding the sparse alternate.

**Rationale**: reusing the mutated reactions and the identical per-sweep
sequence is what makes both basins live under the same rulebook; a separate
"clean" dynamics would find fixed points of a model we do not serve.
Pool-managed components are searched, not excluded: `update_pools!` runs in
the re-solve exactly as in the canonical sweep, and this matters because the
RAF component — the serving bar — carries pools. The specs/076 restart
machinery already demonstrates that re-sweeping a component from an
arbitrary iterate is sound.

**Alternatives considered**: excluding pool components (would exclude RAF —
fails the bar); LM from random starts (specs/076 showed LM's answer depends
on its start; basin inits are the principled starts); refactoring the main
loop into a shared closure (a 700-line high-risk refactor for no behavioural
need — the post-pass duplicates only the ~30-line sweep core).

## R3. The distinctness floor and the readout switch rule

**Decision**: two fixed points are distinct when any component node differs
by more than `DS_SWITCH_FLOOR` = 1e-3 on the internal [0,1] scale (0.1 UI
units). A READOUT is a switch when its two folds fall on different sides of
the served call cutoffs (the 0.85/1.15 the engine already serves), i.e. the
two states would produce different calls — the readout-level rule reuses the
existing cutoff constants, never new magic numbers.

**Rationale**: numeric noise is orders below 1e-3 (specs/013 recorded
relabelling diffs up to 6e-16; the regeneration record's tied-order noise is
the same scale), while the real separations are orders above it (RAF high
state ~0.5 internal vs low ~1e-4, specs/054; the 100× rail vs 0.01× floor is
the full scale). 1e-3 sits four orders above noise and two below the known
real separations. The call-cutoff rule for readouts makes "switch" mean
exactly "the answer class depends on which state you are in".

## R4. Env flags and defaults

**Decision**: `DS_SWITCH_DETECT=on` (default) | `off`;
`DS_SWITCH_MAX_SWEEPS` default 250 (= the specs/076 restart sweep budget,
`DS_NONCONV_RESTART_SWEEPS`' default — the recorded cost of re-sweeping a
component once). Both validated by the config layer like every DS_* flag;
the solve reports `switch_detect`, `switch_components_searched`,
`switch_components_found`, `switch_sweeps_capped`.

**Rationale**: Constitution IV's default-off is for behaviour changes; this
feature's FR-001 parity replay PROVES no behaviour changes. Precedent:
specs/142's `railed` and specs/146's `confidence` shipped default-on as
additive fields. The `off` escape exists for latency-sensitive callers.

## R5. Tests and fixtures

**Decision**: new `test/test_switch_states.jl`, registered in
`test/asserting_suites.txt` and the CLAUDE.md table. Fixtures built fresh
(the existing loop fixtures in `test_loop_elasticity.jl` pin the loop values
they reason about — specs/014's knife-edge nets are close but tuned for
elasticity): a two-node mutual-activation loop with hill gain above the
knife-edge (both basins converge to distinct points; a pin to either side
selects it), the same loop with gain below the edge as the monostable
control (second solve converges back to the canonical point → one point),
and a floor-rejection case (two convergences within 1e-3 → `switch: false`).
No fixture pins the intermediate the detection reasons about
(feedback-fixture-must-compute).

## R6. The coverage report (FR-006)

**Decision**: `bench/analysis/154/` in the specs/140-142 pattern:
`dump_switches.jl` runs in the dev container with the served catalog mounted
(the specs/152 g4 replay pattern, RW gate mount), solving every benchmark
(pathway, perturbation) from the specs/142 rails dump's case list with
detection on, and writing per-readout switch flags; `report_154.py` joins
them against the recorded 3,922 railed folds and prints: the share of railed
folds on flagged readouts; the RAF check (every railed MAPK readout under
RAS/RAF perturbations flagged, high ≈ the recorded ~50×, low near
baseline/zero); the PTEN/MKRN1 PIP3 check (specs/082). The two named checks
are the serving bar; the share is reported whatever it is.

**Label-invariance check (FR-005)**: the detection pass re-run on the
`20261113-1235_baf608b_hgncfix97` build (identical content to the frozen 97
catalog, different uuids — it exists already) must flag the same
(pathway, readout) set with the same states to served precision.
