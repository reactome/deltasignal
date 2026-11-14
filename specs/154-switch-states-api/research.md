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

## The FR-006 coverage measurement: the bar FAILS as registered, and the
## failure is a finding (2026-11-13)

The full pass over the benchmark record (every (pathway, gene, direction)
solved through the SERVED stack, `bench/analysis/154/dump_switches.py`,
`report/switch_cases.tsv`):

- **Parity re-proven end to end**: all 24,356 recorded folds reproduced
  through parse→pins→solve→readouts with detection on — zero mismatches.
- **605 of 3,922 railed folds (15.4%) flag as switches.** The RAF check fails
  as written: knockouts flag (ARAF/BRAF/HRAS/KRAS/NRAS KO — the loop, freed
  of the perturbed paralog's drive, holds both states, and the served states
  match the specs/054 record), but over-expressions do NOT, and neither do
  the PTEN/MKRN1 cases. Traced, not assumed: the unflagged components WERE
  searched, nothing hit the sweep cap, and both basin re-solves converged —
  to the same point.
- **The finding**: a rail under a FORCING pin is not a two-state situation.
  Under HRAS at 80x the cascade's low state is not a fixed point (the drive
  carries any start back up); under a PTEN pin in either direction the MKRN1
  loop is slammed to one side. specs/054's own threshold says why: the flip
  point is ~3.5% from baseline, so any benchmark-dose perturbation of the
  switch's own driver forces it. The recorded rails therefore split into two
  classes this feature now distinguishes mechanically: **free-switch states
  (15.4%, served as switches with both states) and perturbation-forced
  saturations (84.6%, correctly served as railed, not switches)**. The
  registered bar ("every railed RAF MAPK readout under RAS/RAF perturbations
  is flagged") conflated them; the detector separated them.
- RAF1 knockouts (3 folds) also do not flag: unlike ARAF/BRAF, the RAF1 copy
  is structurally required by the remaining loop's complexes, so its KO
  leaves no free bistable component. Consistent with the same rule.

## Consequence (the registered one)

The serving bar fails, so the feature does not serve on this spec's
authority. The proposed amendment, Adam's to accept or reject:

- Bar (a): the detector finds both RAF states where the component is free
  (the KO cases) with states matching specs/054 — PASSES on this measurement.
- Bar (b): no false second state anywhere (the monostable control and the
  searched-and-converged trace above) — PASSES.
- Bar (c): forced rails keep today's honest `railed: true, switch: false` —
  by construction, with the 15.4% / 84.6% split documented and served
  semantics stated in docs: "switch means two states exist UNDER YOUR
  CONSTRAINTS; a forced rail has one".

Everything else stands regardless: FR-001 parity (unit fixture, the two
dump replays, and the 24,356-fold end-to-end reproduction), the 47-assert
suite, label-invariant detection by construction (sorted component order,
no uuid-dependent choices — the identical-content rebuild check runs with
the serving decision).

## The informativeness check (registered 2026-11-14, before measuring)

Adam: "can you run something to show if this works better or not?" The
feature cannot move accuracy (parity-proven); the measurable claim is that
the split it serves is INFORMATIVE: free-switch rails and forced-saturation
rails should match the recorded ground truth at different rates. Question,
registered before looking: on the benchmark record's railed calls, what is
the share correct (vs curators on the curator axis, vs experiments on the
experimental axis) for switch=true vs switch=false cases? If the rates
differ, the served distinction carries information a user (and a finer
specs/146 confidence stratum) can use; if they are equal, the flag is
honest but uninformative and the PR should say so. Descriptive stratum
rates on recorded data — the specs/146 footing; no call changes, no gate.

**Measured (same day):** on the record's railed calls, free-switch rails are
right 95.4% vs forced-saturation rails 93.2% against curators (+2.3pp,
z = 2.03 — marginal), and 94.1% vs 95.1% against experiments (−0.9pp,
z = −0.23, n = 34 — nothing). Verdict, honestly: the flag does NOT
meaningfully separate hit rates. Both rail classes are highly reliable
DIRECTION calls; the distinction's value is purely semantic — what the
served number MEANS (a state of two, or the only state there is) — not
statistical. A finer specs/146 confidence stratum on this split would buy
almost nothing. The serving decision should rest on the semantics alone.

## Serving decision (Adam, 2026-11-14)

The amended bar is accepted: "lets keep it. I am not sure of how you
implemented it exactly but I am sure we can modify it later." Merged with the
informativeness verdict on the record (the flag's value is semantic, not
statistical). Later modification is cheap by design: DS_SWITCH_DETECT=off
disables it whole; the floor and sweep cap are env-served; the detection is a
self-contained post-pass.

## The FR-005 label-invariance check (run 2026-11-14, after serving, as specified)

Detection re-run on two Release-97 catalogs of identical content and different
uuids — the frozen catalog `20261026-1349_8b6527e` and the rebuild
`20261113-1235_baf608b_hgncfix97` — through the merged, served engine
(`bench/analysis/154/dump_switches_inv.py`, which re-resolves pins per catalog
from the specs/140 dump's stable_ids; `compare_inv.py`; the specs/140-142
dumps supply the 24,356-fold case list).

- **The registered check PASSES**: both catalogs flag the identical
  (axis, pathway, gene, direction, readout) switch set (2,300 flagged folds,
  zero flag differences), with identical states (zero differences), and the
  canonical folds agree to 1.3e-12 (both reproduce the recorded folds to the
  dump's 5e-7 rounding).
- **One traced, documented finding on the alternate COORDINATE**: 10 of the
  2,300 flagged folds — all of them SMAD2 and SMAD3 knockouts in TGF-beta —
  serve an `other_state_fold` that differs between the catalogs by up to 0.5
  at ~50× (≤1% relative), and the values SWAP between the two genes'
  perturbations (frozen SMAD3 reads what rebuild SMAD2 reads, and mirrored).
  The mechanism: SMAD2/SMAD3 are structurally symmetric paralogs, so the
  component holds two near-degenerate alternate fixed points ~1% apart, and
  which one a basin re-solve lands on follows the sweep order's tie-break —
  which is label order exactly where structure is tied (the specs/013 record:
  "tied positions keep label order"). The switch's existence, its state
  labels, and the canonical answer are unaffected; only the second state's
  last 1% of magnitude is, in exactly the component class where the two
  states' coordinates are near-mirror-images anyway.
- Served honestly as is: the note's claim is "two states exist and yours is
  the low one; the other reads ~N×" — a ≤1% label dependence on N in one
  symmetric pair does not change any claim the panel makes. Recorded here so
  a future exact-replication harness knows to compare the alternate at 1%
  tolerance, not bit for bit.
