# specs/156 research: the reconciliation's form, resolved from the solver

Resolved from src/solvers/steady_state.jl, 2026-11-14.

## R1. The iteration form (FR-002)

**Decision**: a soft observation on node t (internal value y, weight w)
redefines t's update target everywhere the solver computes it: wherever a
site computes `fwd = compute_reaction_output_vec(x, r, ...)` for t, the value
written/tested becomes `g = w*y + (1-w)*fwd`. The fixed point then satisfies
`x_t = w*y + (1-w)*F_t(x)` — weight 1 is the hard pin's equation, weight 0 is
today's, and between them the blend is monotone in w for fixed F (the fixture
computes it through a chain where F also moves, asserting monotone end
values). The sites, exhaustively (the pin-skip census `t in obs_set`):

- the acyclic single-pass write (steady_state.jl ~2054) — blend and write;
- the damped component sweep (~2166) — blend `fwd` BEFORE the damping and
  BEFORE the residual test, so the undamped residual `|g - x_t|` measures
  self-consistency of the soft-observed system (the convergence verdict and
  the stopping rule stay on one scale, the recorded lesson at that line);
- `_component_residual` (~1140) and `_component_residual_pools` (~1264) — the
  final consistency check tests `|g - x_t|`;
- the LM fallback's `ComponentProblem` residual rows (~854): a soft node's
  row becomes `g(x) - x_t` (the row stays smooth in x; w, y are constants);
- the specs/154 switch-detection sweep (`_sweep_component!`): the same blend,
  so an alternate stable state also honours the soft observation.

A soft-observed ROOT (no producing reaction) has F = its resting value, so it
is set once to `w*y + (1-w)*baseline_t` — the same formula with the only F a
root has. Seeding: a soft-observed node seeds at its blend-with-baseline (the
hard-pin seeding precedent), which the iteration then refines.

**Rationale**: this is the smallest rule that makes the served fixed point the
spec's reconciliation under EVERY solve path (sweep, pools' neighbours, LM,
restarts, switch re-sweeps), with convergence honestly measured against the
blended system. It is specs/003 FR4's intent without the minimiser: the
`Σω(x−y)²` term's stationary condition under a quadratic model-consistency
term is exactly a weighted pull toward y, but the fixed-point form needs no
objective and inherits every adopted rule unchanged.

**Alternatives considered**: implementing FR4 inside `DS_SCC_METHOD=minimize`
(the recorded −71/−87/−148 held-out loser; soft pins would inherit a worse
solver); post-hoc blending after the solve (downstream nodes would not read
the reconciled value — the whole point is that they do); reinterpreting the
existing confidence (the registered trap: every recorded pin says 0.9).

## R2. Precedence and scope (the edge cases)

- A node both hard-pinned and soft-observed: the hard pin wins (soft entries
  are filtered against `obs_set` up front) and the response lists the node
  under `soft_overridden` with reason `pin`.
- A POOL-MANAGED node (`t in managed`): its value is written by
  `update_pools!`, and v1 does NOT blend inside the pool solve — the node is
  listed under `soft_overridden` with reason `pool`, never silently averaged.
  Blending inside the cycle balance is its own design (the pool's stationary
  distribution with an observation term) and is out of scope until the arm
  shows conditioning carries signal at all.
- Cofactor/drug/pathogen holds already outrank hard pins' absence the same
  way they do today; a soft observation on a held node is `overridden: held`.
- Duplicate soft entries with different values: 400, as hard contradictions
  are. Weight outside [0,1] or value outside [0,100]: 400.

## R3. The API and env surface

- `POST /api/solve` gains optional `soft_observations`:
  `{uuid: [value_ui, weight]}` (value on the 0–100 scale like observations;
  weight in [0,1]). Present while `DS_SOFT_OBS=off` (the default): a 400
  naming the flag — a client must never believe it conditioned a solve that
  ignored it. `DS_SOFT_OBS=blend` enables.
- CLI `solve` accepts a `--soft-observations` CSV (node, value, weight).
- Diagnostics: `soft_obs_rule` ("off" | "blend"), `soft_observed` (count
  applied), `soft_overridden` (list of {uuid, reason}). With no soft
  observations in the request the solve path is UNTOUCHED (FR-001's
  bit-identity is by construction: the soft map is empty and every site's
  check is `isempty`-guarded).

## R4. The registered arm (FR-004) — how it runs

`bench/analysis/156/arm_156.py`: for each korkut (drug set, readout) scored
case of the canonical record, build the solve as the bench does (hard pins
from the recorded pin uuids on the serving catalog), add soft observations
for the drug set's OTHER measured readouts — each at its measured fold
(2^log2, the specs/155 pair machinery; readout → its resolved node set via
/api/readouts groups, soft-observing each unpinned direct node), weight w —
and read the held-out readout through /api/readouts. Weights {0.25, 0.5,
0.75}, registered, no fitting. Score balanced accuracy of the held-out calls
vs the recorded measured calls; baseline = the same cases' unconditioned
recorded calls. Gate: some registered weight beats baseline with a 95%
bootstrap interval over drug sets excluding zero. Role guard: korkut2015 FIT.
The engine runs with `DS_SOFT_OBS=blend` via run_arm.sh (`--server`), on a
pinned worktree, per the house arm rules.

## R5. Tests

New `test/test_soft_observations.jl` (bare name in asserting_suites.txt,
CLAUDE.md row): the US1 chain fixture (weight-1 equality with the hard pin,
weight-0 equality with no observation, monotone between, downstream reads the
blend), the root rule, precedence (`pin`/`pool` overrides reported), the 400s
(flag off, bad weight, contradiction), FR-001 bit-identity (same request,
flag on vs off, no soft entries), and a cyclic-component fixture where the
blend's fixed point is computed by hand. `test_api_errors.jl` gains the
request-shape cases.

## FR-001 parity, end to end (P3b, 2026-11-14)

The whole recorded benchmark (the tracked 98-record dumps,
`specs/152-release-98-redo/g4/`) replayed through the feature branch's served
stack with `DS_SOFT_OBS=blend` LIVE in the container (verified positively: a
soft-carrying solve succeeds and reports `blend`/1 — never inferred from the
env alone, the container-env rule) and NO soft entries in the replay's
requests: **all 24,356 recorded folds reproduce to the dump's rounding (max
|diff| 5e-7, zero past 1e-5)** (`parity_folds_156.tsv`). The flag alone
changes nothing; only a request's own soft observations do — exactly FR-001.
