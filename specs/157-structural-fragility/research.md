# specs/157 research: the ensemble registry, the labels, the serving shape

Resolved from CLAUDE.md's rule record, the measuring specs, and the code,
2026-11-15 — BEFORE any measurement runs.

## R1. The registry and its quantitative bar

**The inclusion bar, registered here**: a variant enters the ensemble only if
it is a recorded measured arm (clause a: a not-adopted near-equivalent; clause
b: the immediately-preceding default an adopted rule replaced) whose recorded
case deltas on BOTH record axes — curator held-out and experimental — are
within ±20 cases of the canonical. The tuning (in-sample curator) figure does
not gate membership (it is FIT-era, and several faithful adoptions moved it
without moving either record axis). Arms that fail the bar are listed with
their numbers; REFUTED arms with large losses are out by the spec's own rule.
Two exclusions by KIND, registered up front:

- **Convergence machinery** (specs/066 fallback, specs/073 pools, specs/076
  restarts): excluded. They change HOW the fixed point is found, not WHICH
  structure is modelled; a flip under them is cap-dependence, which specs/066
  already treats as its own honesty problem (and the solve reports it).
- **Generator-side and bench-side flags** (LNG_*, DS_PIN_*, DS_KO_AGG,
  DS_TOTAL_READOUT): excluded. The ensemble runs on the SERVED catalog with
  the recorded scoring protocol; a variant must be a solver `DS_*` env the
  served container can carry.

**IN (20 variants)** — each row: setting · measuring spec · recorded deltas
(held-out / experimental unless noted):

| # | variant | spec | recorded record-axis deltas |
|---|---------|------|------------------------------|
| 1 | `DS_DEPLETION_OWN_PRODUCT=suppress_only_entity` | 060 | held-out 0 / exp −17 |
| 2 | `DS_DEPLETION_OWN_PRODUCT=full` | 061 (pre-adoption default; bundle re-read) | the 061 bundle's +32/+7 is the THREE-part bundle vs old defaults; `full` alone vs adopted, holding twins and pins: resolved as within-bar from specs/060's suppress_only_source-vs-full reading (held-out 0 alone); accepted under clause (b) with the caveat recorded |
| 3 | `DS_DEPLETION_COMBINE=share` | 103 | held-out 0 / exp 0 (PIP3 −13 inside tuning) |
| 4 | `DS_SET_POOL_MODE=mean` | 033/038 | "mean does not gain" — ~0 / ~0 |
| 5 | `DS_UBIQUITIN_MODE=inert` | 098 | held-out +6 / exp −4 |
| 6 | `DS_POOL_SUPPLY_AND=on` | 090 | held-out 0 / exp 0 |
| 7 | `DS_POOL_DRIVE_TWINS=hold` | 104 | no class change on either axis |
| 8 | `DS_POOL_BASE_WEIGHT=full` | 091 (pre-adoption default) | held-out 0 / exp +1 |
| 9 | `DS_CARRIER_LOOP_MOVES=keep` | 093 | held-out −10 / exp 0 |
| 10 | `DS_CARRIER_LOOP_OPEN=fed_single` | 088 (pre-adoption default) | held-out 0 / exp +1 |
| 11 | `DS_CARRIER_COMPONENT_TWIN=off` | 084 (pre-adoption default) | held-out 0 / exp +4 |
| 12 | `DS_METABOLITE_POOLS=on` | 107 rule 1 | no class change on either axis |
| 13 | `DS_SELF_FED_MODE=catalyst_loop` | 109 | no class change on any axis |
| 14 | `DS_SELF_ACTIVATOR_LEAVES=1` | 136 arm L | no call changes |
| 15 | `DS_SELF_ACTIVATOR_WEIGHT=off` | 083 (pre-adoption default) | no class change (adopted on faithfulness) |
| 16 | `DS_SELF_INHIBITOR_COPIES=off` | 095 (pre-adoption default) | held-out 0 / exp +9 |
| 17* | `DS_CYCLE_REEVAL_CARRIERS=full` | 075 (pre-adoption default) | canonical identical (0 pairs) |
| 18* | `DS_RNA_TWIN=off` | 100 (pre-adoption default) | held-out 0 / exp +4 |
| 19* | `DS_DRUG_MODE=propagate` | 055 (pre-adoption default) | held-out +8 / exp +2 |
| 20* | `DS_PATHOGEN_MODE=propagate` | 048 (pre-adoption default) | held-out +1 / exp +7 |

(*) rows 17–20 numbered through 20; the table is the registry's source; the
frozen machine-readable copy is `registry.tsv` beside this file, committed
before the run (FR-001).

**OUT, with the number that excludes each** (beyond the two kind-exclusions):

- specs/063 sequestration (all forms): held-out −1,121 to −1,686 — refuted.
- specs/056/132/133 context arms: held-out −164 and worse; a different
  catalog besides — refuted and out of kind.
- specs/040 rule A `DS_SELF_FED_MODE=entry`: TP53 −190, PIP3 −126 — refuted.
- `DS_SET_POOL_MODE=max`: held-out −120 / exp −57 — refuted.
- `DS_SET_POOL_MODE=geomean`: no recorded measurement found in specs/033/038
  — not a RECORDED arm, so not admissible (the registry holds measured arms
  only).
- `DS_SELF_INHIBITOR_LEAVES=0`: +48 / +20 — over the bar.
- `DS_SELF_INHIBITOR_WEIGHT=off`: specs/012's −61 arm was hash-seed dependent
  and "must be re-run before citing" (CLAUDE.md) — no citable number, out.
- `DS_RECYCLED_TWIN=off`: only measured inside the specs/061 bundle (+32/+7
  with TGF-beta −12 traced to exactly this flag) — over the bar on held-out.
- `DS_CYCLE_MODE=off`: +77 held-out — over the bar.
- `DS_CARRIER_LOOP_MODE=off/closed/all`: +122/−2 (off), ~+41 to twin
  (closed), −55 (all) — over the bar.
- `DS_ZERO_SUPPLY_MODE=hold`: held-out −20 (p 1e-5) is AT the bar but the
  spec records it as a traced systematic artefact (pins as reach barriers),
  not a near-equivalent — out, by the spec's own diagnosis.
- `DS_CARRIER_MULTISTEP=strict/wide/compartment`: strict's +6 is free
  ubiquitin (duplicating row 5's mechanism); wide −50 over the bar;
  compartment −10 within the bar but needs `leaf_refs.csv` columns the served
  build carries — resolved at registry freeze: strict and compartment are
  admissible but REDUNDANT with rows 5 and each other's mechanism; kept OUT
  to avoid double-weighting one contested choice (ubiquitin) — recorded as a
  judgment, revisitable.
- `DS_COMPOSITION_MODE=limit/limit_novel`, `DS_SCC_BREAK_ROLES`: measured on
  pre-freeze catalogs with welded loops (specs/016/018 era) — the recorded
  numbers do not transfer to the served build; out as not-currently-recorded.
- specs/017 pool modes (`minimize/pool/pool_all`): refuted or (minimize)
  −71/−87/−148 — refuted.
- specs/110 `DS_DEDUP_ACTIVATORS=1`: held-out −60 — over the bar.

## R2. Labels and the non-convergence rule (registered before the run)

Per (case, readout): the CANONICAL call (the recorded one) plus one call per
variant. ROBUST = every contributing call identical to the canonical; FRAGILE
= any contributing call differs. Contribution rule: a variant solve that does
not converge for a case contributes NO call for that case (specs/066: a
cap-decided state anchors no claim); the per-case count of non-contributing
variants is recorded. A (case, readout) where fewer than 80% of the
ensemble's variants converge is EXCLUDED from the gate population and
counted in the report. The canonical side always converges on the record
(the recorded benches report it), so the canonical call is always present.

## R3. The serving shape (resolved from the code)

`resolve_reaction_eval_config()` reads process ENV once per solve with no
per-request override, and the structural rules (twins, carrier drops, holds)
are likewise env-resolved at solve setup — a live per-request ensemble would
need env mutation in the shared served container, which is thread-unsafe and
rejected (the spec's own rejection clause). **Serving, if the gate passes, is
a precomputed per-class fragility table**, the specs/146 confidence pattern
the codebase carries twice (confidence.json, fold_ranges.json): per
(confidence stratum, call) class, the share of the recorded ensemble
population in that class labelled fragile, with n; stratum-walk lookup,
registered minimum cell 30, provenance refusal on catalog mismatch, honest
absence for untouched readouts and silent classes. `/api/readouts` serves
per-entry `fragility {rate, n}` plus a top-level table-status string — all
additive. (Exact field names fixed at implementation; the dumps replay
unchanged either way.)

## R4. The measurement runner

`bench/analysis/157/ensemble_157.py`: for each registry variant, a replay
server (`docker compose run` with the canonical serving env PLUS exactly one
variant setting, build pinned to the serving record) whose variant env is
verified POSITIVELY in-container before any measurement — each `DS_*` here
has a diagnostics fingerprint in the solve response (depletion_combine,
pool_supply_and, metabolite_pools_rule, self_fed_rule, …): the probe asserts
the variant's own diagnostic, never the env. Population: the recorded
curator and korkut case lists on the served 98 record, replayed exactly as
the specs/154/156 dumps replay them (pins from the tracked g4 dumps for
curator; the korkut bench machinery for drug sets). Up to 4 servers, variants
sequenced across them; ~20 × ~4,200 solve-cases ≈ 84k solves ≈ 2–3 h at the
specs/156 re-run's measured rate. Output: one TSV per variant
(case keys, variant call, converged), committed whatever it says. Scoring:
`score_157.py`, pure python — labels per R2, per-axis robust/fragile accuracy,
95% bootstrap over perturbations (curator) and drug sets (korkut), 10k
resamples, seed 157; RAF set aside read beside the headline.

## R5. Tests

If served: `test_fragility.jl` (registered, CLAUDE.md row) on the specs/146
suite's pattern — the lookup (class present/absent), the minimum-cell
refusal, the provenance refusal, additive parity on the existing fixture.
The scorer is fixture-checked on computed known-answer TSVs before the real
run (the specs/156 scorer discipline).

## The measurement (2026-11-15) — GATE PASSES

All 21 runs (canonical + the 20 registry variants) completed on the served 98
build through per-variant replay servers, every probe passing (variant live,
every other registry env canonical). Self-check: the canonical run reproduces
the recorded call on 34,789 of 34,825 truth-joined pairs — curator and
experimental exactly (24,356/24,356), korkut 10,433/10,469 (36 borderline
valuation differences between /api/readouts and the bench's own korkut
reduction; internally consistent since every variant is compared against the
canonical RUN's calls, same valuation both sides).

34,820 labelled pairs (5 excluded below the 80% contribution floor), scored by
`score_157.py` (10k bootstrap, seed 157):

| axis | robust acc (n) | fragile acc (n, share) | delta [95% boot] | verdict |
|---|---|---|---|---|
| curator | 0.9068 (23,058) | 0.6250 (448, 1.9%) | **+0.2818 [+0.1298, +0.4573]** | PASS |
| — RAF aside | 0.9073 | 0.6038 | +0.3035 [+0.1428, +0.4888] | — |
| experimental | 0.7776 (769) | 0.5921 (76, 9.0%) | +0.1855 [−0.0192, +0.4183] | fail (direction holds, small n) |
| korkut | 0.8303 (10,216) | **0.1976** (253, 2.4%) | **+0.6326 [+0.5441, +0.7048]** | PASS |

**The gate passes on two of three axes.** The informativeness statement: a
fragile korkut call is right one time in five; a robust one five times in six.
Fragility is RARE — the model's answers are structurally stable almost
everywhere — and where it appears it marks exactly the contested-choice
dependence the flag exists to name.

**Flip anatomy** (a pair may flip under several variants): own_entity 326,
depl_share 116, open_fed_single 114, drug_propagate 91, ub_inert 79,
rna_twin_off 77, base_full 54, pool_supply_and 40, own_full 12,
comp_twin_off 11, setpool_mean 8, pathogen_propagate 7, drive_twins 4 — and
SEVEN variants with zero flips (metab_pools, catalyst_loops, act_leaves,
act_weight_off, inh_copies_off, reeval_full, moves_keep), exactly the arms
their specs recorded as no-class-change, which is itself a consistency check
of the whole pipeline.

**Serving decision (per the registration)**: the gate passes → the flag
serves, as the R3 shape (per-class fragility-rate table, specs/146 pattern,
min cell 30, provenance refusal). P4 proceeds.

**R3 refinement, fixed before the table is fitted**: the serving table's
classes use exactly the features the readouts endpoint knows per entry AND the
ensemble record carries per pair — call, magnitude bin (the specs/146 bin on
the canonical fold), converged — with the level walk full (call|mag|conv) →
no_converged (call|mag) → call_only (call), min cell 30, pooled across the
three record axes (per-axis counts kept in the table's validation block). The
confidence table's basis/railed features are not in the ensemble record and
re-running 21 arms to add them buys nothing the gate needs; recorded here as
the registered shape before any fitting. An untouched readout (fold exactly
1.0) serves no fragility, the specs/146 amendment-1 rule.

**Additive parity (P4, 2026-11-15)**: the g4 replay through the feature branch
with the fragility table LIVE: all 24,356 recorded folds reproduce (max |diff|
5e-7, zero past 1e-5; `parity_folds_157.tsv`). The flag adds fields; nothing
existing moves.

**Provenance note (2026-11-20, code review): `solver_commit: 7fe3dcc` in
`src/api/fragility.json`.** `7fe3dcc` is not on main, and that is expected:
it is the P3 measurement commit on PR #261's branch
(`origin/specs/157-structural-fragility`, still on the remote), and #261 was
squash-merged as `bbcc44b2`, which replaced the branch's commits with one. The
ensemble the table was fitted from ran at `7fe3dcc`. Checked: `src/solvers`,
`src/core` and `src/io` are identical at `7fe3dcc`, at `bbcc44b2` and at
today's main (`43066b15`); the only `src/` difference between `7fe3dcc` and
`bbcc44b2` is P4's serving code (the table itself, `readouts.jl`,
`server.jl`). The table is unchanged; the recorded commit names the code that
produced it, and the solver it describes is the one main serves.
## Amendment, 2026-11-20 (code review): the untouched rule read real folds exactly

The fitter (`fit_157.py`) and the engine tested `fold == 1.0` on folds replayed straight from solves,
which carry float noise (1.0000000000000218). The recorded case files print 6 figures, so there an
untouched readout reads exactly 1. Of 22,975 record-untouched pairs in `runs/calls_canonical.tsv.gz`,
16,576 had a replayed fold that was not exactly 1.0. The shipped table's moved-but-normal cell
(`not|le1|True`, n 17,312, rate 0.0064) was 96% untouched readouts.

Fix: one rule, the fold rounded to 6 significant figures (`bench/analysis/fold_rules.py`, engine
`record_fold` / `readout_untouched` in `src/api/readouts.jl`), for untouched, railed and the magnitude
bin. Real-fold check (`bench/analysis/review_1120/untouched_parity.py`, replay vs record, 34,825 joined
pairs): untouched mismatches 16,576 → 67; railed mismatches 71 → 35. The residual 67 + 35 are all
korkut, and they are not rounding: 394 korkut replay folds differ from the record beyond print
precision, because this replay sent korkut phospho readouts as `forms` queries (first resolving
stable id) while `korkut_bench.py` takes the union of a readout's sites over unpinned copies. Filed as
its own issue; it touches the korkut side of the labels and the gate below.

Refitted table (`src/api/fragility.json`, version 2026-11-20): 11,780 pairs (curator 9,141,
experimental 707, korkut 1,932; 23,040 untouched skipped). The moved-but-normal class is now
`not|le1|True` 0.1457 (n 693), not 0.0064. Other cells move by at most 0.005. The labels and the
registered gate do not use folds, so the gate result above stands. The specs/165 figures read the
recorded (rounded) folds and stand.
The refitted table keeps `solver_commit: 7fe3dcc` (the code the ensemble runs were solved with, note above),
passed explicitly: `fit_157.py runs src/api/fragility.json 7fe3dcc`.

## Amendment, 2026-11-21 (issue #278): the korkut leg re-run with the record's valuation

**Why.** The ensemble's korkut leg sent phospho readouts to `/api/readouts` as `forms` queries (the
first stable id that resolves, every copy including pinned ones). `korkut_bench.py` values a phospho
readout over the union of its sites' ids, unpinned copies only. 442 of the 3,565 phospho folds in the
canonical run differed from the record at print precision (all 6,904 total-readout folds already
matched under `abundance`). That was the source of the 36 canonical call mismatches in the self-check
above, and of the residual 67 untouched / 35 railed label mismatches in the 2026-11-20 amendment.

**Fix.** `/api/readouts` gains `mode: "sites"` (korkut_bench's union, own filter inside it, no
abundance cut, mean over unpinned copies, the abundance mode's `status` / `occurrences`); `forms` and
`abundance` are unchanged. `ensemble_157.py` sends korkut phospho readouts as `sites`
(`--korkut-phospho-mode forms` reproduces the 2026-11-15 runs).

**Parity gate** (`bench/analysis/157/korkut_parity_157.py`, replay vs the recorded
`results/cff961d/r98fix_korkut/korkut_cases.tsv`): **10,469 of 10,469** korkut folds and calls reproduce at
print precision (phospho 3,565 / 3,565, total 6,904 / 6,904). On the stored canonical run,
`untouched_parity.py` now finds 0 untouched and 0 railed mismatches over all 34,825 joined pairs (from
67 and 35).

**Re-run.** All 21 runs (canonical + the 20 registry variants), korkut axis only, at `a8b69661`, on the
served 98 build (`20261113-r98fix-plus-leads`, pathway directories byte-identical to the record's
`20261113-1311_baf608b_r98fix`), each in its own replay container with the dev compose environment plus
exactly that variant's override, every env probe passing (variant live, every other registry rule
canonical). Two servers at a time, 1,469 solves per run, 10–15 min per run, about 2 h 45 min wall.
Curator and experimental rows do not depend on the korkut readout mode and are carried over verbatim
from the 2026-11-15 runs (`merge_korkut_278.py`, which refuses unless every variant's korkut key set
is identical). The merged record is `runs_278/` (calls gzipped; unpack before `score_157.py`, which
reads `calls_*.tsv`). `runs/` is kept as the 2026-11-15 record, and re-scoring it reproduces the
figures above exactly.

**Gate** (`score_157.py`, same registration, seed 157; `runs_278/score_278.txt`):

| axis | 2026-11-15 (forms) | 2026-11-21 (sites) |
|---|---|---|
| self-check | 34,789 / 34,825 (36 korkut mismatches) | **34,825 / 34,825** |
| curator | +0.2818 [+0.1298, +0.4573] PASS | unchanged |
| experimental | +0.1855 [−0.0192, +0.4183] fail | unchanged |
| korkut robust acc (n) | 0.8303 (10,216) | 0.8268 (10,207) |
| korkut fragile acc (n) | 0.1976 (253) | 0.2252 (262) |
| korkut delta | +0.6326 [+0.5441, +0.7048] PASS | **+0.6016 [+0.5024, +0.6796] PASS** |

**The gate still passes on two of three axes**, and the conclusion stands: a fragile korkut call is
right about one time in four and a half (was one in five), and a robust one five times in six. Flips by variant change only on
korkut: own_entity 326 → 339, depl_share 116 → 110, drug_propagate 91 → 93; the other variants are
unchanged.

**Label changes** (korkut only, 34,820 labelled pairs either way, the same pairs): 25 robust → fragile,
16 fragile → robust; 36 canonical calls change (now equal to the record), of which 34 were right only
under the `forms` valuation and 2 were wrong both ways.

**Table refit** (`src/api/fragility.json`, version 2026-11-21; `fit_157.py runs_278 src/api/fragility.json
"7fe3dcc (curator, experimental) + a8b69661 (korkut)"`, because the two parts of the record ran at two
commits). Between them `src/solvers` / `src/core` changed only in switch detection (#279) and
inverse/soft-observation guards (#280), which do not move a canonical fixed point. The korkut canonical
run at `a8b69661` reproduces the `cff961d` record exactly. The same 17 cells remain. 11,847 pairs
(korkut 1,932 → 1,999; 22,973 untouched skipped, was 23,040).

| cell | 2026-11-20 rate (n) | 2026-11-21 rate (n) |
|---|---|---|
| down\|gt3\|True | 0.0270 (3,551) | 0.0276 (3,550) |
| down\|le1\|True | 0.1339 (1,479) | 0.1168 (1,404) |
| down\|le3\|True | 0.0993 (937) | 0.1085 (1,023) |
| not\|le1\|True | 0.1457 (693) | 0.1531 (725) |
| up\|gt3\|True | 0.0265 (4,340) | 0.0265 (4,340) |
| up\|le1\|True | 0.2181 (188) | 0.2687 (201) |
| up\|le3\|True | 0.2078 (592) | 0.2036 (604) |
| call_only down / not / up | 0.0649 / 0.1457 / 0.0545 | 0.0624 / 0.1531 / 0.0568 |

(no_converged cells equal the full cells: every labelled pair converged.) The served `fragility` rates
move by up to 0.051 (`up|le1`). No class appears or disappears.
