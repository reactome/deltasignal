# Tasks: Derived-edge recycling closures

- [x] T001 `_break_roles_env` + `DS_BREAK_ROLES` allow-list in `src/core/reaction_model.jl`; misspelt role throws.
- [x] T002 Two-pass component detection with closure marking per role; `closure_act` folded into `activator_break`; `depletion_break` field on `IndexedReaction`; `stats` keyword sink.
- [x] T003 `supply` read for closure depletions in `compute_reaction_output_vec`; `use_supply` covers the role list; singleton evaluation passes `supply = x` (`src/solvers/steady_state.jl`).
- [x] T004 Diagnostics (`scc_closures_*`, `scc_cyclic_before/after`, `scc_largest_after`) and additive API keys (`src/api/server.jl`, `docs/API.md`).
- [x] T005 `test/test_scc_break_roles.jl`: guard rail, unset bit-identical, welded fixture splits 1 → 2 with feed-forward assembly kept, depletion closure reads entry value and external depletion untouched, label/edge-order invariance (26 pass); full suite green.
- [ ] T006 Pre-registration P1–P5 in `research.md` §Results, committed before the arms.
- [ ] T007 Arms on `cat_os` from a pinned worktree at the committed SHA (logged per arm), vs `ab_fixed_point_ctrl.tsv` (current tree, cap 40,000): `assembly,depletion`; `catalyst`; `catalyst,assembly,depletion`; curator then experimental; probe prints TP53 `scc` census; no `src/` edits until the queue is empty.
- [ ] T008 Score P1–P5 (held-out/tuning with concentration, false change, per-pathway, experimental, AKT-KO cases, relabel churn on `cat_perm` for the winning list if any); commit; PR.
- [ ] T009 CI list + CLAUDE.md table (`test_scc_break_roles.jl`, 26).
