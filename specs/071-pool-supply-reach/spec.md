# specs/071: a pool's supply reached only by regulation is still its supply

**Status:** pre-registered 2026-10-18; results in `research.md` (prediction 1 failed, 49 of 100).

Adam, 2026-10-18, on specs/070 §5: fix the TP53/MDM2 mechanism first ("ok").

## The defect

specs/039 review finding 6: a cycle-balance pool's supply excludes any producer
the pool's own forms reach, because the pool could feed it. The reach follows
every activator edge, catalysts and positive regulators included.

On the faithful build (`20261018-0815_ddbd828_dsub`, specs/070), MDM2
(R-HSA-69486) is a pool state. TP53 binds the MDM2 gene, and that complex
(R-HSA-3700998) is a positive regulator of MDM2 translation (R-HSA-3700992),
so the pool reaches its only supply and drops it. With no supply the pool
reads s = 1, and an MDM2 knockout or overexpression reaches nothing. 100
curator cases break this way (specs/070 §5). Walking material edges only, the
MDM2 pool does not reach translation (checked on the bundle).

No protein flows back along a regulator or catalyst edge. A producer the
pool reaches only that way is genuine feedback on the supply, not the pool
feeding itself.

## The rule

`DS_POOL_SUPPLY_REACH` = `all` (default, unchanged) | `material`. Under
`material`, the finding-6 supply test and the supply chain's ancestry walk
(which decides the dropped depletion edges, amendment 3) skip (source, target)
pairs whose every positive edge is a catalyst or regulator edge
(`pool_signal_edge_pairs`). A pair that is also an input carries material.
The amendment-5 re-evaluation set is unchanged. The solve reports
`pool_supply_reach` and `cycle_supply_regulated_kept`. Tests:
`test/test_cycles.jl` "supply reached only by regulation stays a supply".

## Arms (solver only; no rebuild)

- `psr`: canonical build `20261017-1246_5e663d9`, `--server
  DS_POOL_SUPPLY_REACH=material`, against canonical scoring `results/6fb4329`.
- `psrd`: faithful build `20261018-0815_ddbd828_dsub`, same override, against
  `dsub` (`results/40adace/dsub`) and against canonical.

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade`; both figures
are reported, on paired cases and on all cases valid in both runs.

**Gate (registered).** This is a defect fix. `psr` is adopted as the default
if, against canonical: curator held-out net ≥ 0, experimental net ≥ −5, and no
pathway below −10, RAF set aside. The usual gain gate (held-out net > 0 at
McNemar p < 0.05) is reported beside it. `psrd` informs the faithful-network
work and is not gated.

**Predictions:**
1. `psrd`: at least 80 of the 100 TP53 MDM2 cases broken by `dsub` return to
   their canonical class.
2. `psrd`: TP53 non-convergence does not fall below `dsub`'s 663 by this rule
   alone. MDM2 supply reintroduces the TP53–MDM2 feedback, and the
   non-convergence is the LM fallback skipping pool components (specs/070 §5).
3. `psr`: few cases move on the canonical build (no pathway beyond ±10). The
   canonical TP53 has no pools, so whether any canonical pool has a
   regulation-only supply is not known in advance.
