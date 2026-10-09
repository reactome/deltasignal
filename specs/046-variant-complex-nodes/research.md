# specs/046 research

## Implementation status (2026-10-09, logic-network-generator `feat/variant-nodes`)

| Commit | What |
|---|---|
| 0ff17fc | `src/variant_keys.py`: key, choice enumeration (D2: one variable per set per reaction; D3: homo), parser, leaves |
| (next) | `variant_parts`: one format-tolerant parser; every member-reading site routed through it |
| 60147f3 | `src/variant_emission.py` behind `LNG_VARIANT_NODES`; per-variant `containment.csv`; variant-aware boundary layer |
| ba672e2 | `nodes.csv` `member_leaves` = a variant's terminal leaves; emitter unit tests |

The emitter:
- output-only slots are bound by isoform-level reference identity, ties by
  sorted rank and counted, replacing F9's positional tie-break;
- negative regulators are expanded (D4);
- over the cap there is one copy with plain stable ids (the first step of D6).

**FGFR2 (R-HSA-5654738), against the canonical build:**

| | Control | Variant nodes |
|---|---|---|
| Required inputs beyond the curated ones | 1,156 | **0** |
| Cut hand-offs | 544 | **0** |

- All 1,782 variant input edges are fed, with **0** missed joins (I3b).
- The build is identical under `PYTHONHASHSEED` 0 and 7 (10,382 edges).
- Size: 9,912 curated and 470 boundary edges, against 4,880 and 156.

Tests: 366 pass, flag off.

**Not yet in:** the pooled cap fallback (D5/D6) and the heteromer arm (D3).

## Pre-registration: catalog arm (before it runs)

- **Build:** `--variant vn` from `feat/variant-nodes` (ba672e2),
  `--env LNG_VARIANT_NODES=1`.
- **Control:** the canonical build `20261004-1959_6990015`.
- **Scoring:** solver main (`flow` order), code defaults.

**Structural checks** (must hold before any score is read):
1. I3 cut hand-offs = 0 and I3b missed joins = 0, catalog-wide.
2. Coverage: valid cases and the `Pinned:` count within 1% of the control on
   both axes. A drop means genes no longer resolve to nodes; that is a
   generator or benchmark bug to fix, not a result.
3. The build completes for all 92 pathways.

**Predictions:**
- **RAF:** largely unchanged. Its large steps are capped either way, and this
  build uses the single-copy fallback.
- **The pathways with the most cut hand-offs** (FGFR2/3/4/1, R-HSA-194315,
  Mitotic G1, DDX58, EPH-ephrin) change most. Signal that was cut at a
  dissolved complex now propagates.
- **Single-isoform catalyst knockouts move toward no change.** Losing one
  catalyst member now stops only that member's copies. Previously a
  `product` pool stopped every copy.
- **Direction on held-out and experimental is not predicted.** That is the
  honest position: earlier faithful fixes went both ways.

**Decision rule:**
- adopt if held-out ≥ −15 and experimental ≥ −15, and checks 1–3 hold;
- a loss beyond the floor on either axis is traced before anything else;
- report both axes, held-out and tuning, the pathways and genes moved,
  McNemar p, and convergence.
