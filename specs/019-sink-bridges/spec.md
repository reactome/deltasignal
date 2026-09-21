# Feature Specification: Acyclic sink bridges

**Feature Branch**: `feat/019-sink-bridges` · **Created**: 2026-09-20 · **Status**: Draft (measurement pre-registered in research.md)

## Why

After the welded-loop fix (specs/018), the largest remaining error class on the
production catalog is **`no_path`: 1,526 of 3,686 wrong cases** — the gene and
the readout are both in the network but no route joins them. Two pathways hold
686 of them: Mitotic G1 (430 of 884 cases wrong) and Interferon α/β (256; 354
of 448 wrong). The curator oracle shows why: the curators' route continues
through a released protein, ours stops at the **dissociation sink** — the
freshly minted readout handle that boundary expansion creates for every
subunit of a terminal complex and that has no outgoing edge (Mitotic G1: 264 of
430 severed routes break at a sink). Bridging sinks was measured harmful twice
(2026-05, −15pp; the silo bridge, −73/−77 held-out) — on catalogs whose giant
components were welded by our own edges, and with no cycle guard.

## The rule

A dissociation sink gains an OR edge (`sink_bridge`) to each consuming copy of
the same entity **only where that copy cannot reach the sink**, so no cycle is
closed. `LNG_SINK_BRIDGES=1` (default off); `LNG_SINK_BRIDGE_MAX_FANOUT=n` skips
sinks with more than n acyclic consumers.

## User Scenarios

1. **P1** A modeller knocks out USP18 in Interferon α/β. Today the readouts are
   unreachable (NORMAL by default). With bridges, the released ISGF3 subunits
   reach their consumers and the readouts respond.
2. **P1** Nothing that was acyclic becomes cyclic: the cycle census is unchanged.
3. **P2** The rule is measured against a same-catalog control (`DS_SKIP_EDGE_TYPES=sink_bridge`), zero churn, both axes, with and without the fan-out cap.

## Requirements

- FR-001 Off by default; `any`-style escape not needed (the flag is the escape).
- FR-002 A bridge is emitted only if the consumer cannot reach the sink over every edge emitted so far.
- FR-003 Bridges are `pos`, `or`, `edge_type = sink_bridge`, joined on the base stable id (variants included).
- FR-004 The emitter logs bridges, skipped cycle-closers, fan-out-capped sinks and fan-out statistics.
- FR-005 The measurement is pre-registered and committed before regeneration.

## Success Criteria

- SC-001 Cyclic components identical before and after (210 → 210); oracle severance falls catalog-wide.
- SC-002 Held-out ≥ +100 vs the same-catalog control (p < 0.05) with the concentration stated; Interferon α/β and Mitotic G1 move as predicted; false change rise smaller than the gain; experimental conditioned ≥ −10.
