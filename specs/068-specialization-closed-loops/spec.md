# specs/068: a specialization bridge is identity, not outside supply, to a carrier loop

**Status:** pre-registered 2026-10-17, before any build or arm. Measured 2026-10-17 (research.md): held-out +18 (18 / 0, p 7.6e-6), experimental 0, every prediction held; passes the gate. **Adopted as the generator default 2026-10-17 (Adam)**: `LNG_SPECIALIZATION_CLOSED=twin` (LNG #118); canonical rebuild follows. Follow-up to
specs/067's largest traced loss.

## Diagnosis (traced on canonical build `20261017-0654_d3ebf38`)

1. specs/067 (`LNG_SPECIALIZATION_BRIDGE=add`) broke 18 held-out PDGF cases:
   PLCG1 knockout and overexpression, 9 readouts each, curator NORMAL. The
   generic PDGF:Phospho-PDGF receptor dimer (R-HSA-186811) has one curated
   producer in Release 97, "Activated PLC gamma dissociates from the PDGF
   receptor", which releases the receptor unchanged after PLC-gamma binds:
   a carrier loop (specs/049). Before 067 nothing outside the loop fed it, so
   the generator marked it `closed` and `closed_twin` held it to its gene;
   PLCG1 could not move it.
2. The generator's closedness test (`find_pools`, `into`) counts every
   positive entity in-edge as supply, so the new `specialization` edge from
   the specific phospho-receptor (R-HSA-389074/76/73) makes the loop open.
   The release edge is then an ordinary OR input, and PLCG1 sets the
   receptor's level (0.67× under knockout, 34× under overexpression).
3. Census, `carrier_loops.csv` of the 067 build against the 066 build
   (`20261016-0006_da08370`): 1,087 listed free forms in both; exactly 10 change,
   all copies of R-HSA-186811 in Signaling by PDGF, `closed` True → False.
   No other loop in the catalog is opened by a bridge.

## Why identity

A specialization edge joins two copies of one species under two entity ids
(specs/067). It is the same relation `closed_twin` already uses: a closed
free form reads its own gene's root copy (`carrier_supply`), which also
does not make the loop open. The specific receptor's autophosphorylation
does supply receptor, and the edge carries it; what it must not do is turn
PLC-gamma's bind-and-release into the receptor's production.

## Rule (generator, `LNG_SPECIALIZATION_CLOSED`, default `off`)

- `off`: no change (closedness counts `specialization` sources as supply).
- `twin`: the carrier-loop closedness test ignores in-edges of type
  `specialization`. The edges stay in `logic_network.csv`; only
  `carrier_loops.csv` changes. The solver is unchanged: under `closed_twin`
  a closed free form loses its listed loop in-edges and keeps every other
  in-edge, so the specialization edge still feeds it.

## Arm and measurement

One catalog build from the generator branch with `LNG_SPECIALIZATION_CLOSED=twin`,
everything else at defaults (arm `sctwin`), benchmarked on solver main
(`a177c6f` or later with identical `src/`/`bench/`) and compared case by case
with canonical `results/fa2c451` (`probe_spec/cmp.py`), with and without
`--set-aside RAF_MAP_kinase_cascade`. Report held-out, tuning, experimental
net with McNemar p, per-pathway net, and every class change traced.

## Predictions

- Content: `logic_network.csv` identical in every pathway; `carrier_loops.csv`
  differs only in the 10 PDGF rows above, all back to `closed` True.
- Class changes only in Signaling by PDGF.
- PLCG1: at least 16 of the 18 broken held-out cases fixed.
- The 80 PDGF cases specs/067 fixed: at least 75 still right.
- Held-out net between +10 and +20; tuning 0; experimental 0.

## Decision rule

The specs/067 gate, RAF set aside: adopt if curator held-out net > 0 with
McNemar p < 0.05, experimental net ≥ −5, and no pathway below −10.
Otherwise do not adopt, and trace the losses. Report both readings.
