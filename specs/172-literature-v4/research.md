# Research: literature bundle v4 under curated pins (specs/172)

Run `lit172` at `caa1fd9e` (registered there, before any arm), every arm with `--bench DS_SOLVE_TIMEOUT=3600
--bench DS_PIN_CURATED_ROOT=<builds>/20261123-lit-v4-single-without`. Each arm log's `Protocol:` line prints the flag.
Results: `builds/<arm>/results/caa1fd9e/lit172`. Scores: `~/deltasignal-catalogs/analysis/172/score_172.json`, from
`score_172.py <builds> be1fbe08 caa1fd9e <out>`.

## Verdict: the gate FAILS (net −10). The strict selection is close to inert on the gated axes.

| check | result |
|---|---|
| FR-002 control parity (v4 control `lit172` vs v3 control `lit170`) | identical on every axis |
| gate: v4-strict-single vs control | **primary net −10** (korkut −10: 0 fixed / 10 broken; experimental 0: nothing moves) |
| shuffles | not run (net ≤ 0, as registered) |
| v4-single vs control (described) | primary net −554 (korkut −470, experimental −84), the same as specs/170 |
| ties: v4-single vs v3-single (described) | korkut and experimental identical; curator held-out −10 (11 / 21) |

The prediction was half right.
- **Right:** v4-single reproduced specs/170. The PIP3 collapse is gone in the strict arm, which has no literature edge
  into MTOR.
- **Wrong:** I expected the strict arm to lose on korkut in proportion to its size. Instead it reaches **no** new
  korkut or experimental readout (0 gained on either axis), so it has nothing there to over-call.

## What the strict arm does

- **korkut, 10 breaks:** all are one AKT readout in one pathway, under the SRC-inhibitor drug sets, moving 1.125 → 1.25
  across the change cutoff. Those drug sets did not engage their target in these cells (specs/138).
- **experimental:** no case moves.
- **Curator harm check** (readouts the curated network already reaches; reported, not gated):
  - held-out: 3 fixed / 38 broken;
  - tuning: 0 / 5.
  - 24 held-out cases lose coverage altogether (they read exactly 1.0).
- **Newly reached curator readouts** (uninformative under the 2026-11-20 rule): 11, change-call precision 0.45.

### Traced: a literature reaction that re-makes a curated complex masks knockouts

The largest harm class is 14 RNF8 knockout cases in DNA Double-Strand Break Repair, plus XRCC4 and XRCC6, 2 each.
Every readout reads exactly 1.0. Traced on trace servers at `caa1fd9e` with the bench's own pins:

- **The knockout is pinned at the same node in both arms.** No hold, pool or structural rule differs between the
  solves.
- **How the edges get in.** In the curated network the gene reaches its first complex only through AND `assembly`
  edges. The literature layer adds a one-input reaction whose output is that same curated complex. Its input is the
  complex's other component (Zn2+ for RNF8:Zn2+, LIG4 for XRCC4:LIG4, XRCC5 for XRCC5:XRCC6), never the knocked-out
  gene.
- **Why the knockout vanishes.** The complex now has an OR producer beside its assembly. Under the default
  `DS_OR_COMBINE=max` it reads max(assembly = 0, literature producer = baseline) = baseline. The knockout dies at the
  first hop: 1 node off baseline instead of 4,192.
- **Causal check.** Removing only those three literature output edges, in a copy of the pathway, restores every one of
  the 18 readouts exactly to its curated value.

This is a general hazard, not a quirk of these three edges. Any literature reaction that produces an existing curated
complex from a subset of its components masks every knockout of the missing components.

## What this means

1. Across specs/166, 168, 170 and 172, no literature selection has improved the primary pool against the curated
   network.
   - **Large bundles:** reach new readouts and over-call them. korkut change-call precision on newly reached readouts
     was 0.009 in specs/168 and 0.015 in specs/170.
   - **The strict selection:** reaches nothing new on the gated axes, so it cannot gain there.
2. Two structural hazards are now traced, each worth a check when a bundle is built:
   - an inhibitory edge against a curated activation, which tips a bistable component (specs/170, PIP3);
   - a reaction that re-makes a curated complex from part of its components, which masks knockouts through OR-`max`
     (this spec).
3. Nothing is adopted, and no default changes.
