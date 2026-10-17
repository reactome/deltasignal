# specs/067 research: arms `sbadd` and `sbrep`

Generator: LNG branch `specialization-bridges` at `0f99244` (LNG PR #117),
default off. Builds, all 92 pathways complete:
- `20261016-2333_0f99244_sbadd` (`LNG_SPECIALIZATION_BRIDGE=add`): 500,965
  edges (canonical 500,762, so +203, exactly the census copy edges);
- `20261016-2334_0f99244_sbrep` (`=replace`): 500,589 edges, 19 nodes fewer
  (component nodes that only fed a replaced decomposition).

The module reproduces the census against live Neo4j: 61 linked pairs and
203 copy edges in 24 pathways, none elsewhere. The content diff
(`probe_065a1/content_cmp.py`) shows exactly those 24 pathways differing.

Arms: `run_arm.sh sbadd --catalog …_sbadd` and `sbrep --catalog …_sbrep`,
solver `0e410b5` (`src/` and `bench/` identical to main `f0aa666`), against
canonical `results/0709a1e`. Results:
`builds/<build>/results/0e410b5/<arm>/`. Case comparison:
`~/deltasignal-catalogs/probe_spec/cmp.py`.

## Results

| arm | held-out net | fixed | broken | McNemar p | tuning | experimental | worst pathway |
|---|---|---|---|---|---|---|---|
| `sbadd` | **+59** | 101 | 42 | 8.7e-7 | 0 | 0 | DNA Damage Bypass −7 |
| `sbrep` | +56 | 122 | 66 | 5.4e-5 | 0 | 0 | DNA Damage Bypass −14 |

RAF has no moved case on either axis, so the set-aside reading is the same.

Curator held-out in release: 89.50% / 0.8646 → **89.81% / 0.8695** (`sbadd`),
89.80% / 0.8695 (`sbrep`). All pathways, every case: 86.70% → 86.95%
(`sbadd`).

Per pathway (`sbadd` / `sbrep`): PDGF +61 / +60, NODAL +8 / +16, Hedgehog
−1 / −2, CD28 family −2 / −4, DNA Damage Bypass −7 / −14. Nothing else moves.

## What moved (`sbadd`)

- **PDGF, fixed 80:** PDGFB, PDGFRA and PDGFRB, both directions, 10 readouts
  each (60), and PTPN12 both directions (20). The receptor now reaches the
  signalling steps through the specific → generic phospho-receptor bridge.
- **PDGF, broken 19:** PLCG1, both directions, 9 readouts each (18), and one
  PDGFRA case. Traced: the generic PDGF:Phospho-PDGF receptor dimer
  (R-HSA-186811) had one producer, "Activated PLC gamma dissociates from the
  PDGF receptor", which releases the receptor unchanged after PLC-gamma
  binds: a carrier loop. With no producer outside the loop it was a closed
  carrier free form, held to its root copy by `closed_twin` (specs/049), so
  PLCG1 could not move it. The bridge gives it an outside producer, so the
  loop's release edge becomes an ordinary OR input and PLCG1 now sets the
  receptor's level (0.67× under knockout, 34× under overexpression; curator
  NORMAL). The bridge is right; the carrier rule stops applying.
- **NODAL, fixed 8:** CER1 and LEFTY2 knockouts.
- **DNA Damage Bypass, fixed 8, broken 15:** RAD18 overexpression fixed at
  six readouts (the PCNA ubiquitin ligase now reaches the lesion-specific
  templates). Overexpressing any translesion polymerase (POLI, POLK, REV1,
  REV3L, USP43) or knocking out RCHY1 now drives the POLH:MonoUb-PCNA
  readouts (5654981, 5654988) to 100× (curator NORMAL): the generic
  MonoUb-PCNA:damaged template feeds every lesion-specific template, so
  one polymerase reaches its siblings' complexes (the over-coupling family).
- **Hedgehog −1, CD28 −2:** PTCH fixed 3; SHH and SMO broken 4; CTLA4 and
  PPP2R1A 2 fixed, 4 broken.

`sbrep` adds both directions where `sbadd` fixed one (CER1, LEFTY2, PTCH: the
max(and, or) objection in the spec, confirmed) and breaks the same
DNA Damage Bypass and CD28 cases in both directions too.

## Predictions, scored

| prediction | outcome |
|---|---|
| class changes concentrate in PDGF, DNA Damage Bypass, NODAL, Hedgehog, CD28 | **held**: exactly these five, nothing else |
| same-sign genes: at least 50 of 70 fixed | **held**: 70 (PDGF 60, RAD18 6, PTCH 3, REV3L 1) |
| inverse-sign genes: at most 15 of 46 fixed | **wrong**: 30 (PTPN12 20, CER1 4, LEFTY2 4, CTLA4 1, PPP2R1A 1) |
| at most 25 of the 38 newly routed NORMAL cases broken | **held**: 19. But 23 of the 42 losses came from outside that set, through routes that already existed (PLCG1 18 through the carrier rule, RCHY1 4, SMO 1), and were not foreseen |
| held-out net +20 to +60 (`replace`) | **held**: +56 (`add` +59) |
| `add` no better than `replace` | **wrong** by 3 |
| tuning: no class change; experimental: none outside RAF | **held**: 0 and 0 |
| no pathway gains more than 10% cyclic nodes | **wrong**: Hedgehog 76 → 136 (+79%), TGF-beta 158 → 237 (+50%), DNA Damage Bypass 580 → 642 (+11%), Mitotic G2 142 → 150 (+6%), the same in both arms |

The cyclic-node prediction failed: an identity bridge into a generic form
that a later step recycles closes a loop (GLI:SUFU in Hedgehog). It is not
part of the decision rule, but specs/018 is why it was predicted, and the
Hedgehog and DNA Damage Bypass losses sit in the grown components.

## Decision (by the registered rule)

- `sbadd` **passes**: held-out +59 at p 8.7e-7, experimental 0 (≥ −5),
  worst pathway −7 (≥ −10).
- `sbrep` **fails** on the worst-pathway floor (DNA Damage Bypass −14).
- Only one variant passes, so by the rule the candidate is
  `LNG_SPECIALIZATION_BRIDGE=add`. **Adopted as the generator default
  2026-10-17 (Adam)**; the losses above are follow-ups (the PLCG1 carrier
  loop opened by an outside producer; TLS polymerases coupled through the
  generic template; the grown Hedgehog and TGF-beta components).
