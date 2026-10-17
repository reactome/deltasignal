# specs/067: specialization bridges — a specific complex is one member combination of a generic one

**Status:** pre-registered 2026-10-16, before any build or arm. Measured 2026-10-17 (research.md): `add` passes the gate (held-out +59, p 8.7e-7, experimental 0, worst pathway −7); `replace` fails it (DNA Damage Bypass −14). **Adopted as the generator default 2026-10-17 (Adam)**: `LNG_SPECIALIZATION_BRIDGE=add` (LNG #117); canonical rebuild follows. From a census of
the held-out misses on canonical `results/0709a1e`: Signaling by PDGF loses 80
held-out cases as `no_path`, 60 of them the ligand and its receptors
(PDGFB, PDGFRA, PDGFRB, both directions, every readout). specs/023 counted that
cluster; nothing has traced it until now.

## Diagnosis (traced 2026-10-16, canonical build `20261016-0006_da08370`)

Probes are in `~/deltasignal-catalogs/probe_spec/` (census, reach estimates)
and `~/deltasignal-catalogs/probe_osc/pdgf*.cypher` (Neo4j, Release 97).

1. **PDGF, verified in Neo4j.** Reactome curates two hand-offs between a
   GENERIC complex (components are sets) and a SPECIFIC complex (components
   are members of those sets), with no shared entity and only a
   `precedingEvent` linking the reactions:
   - "PDGF dimer binds two receptors simultaneously" (R-HSA-186773) outputs
     PDGF:PDGF receptor dimer (R-HSA-186766 = Active PDGF dimers set :
     PDGF receptor dimer set). It is consumed by nothing. The three
     autophosphorylation reactions (R-HSA-186786, 389086, 389083) list
     186773 as `precedingEvent` and consume specific complexes
     (R-HSA-389078 = PDGF-BB : beta homodimer, 389077, 389079).
   - Autophosphorylation outputs specific phospho complexes (R-HSA-389074,
     389076, 389073). Every downstream step (GRB2:SOS1, PLC-gamma, PI3K, SRC,
     STAT, SHP2, Crk, Nck, Grb7: about ten reactions) consumes the generic
     PDGF:Phospho-PDGF receptor dimer (R-HSA-186811 = Phospho-PDGF receptor
     dimer set : Active PDGF dimers set), whose only producer is a later
     dissociation. R-HSA-389074 is exactly one member combination of
     R-HSA-186811 (186805 with 184206).
   - In our network the first gap is bridged by accident (the specific
     complexes are roots, decomposed by boundary assembly from their
     components); the second is not. A PDGFRB perturbation reaches 314 of
     2,434 nodes and none of the signalling readouts.
2. **The relation.** Complex S *specializes* generic complex G when their
   components (with stoichiometry) match one to one and each S component is
   the G component, a member (transitively) of it, a set whose members are a
   subset of its members, or a complex that itself specializes it. Then S
   is one of the species G names; the generator's variant nodes (specs/046)
   already give G one node per member combination, so a copy of S and the
   copy of G with the same `member_leaves` are the same species under two
   entity ids.
3. **Census** (Release 97, the 92 catalog pathways, `spec_census.py`): 117
   produced/consumed pairs in 31 pathways where one complex specializes the
   other. 61 are linked by `precedingEvent` (a reaction consuming the target
   lists a reaction producing the source), in 24 pathways. Every linked pair
   has an exact copy match (same `member_leaves`) in the canonical build.
   Spot-checked by name: GLI:SUFU → GLI1/2/3:SUFU (Hedgehog), NODAL
   ligand:receptor → the ACVR1B / ACVR1C forms (NODAL), CD28:B7-1 →
   CD28:B7-1,B7-2 (CD28 family), MonoUb-PCNA on damaged template → the
   lesion-specific templates (DNA Damage Bypass).
4. **Reach estimate** (`reach_replace.py`, reachability only, no solve):
   bridging the 61 linked pairs copy to copy routes 116 held-out `no_path`
   cases (58 expected DOWN, 58 UP) in five pathways (PDGF 82 cases including
   2 NORMAL, DNA Damage Bypass 28, NODAL 16, Hedgehog 16, CD28 12), and gives
   a route to 38 held-out cases now right as NORMAL. No tuning case changes
   reach. Bridging every pair, linked or not, routes the same 116 held-out
   cases and adds 2 more held-out NORMAL cases and 32 tuning cases in WNT
   (24 of them now right as NORMAL); matching copies by entity instead of
   by identical leaves adds a further 50 held-out NORMAL cases in TFAP2. The
   rule therefore requires both the link and the exact copy match.

## Relation to earlier negative results

- `LNG_HANDOFF_EDGES` (July; three variants, net negative each time) bridged
  `precedingEvent` gaps on a shared COMPONENT, a subunit handed from one
  complex to another, which coupled unrelated readouts through promiscuous
  subunits. This rule bridges only an identity: the same species under two
  ids, copy to copy with identical leaves. No shared-subunit test is used.
- `LNG_DIAGRAM_SET_MEMBER` (default off, never benchmarked) carries the
  diagram's `EntitySetAndMemberLink` between a specific and a generic complex.
  Its recorded objection: the target also has AND assembly in-edges, and the
  solver combines AND and OR inputs as max(and, or), so a decrease through the
  link is discarded. Measured here: 189 of the 203 matched copy pairs land on
  a copy with assembly in-edges, which the boundary decomposition added
  because the target had no producer. Variant `replace` below answers that
  objection; `add` measures it.
- The composition-edge and silo-bridge families (specs/016, 018, 019)
  welded cycles. This rule adds edges only in the curated direction (the
  `precedingEvent` order), never both ways; mutual pairs that are not linked
  (CD28 9909625 / 9909637) are excluded by the link test.

## Rule (generator, `LNG_SPECIALIZATION_BRIDGE`, default `off`)

- `off`: no change.
- `add`: for every pair (A produced, B consumed, in the pathway) where A
  specializes B or B specializes A, and some reaction consuming B lists some
  reaction producing A as `precedingEvent`, emit an edge from each copy of A
  to each copy of B with the same `member_leaves`: `pos`, `or`, edge type
  `specialization`, written to `logic_network.csv` (curated: the relation and
  the order are both in Neo4j).
- `replace`: as `add`, and a bridged target copy loses its boundary
  `assembly` / `composition` in-edges (the decomposition stood in for a
  producer the pathway lacked; now the curated one feeds it).

Solver: no change. A positive `or` edge is an OR activator already.

## Arms and measurement

Two catalog builds from the generator branch, everything else at defaults:
`LNG_SPECIALIZATION_BRIDGE=add` (`sbadd`) and `=replace` (`sbrep`). Each is
benchmarked with `scripts/catalog.sh bench`-equivalent settings on solver
main (`f0aa666` or later with identical `src/`/`bench/`) and compared with
canonical `results/0709a1e` (`holdout_report.py --compare`, with and without
`--set-aside RAF_MAP_kinase_cascade`). Report: held-out, tuning,
experimental net with McNemar p; per-pathway net; distinct readouts; each
held-out class change in the five predicted pathways, and every change
elsewhere traced; content diff (edges added / removed per pathway; nodes in
non-trivial strongly connected components per pathway, before and after).

## Predictions (from the reach estimate)

- Class changes concentrate in PDGF, DNA Damage Bypass, NODAL, Hedgehog and
  CD28. Other bridged pathways change values but few classes.
- Same-sign genes newly routed (PDGFB, PDGFRA, PDGFRB, RAD18, REV3L, PTCH:
  70 cases): at least 50 fixed.
- Inverse-sign genes newly routed (PTPN12, CER1, LEFTY2, SMO, CTLA4,
  PPP2R1A: 46 cases): at most 15 fixed; these need an inhibitory step the
  route may not carry.
- Of the 38 NORMAL cases newly routed: at most 25 broken.
- Held-out net between +20 and +60 for `replace`; `add` no better than
  `replace`.
- Tuning: no class change. Experimental: none outside RAF.
- No pathway gains more than 10% cyclic nodes.

## Decision rule

A route-adding structural change, so the specs/046 gate, with RAF set aside:
- **Adopt** a variant if curator held-out net > 0 with McNemar p < 0.05,
  experimental net ≥ −5, and no pathway below −10.
- If both pass, adopt `replace`, unless `add` is better on held-out by 5 or
  more.
- Otherwise do not adopt; record which pathways carry the loss and trace
  them.
- Either way report both readings (with and without RAF).
