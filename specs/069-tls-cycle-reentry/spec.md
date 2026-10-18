# specs/069: terminal re-entry (rule T)

**Status:** pre-registered 2026-10-17, after the variant catalog was built
(structure only) and before any arm on it. Follows the diagnosis in
`research.md`. Adam, 2026-10-17: derive the rule away from DNA Damage Bypass
and test it before any question goes to Reactome curation.

## Why a rule, and why this one

DNA Damage Bypass's last reaction, "POLD,POLE complete replication after TLS"
(R-HSA-5653840), feeds the pathway's first reaction, so any polymerase
perturbation circulates and rails at 100× (`research.md`). Cutting the two
re-entries fixed 45 cases and broke none in a one-pathway probe. A cut chosen
on those 45 cases would be fitted to them, so the rule below is defined by
Reactome's own event order, and its evidence is read **without** DNA Damage
Bypass and its parent DNA Repair, whose outcome is already known.

Reactome orders a pathway's events with `precedingEvent`. A reaction that
nothing in its pathway follows is where the curated sequence ends. If such a
reaction outputs an entity that an EARLIER reaction (one preceding it,
transitively) consumes, the network gains an end-to-start cycle that the
curated order does not draw.

## Rule T

Within each catalog pathway (its `hasEvent*` reactions; `precedingEvent`
restricted to them), for every reaction L with no follower, every output E
of L, and every reaction S that precedes L (and that L does not precede) and
takes E as an input: drop the network's `input` edges from copies of E into
S's reaction nodes. Exempt, because existing rules already handle them:
- the generator's carrier-loop free forms (`carrier_loops.csv`; closed by
  `DS_CARRIER_LOOP_MODE`, specs/049);
- cofactors (`cofactors.csv`; held inert, specs/007).

Census on Release 97 (`bench/analysis/069/terminal_reentry.py census`,
input `reactions_r97.tsv` from `reactions.cypher`): **291 (entity, earlier
reaction) pairs in 45 pathways**; applied to canonical build
`20261017-1246_5e663d9` it drops **2,249 edges in 43 pathways**
(variant build `20261017-1246_5e663d9_t069`; per-pathway counts in its
`BUILD.json`). The pairs fall in three groups seen in pathways other than DNA
Damage Bypass:
- repaired dsDNA fed to the black-box "damage" event that starts a repair
  pathway (Base Excision Repair, R-HSA-73942, DNA Repair);
- enzymes and machinery released by the last step and reused by an earlier
  one (APEX1, TFIIH, ALKBH2/3, BRCA1, ubiquitin);
- GTPase GDP forms returned by the GAP step to the GEF step (RHO GTPase
  pathways, 194315 and 9012999).

In DNA Damage Bypass it cuts 3 pairs: the two diagnosed re-entries
(5651799 into 5651992; Damaged dsDNA into black box 5652137) and ubiquitin
released by USP1 (5655466) into RAD18's PCNA monoubiquitination (5652009),
which no probe tested.

This is a solver-input transform measured as a variant catalog. If adopted,
it moves into the generator behind a flag before becoming a default.

## Arm and gate

Arm `t069`: `scripts/run_arm.sh t069 --catalog 20261017-1246_5e663d9_t069`,
both axes, against the canonical scoring (`results/6fb4329`).

**Primary gate (blind), with DNA_Damage_Bypass, DNA_Repair and
RAF_MAP_kinase_cascade set aside:** curator held-out net > 0 with McNemar
p < 0.05, experimental net ≥ −5, and no pathway below −10. Otherwise not
adopted, and the losses are traced.

Reported as well, not gating: every-pathway figures with RAF set aside and
with RAF kept, and DNA Damage Bypass on its own.

## Predictions (written before the arm)

1. DNA Damage Bypass gains at least 35 cases (the probe's 45 less a margin
   for the untested ubiquitin cut), with at most 3 broken.
2. Outside DNA Damage Bypass and DNA Repair, no sign is predicted: the gate
   decides. Most of the 43 pathways should move by fewer than 5 cases, because
   the cut edges run from a pathway's last step, and most readouts sit
   upstream of it.
3. Experimental: no sign is predicted. The tuning pathways are exposed:
   PIP3/AKT loses 386 edges (PTEN and p-MAPK dimers returned to earlier
   steps), WNT 109, HDR 80, TP53 62, G2/M 15, Cell Cycle Checkpoints 13,
   RAF 4, S Phase 2, Prophase 1, ERBB2 0. The gate's −5 bound decides.
