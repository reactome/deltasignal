# specs/069: translesion synthesis re-enters its own start (diagnosis, no rule)

**Status:** diagnosis only, 2026-10-17. No rule pre-registered or adopted: every
structural definition tried was shaped on one pathway (below), and the curator
held-out split is SELECT (specs/051). Follow-up to specs/067's DNA Damage
Bypass losses.

Probes: `~/deltasignal-catalogs/probe_tls/` (one-pathway probe harness
`probe.sh` over a copy of the canonical pathway, with a probe-only copy of
`bench/` whose single-pathway guards print instead of exit and which dumps
each solve; `walk.py` follows the most-perturbed input back from a readout;
`mkcut.py` builds cut variants; `census*.py`). Canonical build
`20261017-1246_5e663d9`, solver `6fb4329`; the probe control reproduces the
canonical classes on all 126 DNA Damage Bypass cases.

## The route (traced by value, POLI overexpression)

POLI (80×) → boundary assembly into the root complex DNA polymerase
Y:MonoUb-PCNA:TLS template (R-HSA-5653660, no curated producer) → TRIM25,
USP10 (releases the polymerase; the backbone keeps the fold, mean over the
polymerase variants ≈ 26×), USP43 → POLD,POLE complete replication after TLS
(R-HSA-5653840) → **its outputs** → PCNA loading on damaged DNA (5651992) →
RAD18 monoubiquitination → generic MonoUb-PCNA template (5652003) →
specs/067 bridge → lesion-specific templates → every polymerase branch,
railing at 100×. With specs/067's spec→gen bridges into 5653660 this is a
closed positive loop.

Release 97 curates 5653840 with two outputs that feed the pathway's first
reaction: PCNA:POLD,POLE:RPA:RFC (5651799, consumed by 5651992) and Damaged
dsDNA (5652135), which a black-box event (5652137) turns back into Damaged
DNA Template (5651987, consumed by 5651992). The lesion is bypassed, not
repaired, so the cycle is real biology. MP-BioPath's 2019 network has
5653840 with **no out-edges** (export, not a hand edit; QA/changes has no
row for these ids); whether 2019 curation lacked the outputs or the exporter
dropped them is not established.

The curator answers treat 5653840 as terminal: every polymerase reaches it
(REV3L, RCHY1 KO → 5653840 UP are curator UP), but no polymerase moves the
generic template or a sibling branch (all NORMAL).

## Probes (DNA Damage Bypass only; class changes against the control)

| cut | edges dropped | fixed | broken |
|---|---|---|---|
| PCNA:POLD,POLE → loading (5651799 into 5651992) | 14 | 0 | 0 |
| black box → template (5652137 outputs) | 49 | 10 | 0 |
| both | 63 | **45** | **0** |
| every output of 5653840 | 294 | 45 | 0 |
| specs/040 rule A (`DS_SELF_FED_MODE=entry`) | 2 held | 0 | 0 |

45 fixed: 15 are the cases specs/067 broke, 30 were wrong before it (POLH,
POLI, POLK, USP43 at the generic template and sibling readouts). Pathway
accuracy 57.94% → 93.65% in the probe.

## Rule definitions tried (census over the 92-pathway build)

- "a copy whose producers are all downstream of itself, entity with no root
  copy": 4,260 nodes, 7,116 edges, 66 pathways. Every node of a pure cycle
  qualifies, so it would cut every loop edge, not the re-entry. Rejected
  without a run.
- the same, plus "the producer releases it (its leaves are a subset of an
  input's)": 248 nodes, 733 edges, 22 pathways. In DNA Damage Bypass it takes
  the TLS backbone intermediates (5653772, 5653782), which carry the curator-UP
  REV3L/RCHY1 → 5653840 route, and misses the template re-entry (a
  relabelling, not a release). Rejected without a run.

## Open

A label-free definition of a cycle's re-entry point that takes exactly these
two edges here is not found. Options: a curated-flag question to Reactome
(is TLS completion meant to feed new loading?), which needs Adam's approval
to raise; or a rule derived blind on other pathways and only then scored here.

## Rule T as registered (spec.md): refuted, 2026-10-17

Arm `t069`: `scripts/run_arm.sh t069 --catalog 20261017-1246_5e663d9_t069`
at solver `5e2320e` (results under that build, `results/5e2320e/t069`),
against `20261017-1246_5e663d9` / `results/6fb4329`.

| reading | fixed | broken | net | p |
|---|---|---|---|---|
| curator held-out, DNA Damage Bypass / DNA Repair / RAF set aside (**gate**) | 36 | 298 | **−262** | < 1e-4 |
| curator held-out, all | 36 | 298 | −262 | < 1e-4 |
| curator tuning | 13 | 320 | −307 | < 1e-4 |
| experimental | 1 | 104 | **−103** | < 1e-4 |

Coverage: 86 curator cases valid in the baseline only. Per pathway
(`bench/analysis/069/t069_perpathway.txt`): TP53 −274 (13 / 287), DSB
Repair −93, PTK6 −91, Base Excision Repair −43, HDR −32, TGF-beta −14,
KIT −8, PDGF −8, others ≥ −3; Fanconi +6. Experimental: TP53 −97, HDR −6.

Predictions: 1 failed (DNA Damage Bypass 73 → 73, **no case moved**);
2 failed (TP53, DSB, PTK6 moved by far more than 5); 3 had no sign. The gate
fails on every term.

**Why (traced).** The registered operation drops E's `input` edges into the
earlier reaction S. When E is also the pathway's own supply of S (the usual
case: the factor that starts the process is the one released at the end),
that severs the start, not the re-entry. TP53: the TP53 tetramer
(R-HSA-3209194), released by a terminal step, is cut as an input of five
early reactions, so TP53, MDM2, MDM4 and AKT perturbations no longer reach
their readouts. DNA Damage Bypass: cutting PCNA:POLD,POLE:RPA:RFC into
5651992 is the probe's "recycle only" variant, which moved nothing; and
Damaged dsDNA (5652135) is a set whose network nodes carry its MEMBERS'
ids, so matching the set id cut no black-box edge at all. The probe that
fixed 45 cut the terminal reaction's OUTPUT edges. Amendment 1 registers
that operation.

## Amendment 1 (cut the terminal output): fails the gate, 2026-10-17

Arm `t069o`: `scripts/run_arm.sh t069o --catalog 20261017-1246_5e663d9_t069o`
at `b14350a` (results under that build, `results/b14350a/t069o`), against
`results/6fb4329`.

| reading | fixed | broken | net | p |
|---|---|---|---|---|
| curator held-out, DNA Damage Bypass / DNA Repair / RAF set aside (**gate**) | 23 | 11 | **+12** | **0.058** |
| curator held-out, all | 68 | 11 | +57 | < 1e-4 |
| curator tuning | 30 | 8 | +22 | 0.0005 |
| experimental, RAF set aside (**gate**) | 2 | 6 | −4 | 0.29 |
| experimental, all | 2 | 6 | −4 | 0.29 |

holdout_report pairs only cases whose pin set is unchanged; the cut changes
which nodes are roots, so 2,174 curator and 167 experimental cases are
unpaired (DSB Repair 1,044, TP53 408, RUNX1 168). Counting every case valid
in both (`bench/analysis/069/t069o_perpathway.txt`): curator +83, the same
held-out +12 outside DNA Damage Bypass and RAF (none of the 1,044 DSB cases
changes class); experimental −2, and **−7 with RAF set aside** (TP53 −5,
ERBB2 −2), which would also fail the −5 bound.

Per pathway (curator): DNA Damage Bypass +45 (45 / 0, prediction 1 held:
the probe's cases, as expected, so not evidence), RAF +12 (all NF1 knockout
and overexpression; experimental NF1 KO +5), EGFR +12, ERBB2 +8, Fanconi +6,
TP53 +6, Insulin receptor +2, DNA Damage Reversal −1, **PTK6 −7**.

**Gate: fails** (p 0.058 ≥ 0.05). Not adopted.

Traced loss, PTK6 −7 (ERBB2/3/4 knockout and overexpression at
p-Y565-HIF1A, which reads exactly 1.0 in the arm): "PTPN1 dephosphorylates
PTK6" (R-HSA-8849435) is the pathway's last step and returns free PTK6
(R-HSA-8848000) to "Activated ERBB2 binds PTK6". The ERBB signal reaches
the free PTK6 that phosphorylates HIF1A only through that return, so the cut
disconnects it. That is a curated signalling cycle (phosphorylation and its
reversal), the same kind as the RHO GTPase GDP returns the census already
showed. `precedingEvent` order alone does not separate a signalling cycle
from the TLS re-entry; the class is impure.

## Where this leaves specs/069

- DNA Damage Bypass: cutting the terminal output fixes 45 and breaks none,
  now reproduced on the full catalog.
- A rule derived away from that pathway is not supported: its first form was
  harmful (−262), the corrected form gains +12 outside it at p 0.058, and
  its one traced loss is a cycle the rule should not cut. A third definition
  shaped on these results would be fitted to the SELECT split.
- Not established: whether Release 97 means TLS completion to feed new
  loading. That is a curation question (MP-BioPath's 2019 export had no
  out-edges on 5653840).

## What the diagram says (Adam: diagrams show intent; dashed lines are not in Neo4j), 2026-10-17

Release 97 layout `~/reactome-diagrams/97/R-HSA-73893.json`:

- The loop is DRAWN. Damaged dsDNA is one glyph (id 551): an output of
  5653840 and the input of the black box 5652137. PCNA:POLD,POLE:RPA:RFC is
  one glyph (id 71): an output of 5653840 and an input of loading 5651992.
  The two steps sit side by side on the top row (x 543 and 234, y ≈ 117),
  wired end to start.
- `precedingEvent` does not draw it: 5652137 has no predecessor (a start;
  followed by 5651992 and USP1 autocleavage 5655483), and 5653840 has no
  follower (an end).
- A dashed Interaction (link 552, "Activate") runs from Damaged dsDNA to the
  encapsulated pathway node Base Excision Repair (R-HSA-73884). It is not in
  Neo4j, and no BER reaction takes Damaged dsDNA as input: it marks where
  the bypassed lesion is repaired, not a hand-off the networks can carry.
  Interaction links are not read by the generator (LNG #41, deferred).

Reading: Reactome draws both continuations of bypassed DNA, still damaged
— another round of bypass (the drawn loop) and repair by BER (the dashed
exit). The loop is intended curation, not an artefact. The 45 cases are
then a disagreement between Reactome's model and the 2019 expert answers
(which, like MP-BioPath's 2019 network with no out-edges on 5653840, treat
completion as the end), i.e. a ground-truth question rather than a network
error.

Generalised over rule T's 291 pairs (`bench/analysis/069/diagram_check_r97.tsv`:
does a diagram holding both reactions draw L's output and S's input as one
glyph of E or a member?): **196 shared glyph, 82 separate glyphs, 13 on no
common diagram.** Both DNA Damage Bypass re-entries are shared, and so are
the cycles rule T should not cut (PTK6's PTPN1 return of free PTK6, every
RHO GTPase GDP return, TP53's tetramer returns). A shared glyph marks a
cycle the curators drew; it does not separate the TLS loop from signalling
cycles. The 82 separate-glyph pairs (mostly Base Excision Repair: released
enzymes drawn as fresh copies) are the curator-drawn NON-loops; they do not
include the TLS re-entries, so they are not this pathway's lever.
