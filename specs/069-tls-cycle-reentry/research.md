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
