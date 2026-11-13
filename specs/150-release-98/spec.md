# Feature Specification: The Release-98 Upgrade, as a Measurement

**Feature Branch**: `specs/150-release-98`

**Created**: 2026-11-13

**Status**: Draft (pre-registered before any 98 measurement)

**Input**: Adam approved upgrading to Release 98 "when you think would be ideal" (2026-11-07) and supplied the mechanism (the graphdb image's `Release98` tag). The ideal is now: the leads arc is shipped, and the record has been waiting for this since specs/061 ("Release 98 remains the replication").

## Why this is a measurement, not maintenance

Every structural rule adopted since specs/061 was validated on Release-97 networks. Release 98 rebuilds every network from a year of new curation the rules have never seen. Scoring the FROZEN rules on the 98-built benchmark — before any adjustment — is the replication the record reserved: if the held-out figure holds, the rulebook generalizes beyond the curation it was developed on; if it falls, that is the most important fact available and nothing else proceeds until it is traced.

## Pre-registered replication protocol (R1, run FIRST)

- Networks: the benchmark pathway list built on Release 98 by the **frozen generator tag** (`freeze-2026-10-26`) with the 98 diagram mirror — so the ONLY change from the frozen record is the release.
- Rules: the frozen solver defaults, untouched; the standard protocol (`run_arm.sh`, curator + experimental axes).
- Ground truth: unchanged (the MP-BioPath curator and experimental tables) — the cases are the same; the NETWORKS are new, which is what structural rules are claimed over.
- **Success statement (fixed now)**: the curator held-out accuracy on 98 is within 2.0pp of the frozen 89.77% (identical-content regeneration noise is ~0.02pp; cross-release content change is real, so the band is wide and stated in advance). Experimental reported beside it, RAF set aside both figures. Outside the band, in either direction: trace before any claim (specs/108 discipline for gains too — a jump can be new curation answering cases for free, which is Reactome improving, not us).
- No rule, cutoff, or table is refitted before R1 is recorded.

## The rest of the coordinated move (after R1 is recorded, whatever it says)

- R2: parity harnesses re-run on the 98 build (pins 140, readouts 141, rails 142) — the engine-serves-the-bench invariant must survive the release.
- R3: the arrows pipeline on 98 (census already in hand: 267 arrows), source build, composites, the specs/148 scoring repeated; `src/api/leads.json` regenerated with 98 provenance.
- R4: the confidence table refit on the 98 record (its provenance test forces this before any 98 catalog is served); scope statement updated with whatever R1 says, including a changed number if the number changed.
- R5: serving: a packaged 98 catalog; `current` stays Adam's call.

## Mechanics

- Graph DB: stop `reactome-neo4j` (Release97; kept for rollback), run `graphdb:Release98` on the same ports with FRESH volumes (the 97 container's named volumes would shadow the baked 98 store). Verify `DBInfo.releaseNumber = 98` before anything reads it.
- Diagrams: `LNG_DIAGRAM_DIR` pointed at a Release-98 mirror (fetched from the download server; the generator needs covering diagrams for every built pathway, not only the 92).
- Catalog discipline unchanged: builds under `~/deltasignal-catalogs/builds/`, BUILD.json provenance, `current` untouched until Adam moves it.

## Non-goals

- No new rules, no tuning, no gate-shopping on the 98 numbers in this spec.
- No Lun 2019 (sealed; Adam's alone).
- No change to the frozen 97 record — it remains the published provenance for everything served today.
