# specs/153: rating the unrated arrows on decryptM

## Why

The served leads table (specs/149, re-measured specs/152 G5) carries 16
curator-drawn arrows into 10 centres; only the p-AKT arrow has a measured
rate class (down-call precision 0.24 on korkut). The other 15 say "unrated —
no dataset in the record perturbs upstream of it yet". A lead without a rate
is honest but less useful; decryptM (FIT since its first scoring, specs/021)
perturbs kinases with drugs and reads phospho-sites, i.e. exactly the
upstream-perturbation shape the korkut record lacked for these arrows.

Rating is honesty metadata fit on FIT data — the same footing as the
confidence table (specs/146, fit on tuning). No adoption gate is at stake;
the output is per-arrow rate classes, or the honest finding that decryptM
cannot reach them either.

## Protocol

- **S1 — structural coverage screen (free: uses no measurements).** For each
  arrow (source pathway, entity, centre): does any decryptM drug's measured
  kinobeads target (the specs/051 R4 map) resolve to a pinnable root in the
  arrow's SOURCE pathway, and does any decryptM phospho-site resolve in the
  CENTRE network? Count (drug, arrow, centre-site) combinations. Arrows with
  zero combinations are closed as "unratable on the current record" — a
  stated fact for the leads table, not a failure.
- **S2 — registered measurement (only for arrows S1 finds covered).** Run the
  decryptM axis on the specs/152 composite catalog (arm) and on the sound
  92-build (control), F discipline as specs/148: the centre's own answer
  where it routes; bridged answers only for control-untouched readouts.
  Per covered arrow: the analog of down-call precision — of the bridged
  change calls on control-untouched sites, the share matching the measured
  regulation — with an exact binomial interval. A rate on fewer than 20
  bridged calls is reported but the arrow stays "unrated" in the served
  table (interval too wide to be a rate).
- **S3 — serving.** Arrows gaining a rate class get it in `src/api/leads.json`
  via the usual PR with provenance (arm names, n, interval); everything else
  unchanged.

## Non-goals

No new arrows, no rule changes, no re-scoring of korkut or the benchmarks;
decryptM stays FIT (this is a fitting-side use).

## Amendment 1 (2026-11-13): S1 result, registered S2 arms

S1 (s1_coverage.py, s1_coverage.tsv; anchors expand a set-valued arrow entity
to its Neo4j leaf descendants): **six of 16 arrows are decryptM-coverable**,
with per-centre attribution unambiguous (in every covered centre only one
arrow's source carries the covering drug's targets):

- PTEN Regulation -> PTEN -> PIP3 (AZD8055, Dactolisib via MTOR; MK2206 via
  AKT1-3; 50 centre sites) — unrated today.
- Hedgehog 'off' -> GLI3R -> RUNX2 (Tideglusib via GSK3B, the kinase that
  makes GLI3R; 10 sites) — unrated today.
- DSB Response -> p-S456-ABL1 -> RUNX1 (Dasatinib, Imatinib, Nintedanib via
  ABL1; 6 sites) — unrated today.
- The three p-AKT arrows (ERBB2: 9 drugs/21 sites; Intrinsic Apoptosis: 1/7;
  Mitotic G1: 2/29) — already rated on korkut; decryptM gives an independent
  rating.

The other ten arrows have no decryptM drug reaching them: they stay
"unrated — no dataset in the record perturbs upstream" and S1 is their
closure on the current record.

S2 arms (registered here before running): control `r98fix_dm153` on
`20261113-1311_baf608b_r98fix` and arm `arrows98fix_dm153` on
`20261113-arrows98fix-10c`, both on the SAME case file
`cases_153_centres.tsv` (the kinobeads case table restricted to the ten
centres: 658 cases, 9 centres with cases, 15 drugs) so denominators are
equal by construction. F discipline as specs/148: the centre's own answer
where it routes; a bridged answer counts only on control-untouched cases.
Per covered (arrow, drug set): bridged change calls on control-untouched
sites vs the measured regulation, exact binomial interval; fewer than 20
bridged calls -> reported, arrow stays unrated in the served table.
