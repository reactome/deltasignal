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
