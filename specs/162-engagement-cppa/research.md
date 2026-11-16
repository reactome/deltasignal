# specs/162 research: the verdict (2026-11-16)

Run on the served 98 build (four shards, every server probed, every registry
rule canonical). Scored by `bench/analysis/162/score_162.py` (fixture-checked
and committed before any arm finished). Numbers are kept outside the
repository under CPPA's terms (`~/deltasignal-catalogs/analysis/162/score.json`).

**Registered gate: FAIL.** The paired part passes (many more moved calls are
fixed than broken), but the sensitivity guard fails: reachable sensitivity
drops by more than the registered 0.05.

**The reading.** On CPPA, engagement-scaled pins remove correct and false
change calls at the same rate. Precision of change calls does not improve,
and ranking skill (AUROC on predicted |log2 fold|) falls. The fixed-vs-broken
excess reflects that most units are true "no change", so calling fewer changes
fixes more units than it breaks. That is a blanket shift toward "no change",
not added information. The korkut result (precision and AUROC both up) does
NOT replicate on an independent dataset.

**Consequence.** No validated claim. The docs note on pinning at measured
engagement is amended to say the replication failed, and the protocol stays
an option, not a recommendation.
