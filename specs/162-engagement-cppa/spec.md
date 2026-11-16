# Feature Specification: Engagement pins, replicated on CPPA

**Feature Branch**: `specs/162-engagement-cppa`

**Created**: 2026-11-16

**Status**: Draft, registered BEFORE any CPPA value is read for this spec

**Input**: Adam, 2026-11-16 ("ok" to the recommendation after specs/161): replicate engagement-scaled pins on an independent dataset under the new false-alarm gate rule (CLAUDE.md, the 2026-11-16 protocol rule).

## Why CPPA

specs/161's engagement came from korkut's own markers, on korkut. CPPA (Zhao 2020) is an independent dataset: different drugs, 317 cell lines, and per-(drug, cell line) samples. It also has a REGISTERED direct-marker definition fixed before any value was read, the specs/139 `markers.py` (CPPA phospho antibodies on sites Reactome curates as made by a reaction the drug's primary target catalyses). CPPA is FIT after specs/121's single test scoring, so this is a replication on new cases, not a held-out test. Under CPPA's terms, data and derived numbers stay outside the repository, which carries the design and the verdict only.

## The registered design

**Population**: single-compound CPPA samples whose drug has at least one specs/139 marker and an inhibitor role in the specs/121 drug-target table: 61 drugs, 530 (drug, cell line) pairs, 3,938 samples.

**Engagement** (per (drug, cell line), the protocol's "you measured your drug's direct marker"): e = min(1, 2^m), with m the mean of the drug's marker delta-p (log2) over that pair's samples (all times and doses). Computed by a committed script before any non-marker value is read.

**Arms** (the same runner and the same /api/readouts valuation for both):
- baseline: the recorded protocol, the drug's target pins at 0 (one solve per (drug, pathway), shared across cell lines);
- engagement: pins at e(drug, cell line) (one solve per (drug, cell line, pathway)).

**Scoring unit**: (sample, antibody) for antibodies in the CPPA readout map, EXCLUDING the drug's own markers in both arms. Measured call: delta-p ≥ log2 1.5 is up, ≤ −log2 1.5 is down, otherwise not (CPPA's existing 50% secondary convention, specs/137). Predicted call at the bench's 0.85 / 1.15. A prediction is the mean over the pathways that score the (drug, readout), as specs/121 does.

**Gate** (the 2026-11-16 rule for false-alarm fixes, registered here):
1. paired: among (sample, antibody) units whose call moves, fixed > broken by a one-sided sign test at p < 0.05; AND
2. guard: sensitivity on reachable units (baseline fold ≠ 1; real measured changes caught with the right direction) drops by at most 0.05 absolute.

Reported whatever they show (verdicts in the repo, numbers private): balanced accuracy, specificity, precision, AUROC on predicted |log2 fold|, the specs/121 primary endpoint (per-sample correlation averaged per drug), and the anatomy by engagement level.

## Requirements *(mandatory)*

- **FR-001**: the engagement script, the arm runner, and the scorer are committed and fixture-checked before the run; role guard on the CPPA datasets (FIT).
- **FR-002**: the env is verified by the specs/157 positive probe; the arms differ only in pin values.
- **FR-003**: all CPPA-derived numbers are written under `~/deltasignal-catalogs/analysis/162/`; the repository receives the verdict only.
- **FR-004**: a pass licenses a scope-statement line for the engagement protocol (Adam's call to serve); a fail records the negative and the protocol stays documented guidance.

## Success Criteria *(mandatory)*

- **SC-001**: the verdict is on the record against the registered gate, with no number in the repository.
- **SC-002**: any served claim cites this replication.

## Non-goals

- No engagement PREDICTION; no CPPA number in the repository; no new solver behaviour.
