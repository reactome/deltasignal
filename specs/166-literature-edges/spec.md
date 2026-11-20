# specs/166: literature edges from FullTextMining, registered comparison

Status: REGISTERED 2026-11-18 before any arm was scored. Only the smoke test had run: every pathway
parses and solves in the `with` and `without` arms. **RESULT: FAIL** (see research.md).

## Question

Do edges text-mined from the literature, added to Reactome's curated networks, improve the
predictions? Or do they only add connectivity, the way any extra edges would?

## Material

- **Bundles.** FullTextMining's `composite-h12b`: 10 pathways (R-HSA-1257604, 194138, 3700989,
  388841, 453279, 5693532, 6785807, 75153, 8878166, 8878171).
  - Built with logic-network-generator `d6d026b` (code-identical to the served build's `baf608b`)
    on Reactome Release 98.
  - Literature from FullTextMining `d20df9e`, run label `h12b`, 256 papers, mostly innate immunity
    and metabolism.
  - 15,623 literature edges: 6,683 verified, 8,940 extracted.
  - `MANIFEST.sha256` verifies for all 300 files.
  - Stored at `~/deltasignal-catalogs/external/literature-h12b/` (outside the repository).
- **Excluded papers.** The papers were not chosen from our benchmarks. Korkut 2015 (PMID 26284497),
  Zecha 2023 (decryptM, 36926954) and Lun 2019 (31101498) are excluded on their side, and none is
  among the 256.
- **Dataset roles.** korkut and MP-BioPath experimental are FIT, and curator held-out is SELECT, so
  this is a measurement, not a test.
- **The `without` arm matches our served build.** All 10 pathways have the same edge counts as the
  served Release 98 build (`20261113-1311_baf608b_r98fix`). Three pathways list a few extra ids in
  `stid_to_uuid_mapping.csv` and `nodes.csv` (5, 3 and 6; simple entities and dissociation sinks),
  but no edge uses them, so they are isolated and inert. Every arm is built from their bundle in
  any case.

## Arms

All arms are built by `bench/analysis/166/build_arms_166.py` and served by `scripts/run_arm.sh
--catalog` with the default solver at this branch's commit.

| arm | content |
|---|---|
| `without` | Their bundle minus the literature rows (matched on source, target, edge_type and pos_neg from `literature_origin.csv`) and minus the literature nodes (`context_origin.csv`). This is the generator's bundle as built. |
| `with` | The bundle as delivered. |
| `shuf01`..`shuf19` | The literature edges rewired at random, 19 seeds. Within each class (edge_type, pos_neg, and_or), targets are permuted among that class's edges, so each node keeps its in- and out-degree per class and each edge keeps its sign and type. Self-loops, duplicates and collisions with curated pairs are redrawn. Curated rows are untouched. |
| `only-<slice>` | Descriptive ablations: `without` plus only the literature edges in one slice. Slices are layer (verified, extracted); sign_source (stated, inferred, na); group (binding/other, regulation, expression); and prior_confidence (< 0.4, 0.4–0.6, ≥ 0.6). |

Literature nodes take part under the bench's existing protocols. **Correction (2026-11-18, after
registration, before any arm was scored).** The registered text said literature nodes would be
pinned and read as readout copies. They are neither: both benches resolve genes and readouts to
Reactome stable ids from Neo4j (`neo4j_gene_to_stids`, readout stids), and literature nodes carry
`UNIPROT:`/`LIT:` ids, so they only add routes between curated nodes.

Literature edges still change the pins in one way. An edge into a curated root raises its in-degree
(`root_indegree`), so that root is no longer pinned, and the gene is pinned by the specs/024
fallback or its case becomes invalid. The shuffles do the same. The intersection keeps only cases
valid in every primary arm. The number of cases whose pins differ from `without` is reported per
arm.

## Cases

Cases are scored on three axes:
- MP-BioPath experimental, which covers only 3 of the 10 pathways (PIP3, TP53, Mitotic G1);
- korkut, with each drug set's own markers excluded (`marker_names()`, specs/161);
- curator, on all 10 pathways, with held-out and tuning reported separately.

The case set is the cases scored in every primary arm (`with`, `without`, and all shuffles). Each
arm's denominator before the intersection is reported.

## Registered gate (primary)

The primary pool is korkut plus experimental cases. For each arm X, net(X) is:
- the number of cases correct in X and wrong in `without`,
- minus the number wrong in X and correct in `without`.

The literature edges PASS only if both of these hold:
1. **net(with) > 0**: the real edges help.
2. **net(with) > net(shufN) for every one of the 19 seeds**: they beat random edges with the same
   degrees and signs. This is a rank test: p = (1 + #{seeds with net ≥ net(with)}) / 20, which
   must be ≤ 0.05.

Curator numbers do not enter the gate. They are reported with the SELECT caveat.

## The three numbers (reported, `with` against `without`, per axis)

1. **Coverage gained.** Cases whose `without` prediction is exactly the baseline fold (1.0) and
   whose `with` prediction is not. Gains on the korkut and curator/experimental unreached lists
   handed over on 2026-11-18 are counted separately.
2. **Newly covered.** Among the coverage gained: accuracy, change-call precision, and the calls'
   distribution. These are reported as their own labelled class.
3. **Previously covered.** For cases already off 1.0 in `without`: fixed, broken, net, and accuracy
   in each arm. This is the damage check.

Each shuffle's three numbers are reported beside the real arm's, as context.

## Descriptive (no gate)

- The `only-<slice>` arms' net and three numbers. With 11 slices, some will look good by chance, so
  none of them is a finding.
- **The 105-pair question** from FullTextMining's 002 message. Some korkut pairs (drug target,
  readout) have a directed curated path from a node containing the target to a node containing the
  readout, but are predicted at exactly 1. On the `without` arm, each such pair is classified by
  where the perturbation stops. The classes:
  - no pinned root form on any path;
  - the path runs only through edge types that do not carry signal;
  - a held node (cofactor, drug, carrier rule);
  - the readout form or copy is not on the path;
  - other, traced.

## Outputs

- Per-arm case dumps are written under each build's `results/<sha>/<arm>/`.
- The scorer is `bench/analysis/166/score_166.py`. Its output goes to `research.md` in this
  directory, whatever it says.
- Results are posted to FullTextMining as numbered files in the shared nibi mailbox.

## Amendment 1 (2026-11-19, during the run): a time cap on shuffle arms

Some shuffled networks solve pathologically slowly.
- shuf02 failed twice on a single Mitotic G1 solve: first past the 600 s solve timeout, then past a
  3,600 s one.
- shuf01 spent 2 h 24 min on the curator axis's Mitotic G1 cases without finishing.
- shuf03 cleared the same pathway in 142 s, so the slowness depends on the seed.

Rule from here on: a shuffle arm still running after 9,000 s is stopped and recorded as not
finished.

**Revised at 03:15.** shuf03 was stopped at 9,000 s with only TP53 left on the curator axis (DSB
Repair alone took 5,393 s), so a shuffle arm needs about 4–6 h and a 9,000 s cap would let none
finish. The cap is now 21,600 s. The run order is now:
1. the 11 descriptive slices, which are cheap;
2. shuf08–19;
3. a retry of shuf01–03, which were stopped under the first cap.

At most four arms run at once (scheduler log: `analysis/166/run2/sched.log`). The gate's rank test uses the seeds that finished, and every unfinished seed is listed.

Why this cannot change the verdict:
- Condition 1 (net(with) > 0) was decided before this amendment. It needs only `with` and
  `without`, which both finished.
- The shuffles answer only the second question: are the real edges better or worse than random
  ones?

What this does change: the seeds that finish are a selected subset, the ones the solver can handle.
This is recorded beside every shuffle number.

**Revised at 06:05.** No shuffle finished the curator axis within 3–5 h: shuf04 stalled on Mitotic G1
for over 4 h, and shuf05 and shuf06 were still on TP53. Curator is reported, not gated. From here on,
shuffle arms run only the gate's axes, experimental and korkut, under run name `lit166p`. Curator
results from the shuffles are kept as partial, descriptive per-pathway numbers in `lit166/curator.log`.
shuf04–07 were stopped mid-curator and requeued under the new rule. The 6 h cap still applies.

**Ended at 21:14 on 2026-11-19 (Adam, option 2).** Four seeds finished: shuf08, shuf10, shuf12 and
shuf13. The rest are unfinished or were not run, and `research.md` lists each one. Condition 1
(net(with) > 0) fails, which settles the verdict. Condition 2 is reported with the four seeds:
`with` beats every one, at rank p 0.20. The orchestration that ran the arms (`sched_166.sh`,
`cap_166p.sh`) is committed beside the arms.
