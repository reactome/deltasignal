# specs/166: literature edges from FullTextMining, registered comparison

Status: REGISTERED 2026-11-18, before any arm was scored. Only the smoke test had run: every pathway
parses and solves in the `with` and `without` arms; no case was benchmarked.

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
- **A difference from our build.** In 7 of the 10 pathways, the `without` arm has a few curated
  nodes that our served build lacks; nothing is missing. The cause is probably the build
  environment (ours sets `LNG_EMIT_ONE_SIDED=1` and `LNG_DIAGRAM_DIR`). It does not bias the
  comparison, because every arm is built from their bundle.

## Arms

All arms are built by `bench/analysis/166/build_arms_166.py` and served by `scripts/run_arm.sh
--catalog` with the default solver at this branch's commit.

| arm | content |
|---|---|
| `without` | Their bundle minus the literature rows (matched on source, target, edge_type and pos_neg from `literature_origin.csv`) and minus the literature nodes (`context_origin.csv`). This is the generator's bundle as built. |
| `with` | The bundle as delivered. |
| `shuf01`..`shuf19` | The literature edges rewired at random, 19 seeds. Within each class (edge_type, pos_neg, and_or), targets are permuted among that class's edges, so each node keeps its in- and out-degree per class and each edge keeps its sign and type. Self-loops, duplicates and collisions with curated pairs are redrawn. Curated rows are untouched. |
| `only-<slice>` | Descriptive ablations: `without` plus only the literature edges in one slice. Slices are layer (verified, extracted); sign_source (stated, inferred, na); group (binding/other, regulation, expression); and prior_confidence (< 0.4, 0.4–0.6, ≥ 0.6). |

Literature nodes take part under the bench's existing protocols. A literature node that contains a
perturbed gene and has no in-edge is a root, so it is pinned; and the readout rules read it like
any other copy. That is how the generator's tables present those nodes, and it is the same in every
arm that carries them.

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
