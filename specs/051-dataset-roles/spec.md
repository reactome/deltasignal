# specs/051: dataset roles — what fits the model, what judges it

**Status:** plan, 2026-10-12. Nothing is fitted here.

**Origin:** Adam, 2026-10-12: "At some point we want to use datasets to
learn parameters in the pathways to improve the inference. Other datasets
we want to use to show how accurate the pathways are."

## Why this has to be written down first

The project already learns parameters. It does it by hand, one `DS_*`
default at a time, and specs/009–049 are the record. What it has not done
is keep the data that *chooses* a parameter apart from the data that
*reports* accuracy. Three things from the last weeks show the cost.

1. **The curator "held-out" split is no longer a test set.**
   - It is held out of the *tuning pathways*. But since 2026-09 every
     adopted rule has been chosen partly on its held-out number:
     specs/023, 039, 040, 046, 049 and others.
   - A number used to pick among arms is a *selection* number. Reported as
     a test result, it would be optimistic, by an amount we cannot state.
2. **decryptM was clean for exactly one scoring.**
   - Specs/021 scored the frozen Release97 rules on it once.
   - The trace that followed (the EGFR CBL-cycle rail) makes any rule
     written from that trace a rule informed by decryptM.
3. **Cases are not independent.**
   - The ERK1/2 activation loop appears in dozens of pathways. The 1,862
     per-pathway decryptM cases are 1,165 distinct (drug, site)
     measurements.
   - A random split over cases would put the same measurement on both
     sides.

## Roles

Each dataset (or a declared slice of one) has exactly one role at a time.
The register below records the role, and every change of role with its
date and reason.

| role | may be used to | may not |
|---|---|---|
| **FIT** | fit any parameter, trace any case, motivate any rule | be reported as accuracy |
| **SELECT** | choose between registered arms or fitted candidates; trace cases | have parameters fitted to it directly; be reported as an unbiased test |
| **TEST** | be scored once per frozen release, under a protocol committed before the score | be read case by case before the score is committed; be used to choose anything |

**Demotion is one-way.** After a TEST score is committed, cases may be
traced (a negative result has to be explained). But the first rule or
parameter change motivated by that dataset demotes it to SELECT for the
next release. A TEST set is a consumable, so the register keeps a reserve
of sets that have never been read.

**Split units, against leakage.**
- A TEST case must not share its *(perturbation, readout)* pair with any
  FIT or SELECT case. A perturbation is a gene, or a drug's target set; a
  readout is an entity or phosphosite.
- Overlap in cell line and lab is recorded and reported, not forbidden.
- Counts are given in distinct (perturbation, readout) units as well as
  per-pathway cases.

**Data that can be used at all.**
- Cell lines and other non-human-subject sources only.
- No patient-derived, personal or health data. That excludes tumour
  cohorts such as CPTAC and TCGA.
- Openly downloadable, with the licence recorded.
- SHA256 of every downloaded file in the register's data directory.

## Register (current)

| dataset | slice | role | since | note |
|---|---|---|---|---|
| MP-BioPath curator | 10 tuning pathways | FIT | 2026-05 | the paper's ten |
| MP-BioPath experimental | 849 cases, tuning pathways | FIT | 2026-05 | read on every arm |
| MP-BioPath curator | 70 held-out pathways | **SELECT** | 2026-10-12 (relabel) | used for adoption decisions since 2026-09; still reported, but as selection accuracy |
| decryptM 1 (PXD037285) | all | **FIT** | 2026-10-12 | scored once as TEST (specs/021, 9cb8983), then traced; dose response makes it the richest fitting signal we have |
| *TEST reserve* | — | TEST | — | to be filled from the survey in research.md. Nothing is read before it is registered |

Releases change content, so a dataset scored against Release 98 is scored
on the frozen Release 98 model under the specs/050 procedure.

## Protocol (Adam agreed, 2026-10-12)

These rules bind every fit and every score from here on.

1. **Split by group, never by case.** The groups are:
   - perturbation (gene, or drug target set);
   - readout (entity or site);
   - pathway;
   - cell line.

   Each claim names the group it generalises over. "Works on unseen
   pathways" needs a pathway-grouped split. Counts are given in distinct
   units as well as cases.
2. **Selection by grouped cross-validation inside FIT + SELECT.**
   - A fit is chosen by leave-group-out cross-validation: by pathway group
     and by perturbation, with the folds fixed by a seeded, committed fold
     file (`bench/splits.py`).
   - A single fixed selection set wears out; the curator held-out split
     already has.
   - SELECT data enters the cross-validation as folds. It is never fitted
     whole.
3. **TEST is opened once per frozen release,** under a protocol committed
   first.
   - The strongest TEST is prospective: newly curated pathways in the next
     Reactome release (the specs/050 procedure), and data published after
     the model froze, such as decryptM 2.0.
4. **Every score carries these baselines, on the same cases:**
   - always-majority;
   - reachability (graph only);
   - the unfitted defaults, for any fitted model;
   - a **rewired-network null**: degree-preserving rewiring, edge type and
     sign kept, several seeds. It shows the curated wiring carries the
     signal.
   - Every fit also runs a **label-shuffle control**: fitted to shuffled
     FIT labels, it must gain nothing on held-out folds. If it does, the
     procedure leaks.
5. **Few, shared, anchored parameters.**
   - Parameters are shared by edge type or motif before any per-edge
     parameter is allowed.
   - Each fit is regularised toward the current default.
   - Parameter pairs that trade off are reported jointly.
   - Fitted values carry bootstrap intervals over perturbations, and are
     checked for biological plausibility before adoption.
6. **Statistics.**
   - Macro-F1 on direction is the objective and the headline.
   - McNemar for paired comparisons.
   - **95% intervals by bootstrap over perturbation groups**, never over
     cases.
   - Concentration: pathways and distinct readouts moved.
   - A count of every arm and fit evaluated against SELECT data, kept in
     the spec.
7. **Data definitions are fixed before results.**
   - "Regulated" comes from each dataset's own statistics.
   - Drug → target comes from measured affinity (Klaeger 2017 kinobeads),
     with the cutoff registered before use. The hand-written
     `DRUG_TARGETS` table is retired once that map exists.
   - Entity and site mappings are pinned to the Reactome release.
   - Licence, checksums and the processing script are kept for every file.
   - Cell lines only.
8. **Cell context.**
   - A flat readout in a cell line whose pathway is off is not a model
     error.
   - Until context is learned (L3, from FIT data only), unchanged cases are
     secondary.
   - Scoring may be restricted to cell lines where the pathway is shown
     active, by a criterion registered beforehand.
9. **Reproducible.**
   - One command regenerates each reported number from pinned data, code
     and catalog build.
   - Negative results stay in the record.

**Gate before any fit (L1 onward):** confirmation that the organisation's
data-use policy permits fitting model parameters to public cell-line data.

## Stages of learning

Each stage is registered before it runs. It is fitted on FIT, chosen on
SELECT, and scored once on TEST.

- **L0. Prerequisites.**
  - Done: label invariance (specs/013), deterministic sweeps, the bench
    refusals.
  - Missing: this register, and a guard that refuses to fit on a non-FIT
    dataset.
- **L1. Global parameters that already exist as `DS_*`.**
  - Examples: `DS_DEPLETION_H_MIN`, `DS_SELF_INHIBITOR_WEIGHT`,
    `DS_CYCLE_PHI`, the loop elasticity ε.
  - This formalises what specs/009–049 did by hand. Expected gain: small.
  - Its value is turning the hand search into a reproducible fit.
- **L2. Per-edge-type parameters.**
  - Inhibition strength, the depletion ceiling, AND exponent, and loop
    elasticity, each by edge type and motif (specs/014 §6, specs/050
    register).
  - Tens of parameters, fitted with regularisation toward the current
    defaults.
- **L3. Cell context.** Which pathways are active in a cell line. Without
  it a flat phosphosite cannot tell "no effect" from "pathway off"
  (specs/021, unchanged cases). It needs a context source from FIT data,
  such as baseline phospho or expression.
- **L4. Per-edge weights.**
  - Hundreds to thousands of parameters. They need gradients or a
    derivative-free search, and much more FIT data.
  - `DS_SCC_METHOD=minimize` is kept reachable for this (specs/003).

**The fitting objective is macro-F1 on direction**, as now. Where a
dataset carries magnitudes (decryptM dose response), a magnitude loss is
added at a weight registered before the fit.

## First steps

1. Survey public perturbation datasets for the TEST reserve and for FIT
   volume (research.md). *Done, R1.*
2. A machine-readable register, `bench/datasets.tsv`, and
   `bench/dataset_roles.py` with `require_role(dataset, purpose)`.
   - Any fitting script calls it, and it refuses a dataset whose role does
     not allow the purpose.
   - The scoring bench writes the role it scored into its summary.
3. Pick and register two TEST sets from different labs or technologies.
   Download them with checksums, without reading their values.
4. The scoring tools the protocol needs:
   - perturbation-clustered bootstrap intervals;
   - the grouped fold file;
   - the rewired-network null.
5. L1 as the first registered fit, behind the policy gate.

## Rewired-network null: design (registered 2026-10-12, before any build)

Protocol point 4 asks whether the curated wiring carries the signal. A
rewired network keeps every node, every edge class and every node's degree,
and moves only which reaction each edge reaches. If the real network does
not beat it, the score comes from degrees and pinning, not from Reactome.

**What is rewired.**
- Only the reaction-level classes in `logic_network.csv`: `input`,
  `output`, `catalyst`, `regulator`, `depletion`.
- A class is (`edge_type`, `pos_neg`, `and_or`), so an inhibition stays an
  inhibition and an AND input stays an AND input. There are 7 such classes in
  the build.
- **Kept fixed:** the entity-structure edges, meaning `boundary_edges.csv`
  (assembly, dissociation, set_member) and `variant_pool` / `variant_split`.
  They say what an entity is made of, not what reacts with what. Rewiring
  them would test a different claim.

**How.**
- Repeated double-edge swaps within a class (Maslov–Sneppen):
  (a→r1), (b→r2) becomes (a→r2), (b→r1).
- Every node keeps its in- and out-degree in every class. A swap is
  rejected if it makes a self-loop or repeats an edge already in its class.
- The target is 10 × |class| accepted swaps, tried at most 100 × |class|
  times. The achieved count and the share of edges whose target moved are
  recorded per pathway.
- The RNG is seeded from the string `"{seed}:{pathway}:{class}"` over edges
  sorted by (source, target), so a build is reproducible from its seed.
- A moved edge keeps its source, stoichiometry and attributes. Its
  `edge_reaction_id` is blanked, because that column names the reaction and
  would now be stale. The solver and the bench do not read it.

**Side tables.**
- `pools.csv`, `pool_transitions.csv`, `pool_carriers.csv` and
  `carrier_loops.csv` name specific steps, built from Reactome chemistry
  (`find_pools`). They are **deleted** from the rewired build, not emptied,
  so a default-config run refuses instead of silently solving with zero
  pools.
- Entity-level tables (`containment`, `cofactors`, `drugs`, `pathogens`) and
  the id maps stay. Node identity is unchanged.

**Matched control.**
- Every rewired build is scored with `DS_CYCLE_MODE=off` and
  `DS_CARRIER_LOOP_MODE=off`. The real build is scored under the same two
  overrides (arm `wire0`), so the only difference is the wiring.
- The canonical score is reported beside it, for scale.

**Pins and readouts come from the real network.**
- In-degree is preserved, so `root` pins are identical by construction. The
  `root_cycle` supplement depends on cycles and would change.
- The bench reads the pin network from the source build
  (`DS_PIN_NETWORK_ROOT`), so every arm perturbs the same nodes.

**Builds.**
- Five seeds (1–5), each a variant build `<source>_rewire<seed>` with a
  `BUILD.json` naming the source build, seed, classes, swap counts and the
  rewiring script's commit. `current` never moves to them.

**Axes and statistics.**
- Curator tuning and held-out (SELECT), and experimental, by `run_arm.sh`.
  decryptM is secondary, if the drug bench accepts the overrides.
- Per axis: accuracy and macro-F1 for `wire0`, each seed, and their mean.
- The paired difference `wire0` − seed is computed case by case, with a
  95% bootstrap interval over perturbation groups (gene; drug for decryptM).

**What the result can say (fixed now).**
- On an axis, the curated wiring **carries signal** if `wire0` macro-F1
  exceeds every seed's, and the interval for `wire0` − seed mean correctness
  excludes 0.
- If `wire0` is inside the seeds' range, the axis does not show that the
  wiring matters. That is reported as a headline negative, not explained
  away.
- The curator axes are expected to favour the real wiring strongly, because
  the curators reasoned from the same pathways. The experimental and
  decryptM axes are the informative ones.
- Nothing is adopted from this. It is a baseline, and it is counted in the
  SELECT arm tally (point 6).

## SELECT arm tally (point 6)

Every arm or fit scored on SELECT data (curator held-out) from 2026-10-12,
when the role was assigned. Earlier arms are in specs/009–049 and are not
recounted.

| date | spec | arm | purpose |
|---|---|---|---|
| 2026-10-12 | 051 | `wire0`, `rw1`–`rw5` | rewired-network baseline; nothing selected |
| 2026-10-13 | 052 | `zs0`, `zs` (f65618a; a void pair at 6273941) | zero-supply rule; not adopted |

## Non-goals

- No fitting in this spec.
- No change to what any existing number means. The curator held-out
  numbers stay as they are, renamed for what they now are.
