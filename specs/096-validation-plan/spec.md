# specs/096: how DeltaSignal's accuracy will be validated

**Status:** plan, agreed 2026-10-23, committed before the rules and the test it governs.

Adam, 2026-10-23: "I want to continue with new rules until we get it handling everything we can
find with our tests, and then validate this with an external dataset. Curators could also make
one or two more pathways and see how well this works." And: "look through all the rest of the
pathways and make some sort of observation on which loops exist and see if we are handling that
type of loop." Agreed in that order ("ok").

## Why a plan

Every rule since September was chosen by looking at errors on the benchmark. The curator
held-out split has been a selection set since then (specs/051); the 849 experimental cases have
chosen rules too, and korkut2015 and decryptM 1 have each been scored once and are now fit
material. So the current figures (curator 88.28%, experimental 74.68%, docs/RESULTS.md) are
**development-set accuracy**. They are the right thing to optimise and the wrong thing to claim
as generalisation. This spec fixes, in advance, what will be claimed and from what data.

## 1. Development (continues now)

Rules continue as before: traced case, Neo4j check, pre-registration, catalog-wide rule,
both axes, gate, Adam decides. Three additions:

- **Targets.** Errors on which the curators and the experiment agree (78 experimental at
  specs/094), curator-only errors with a verified mechanism, and unhandled loop types from the
  loop audit (section 2). **Not targeted:** the 142 experimental errors where curators and
  experiment disagree. Fixing those chooses one reference over the other.
- **Reach.** Each adopted rule's research.md records whether it changes any network outside the
  pathway that prompted it (census), so its evidence base is visible.
- **Stop rule.** Development stops, and the model is frozen, when either:
  1. every loop type in section 2 that occurs in at least 5% of the audited pathways has an
     adopted rule or a recorded refutation, **and** the remaining agreed errors have no traced
     mechanism verified in Reactome; or
  2. four consecutive pre-registered arms fail their gate.
  Adam may stop earlier.

## 2. Loop audit (structural coverage)

Every cyclic component of every network is labelled by the handled loop types that explain its
edges: interconversion pools (specs/039), closed and open carrier loops (049, 077), release and
recycle pools (081, 086), Goodwin rings (066), self-contained inhibitors and activators (022,
083, 095). Cycles left once those are removed are unhandled, grouped by shape. For each type:
its frequency in the 92 benchmark pathways and in a sample of other Reactome pathways (no ground
truth, a variant catalog), and the accuracy of benchmark cases routed through it.

This measures **structural coverage**, not accuracy: a handled type can still be wrong in a new
pathway. It says whether the benchmark is representative and where rules are missing. Its own
spec (specs/097) carries the method and numbers.

## 3. Freeze

At the stop: tag the deltasignal and generator commits, build one canonical catalog from the
tagged generator, record every default (`DS_*`, `LNG_*`) and the benchmark protocol
(`DS_PIN_SCOPE`, `DS_PIN_ROOT_IGNORE`, thresholds DOWN < 0.85 ≤ NORMAL < 1.15 ≤ UP, the
perturbation levels). No rule, default or threshold changes after the freeze for the test.

## 4. External dataset (scored once)

- **Choice before looking.** Candidates are public perturbation datasets (knockdown, knockout,
  overexpression or inhibitor; proteomic, phosphoproteomic or transcriptomic readouts mappable to
  Reactome entities). A shortlist is drawn up from descriptions only. Adam approves the dataset
  before anything is downloaded (no outside contact without approval).
- **Pre-registration before scoring:** the mapping of perturbations to pinned nodes and of
  readouts to network nodes, the fold-change call rule, the metric (accuracy and macro-F1 on
  DOWN/NORMAL/UP), the baselines (reachability, always-NORMAL, naive signed propagation), and the
  comparison to MP-BioPath where possible. Committed before the data is scored.
- **One scoring** of the frozen model. Reported whatever it shows. Afterwards the dataset becomes
  development material.

## 5. Blind curator pathways (prospective)

- One or two Reactome pathways **not among the 92**, chosen with the curators.
- The curators write predictions in the existing curator format (gene, perturbation,
  readout, expected DOWN/NORMAL/UP) **before** any model output for those pathways exists.
- The frozen model is run once and scored. A few hundred cases at most, so this tests whether
  accuracy holds, not its precise value. Any contact with curators is Adam's.

## 6. What will be claimed

- Development figures are reported as such ("development set"), never as generalisation.
- Generalisation claims come only from sections 4 and 5, with their case counts and intervals,
  next to the development figures.
