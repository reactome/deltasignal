# specs/040 — Fable derivation (blind; derivation-opus.md not read)

Prototype: scratch copy of `src` (not committed) in the session scratchpad
(`f040/src`), arms and census scripts in `raf/` (`f040_run*.sh`,
`census040_fable_v2.py`, `census_fable_v2.jl`). Build 20260928-1110_06ccb63.

## 1. One principle, two rules

specs/039: a pool cannot count its own output as supply. Generalised: **a
reaction must not read its own signal twice.** The RAF trace does it two
ways, so there are two rules; neither changes the network.

**Rule A — a self-fed input reads baseline (`DS_SELF_FED_MODE=entry`).**
- Scope: a strongly connected component C when it is iterated. Nodes owned by
  a balanced pool (states, intermediates, step copies) are 039's: never
  self-fed, and an entry iff the balance left them off baseline.
- *Signal-carrying entries*, fixed when C's iteration starts: nodes of C
  pinned off baseline, and unpinned targets in C with an activator outside C
  that is off baseline (upstream is solved). A pin **at** baseline (an inert
  cofactor inside C, a held drug) is a constant: not an entry, and the search
  does not pass through it. Two static definitions fail on RAF: "any outside
  activator" makes a binding-step copy an entry through ATP and flags
  nothing; "OR-alternative producers only" finds zero entries, because the
  loop is entered only through AND inputs of reaction nodes.
- Out(u): the targets in C that u feeds; a set-pool node is transparent, so
  the 80 variant copies of "Phosphorylation of RAF" are one consumer.
- u is **self-fed** iff every activator path from an entry to u inside C
  passes through Out(u): fed only by what it feeds. Covers an AND factor and
  a set-pool member.
- Every edge from a self-fed u into Out(u) reads the component-entry value
  (specs/018's `activator_break`/`supply` plumbing), which is **exactly
  baseline**: amendment 5's "re-evaluate with the pool at baseline" with the
  answer in closed form. A pin between consumer and u (MEK knocked out)
  makes u entry-fed, so it reads live.
- Per solve, deterministic, label-free (dominator sets). Under two
  simultaneous perturbations with one on the loop, the loop reads live
  again; the protocol is one root pin per case.

**Rule B — self-contained inhibitors at the leaf level
(`DS_SELF_INHIBITOR_LEAVES=1`).** Traced: the PEBP1 inhibitor
(R-HSA-5675413) is made from R-HSA-5675414, the RAF1-only sibling of the
step's input R-HSA-5672718; the input does not reach it. specs/022's pair test
(I contains a's stId, and a reaches I) gains a second clause: **leaves(I) ∩
leaves(a) holds a non-cofactor leaf ℓ, and a carrier of ℓ reaches both a and
I.** "One generation of common cause" was tried and missed PEBP1 (the shared
precursor is two up). Leaves are the containment table's childless stIds.
Formula, w = 0.1, weaken-only: unchanged; 022's pairs are a subset. It also
covers amplifier 4 without naming drugs: R-HSA-9657583 tracks the MEK step's
input exactly and R-HSA-9653109 tracks RAS:GTP at 100x; both are damped to
their untracked part, not held (the refuted specs/032 arm).

Both are exact at baseline (verified: max deviation at rest 2e-12 across 92
bundles) and byte-identical off (flags and pairs only).

**Relations.** 039: per node instead of per pool. 018: same read, structural
selection instead of an edge role. 022: same formula, containment where the
bundle decomposes. 013/014: elasticity damped every closure edge and churned
under relabelling; A opens only edges carrying no exogenous signal. **035**
(−122 experimental) held inputs *released unchanged* in every case and held
PIP3 itself; A holds nothing that carries signal — PIP3 is flagged under 1 of
322 root pins, PI(4,5)P2 under 118, each a pin that cannot reach it except
through its consumers. **036/037** rebundled nodes; **012** deleted a
catalyst-substrate term; A and B add flags to existing edges and can only
remove a duplicate reading.

## 2. Genuine feedback

MEK → RAF is real. At steady state a product loop enters as fold^k with k =
3 (pool) + 22 (leaves) − 1 (tracking inhibitor); k ≥ 1 has no finite fixed
point off the knife-edge, so the iteration finds the all-zero root, not
feedback. After A the self-fed exponent is 0: the endogenous part is bounded
at 1, the exogenous part (a MEK pin) reads live. A learnable ε ∈ (0, 1) on
exactly these edges is the future version; any ε now would be fitted.

## 3. Detectability

Solver only: SCCs (existing), entries from x vs baseline, one BFS per node
of C with Out(u) removed, O(|C|·|E_C|) per component. No generator input.

## 4. RAF worked (harness, canonical build, balance on)

Folds at the activated dimer (R-HSA-5672718) / binding step / dissociation /
p-T,Y MAPK monomer:

| case | ctrl | A | B | **A+B** |
|---|---|---|---|---|
| KRAS 80x | 9e-8 / 1e-9 / 5e-7 / 7e-26 | 100 / 2.47 / 0.006 / 1e-9 | 8e-8 / 1e-9 / 1e-7 / 3e-28 | **100 / 69.1 / 44.7 / 100** |
| BRAF 80x | 8e-8 / 5e-10 / 5e-7 / 7e-26 | 100 / 1.0 / 0.001 / 1e-12 | 9e-8 / … / 3e-28 | **100 / 63.1 / 45.5 / 100** |
| NF1 KO | 9e-8 / 4e-8 / 5e-7 / 7e-26 | 1.47 / 1.14 / 0.88 / 0.60 | 1e-7 / 4e-8 / 1e-7 / 3e-28 | **1.47 / 1.43 / 3.48 / 100** |
| KRAS KO | 5e-8 / 1e-7 / 5e-7 / 7e-26 | 0.354 / 0.71 / 1.41 / 4.0 | 35.4 / 100 / 63 / 100 | **0.354 / 0.379 / 0.033 / 1e-6** |
| HRAS 80x / KO | collapse | up / mixed | up / up | **92.8 → 100 / 0.65 → 0.004** |
| NF1 80x | collapse | 0.012 / 0.12 / 0.23 / 0.003 | up | **0.012 / 0.017 / 6e-7 / 1e-25** |

- Amplifier 1: the pool reads 1 under every RAS/NF1 case (80 under BRAF 80x,
  a pinned member). Amplifier 2: 21 leaves read baseline at the binding step.
  Amplifier 3: PEBP1 40.5 → damped, binding step 2.47 → 69. Amplifier 4:
  the drug copies, damped by B.
- **A alone** fixes the dimer only (the drug copies divide the MEK step by
  ~400). **B alone rails every KO to 100** (the loop runs to the other root).
  No single rule does it, as the problem predicted.
- The RAS pool is byte-identical to ctrl in every arm (pool-owned exclusion).
- **Open defect:** 3 of 8 cases end with residual ≈ 1 (HRAS 80x, KRAS 80x,
  NF1 KO), 5 with 0.003–0.01; ctrl converged. Something still oscillates,
  untraced. An arm must log convergence per case.

## 5. Outside RAF: census

Rule A (Python, activator SCCs, every single root pin; upper bound, since
"reached by the pin" stands in for "off baseline"): **67 of 92 pathways, 175
components, 1,803 nodes** flagged under some pin, 432 of 11,328 pins touch a
component; 165 flagged nodes are set-pool members. Largest: Class I MHC 286
nodes, Cell Cycle Checkpoints 173, Immune System 163, RAF 98 (14 of 325
pins, max 112 edges). Overlap with 035's held set by identity: 340 of 535
nodes, 256 of them Class I MHC — but by *case* the overlap is small (PIP3
above). Rule B (Julia, general clause): **538 → 1,406 pairs, 41 pathways
gain**; RAF 0 → 41, Chromatin modifying enzymes 0 → 308 (histone leaves;
promiscuous), ERBB2 0 → 75.

## 6. What can go wrong, and the gates

- A removes a railed switch's closing edge: gate, no pathway loses > 10.
- Multi-pin fallback to live: log `self_fed_edges` per case; report cases
  with 0 flags in a flagged component.
- B promiscuity (2.6× 022's pairs; Chromatin 308): B is its own arm; report
  the discordant cases per pathway; refuse adoption if one pathway carries
  > half the gain.
- Label dependence: the relabel check must not raise churn above ctrl's.
- Every pathway with zero flags (A) and zero new pairs (B) must be
  bit-identical to ctrl, or the change leaked.
- Class I MHC: 286 nodes and 151 pool members flagged; run it as a named
  concentration check.

Arms on one solver commit: ctrl, A, B, A+B. Decides: curator held-out and
experimental with McNemar p, noise floor 15, ≥ 2 pathways and ≥ 5
perturbations, net ≥ 0 outside RAF.
