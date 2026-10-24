# specs/098 research

## Arm `ubi`: void

`results/2d0c227/ubi` ran on a branch cut from a local `main` that had not been pulled after
PR #163, so it **lacked specs/095's default** (`DS_SELF_INHIBITOR_COPIES=on`), while its baseline
`results/b61b365/dflt095` has it. Its differences mix the two rules: the 9 PTEN overexpression
"breaks" it shows are specs/095's absence, not this rule (`pten_ub.jl` reproduces 83.5× with free
ubiquitin propagating). Not scored. The branch was merged with main (`specs/095` default on) and
the arm re-run as `ubi2`.

Lesson recorded: before an arm, check that the arm's commit contains the baseline's defaults
(`git merge-base --is-ancestor <baseline sha> HEAD`).

## Arm `ubi2`

Solver `dddb1b5` (contains specs/095's default), `--server DS_UBIQUITIN_MODE=inert`, against
`results/b61b365/dflt095`. Results `results/dddb1b5/ubi2`; 5,922 curator and 1,402 experimental
node-solves held.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | +6 (6 / 0) | 0.031 |
| curator tuning | +3 (9 / 6) | 0.61 |
| experimental (RAF set aside) | −4 (1 / 5) | 0.22 |

**The gate passes** (held-out +6, experimental −4 ≥ −5, no pathway below −10), narrowly.

### What moved

- **Gains.**
  - USP1 in Fanconi Anemia (held-out, 6): USP1 deubiquitinates FANCD2. Its knockout read DOWN
    because it starved the free ubiquitin pool and so lowered FANCD2 monoubiquitination; with the
    pool held, the knockout leaves ubiquitinated FANCD2 to accumulate and reads UP, as the
    curators say.
  - MYC and MAX overexpression (Mitotic G1, 4) and CCNE1 (S phase, 4 curator, 1 experimental) lose
    a false change (14.2× and 45× → 1.0×).
- **Losses.** CCNE1 overexpression at CCNA1 (Mitotic G1, 3: 14.2× → 1.0×), CCNE1 knockout in S
  phase (0.56× → 1.0×), and threshold shifts for CDK4, CCND1 and CDKN2A (0.82× → 0.87–0.93×).
  `ccne1_ub.jl`: CCNE1 reaches CCNA1 **only** through free ubiquitin (`R-HSA-68524`), made by
  "Degradation of ubiquitinated p27/p21 by the 26S proteasome" (R-HSA-187574) and used by
  "Ubiquitination of phospho-p27/p21" (R-HSA-187575). Degradation also releases the active
  cyclin E/A:CDK2, but in this network the only edge from the degradation step back into the p27
  ubiquitination cycle is the free ubiquitin node. Holding it cuts a curated route, so these
  losses are a real cost of the rule, not right-for-the-wrong-reason cases.

### Predictions

1. Only the 10 pathways with a ubiquitin-closed loop move a class: **holds** (Fanconi, Mitotic G1,
   S phase).
2. Of the curator cases crossing a ubiquitin-closed component, more fixed than broken:
   **holds** (15 / 6).
3. No prediction for held-out or experimental: held-out +6, experimental −4.

## Reading

Mixed. The rule is biologically sound for the pool (USP1 is a clean, mechanistic gain), but in
the p27 cycle the free ubiquitin node is also the network's only link from degradation back into
ubiquitination, and holding it cuts CCNE1's curated route. Net small and split (curator +9,
experimental −4). Not recommended as it stands. A version that holds free ubiquitin only where it is
not the sole link of a curated degradation cycle would need its own pre-registration.

## Decision

Adam, 2026-10-23: "lets not adopt then". `DS_UBIQUITIN_MODE` stays `propagate` (opt-in only).
