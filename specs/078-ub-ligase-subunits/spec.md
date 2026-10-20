# specs/078: a ubiquitination does not deplete the ligase that carries it, re-measured under specs/077

**Status:** pre-registered 2026-10-20, before the build and the arms; amendment 1 (structure only, before any arm) gates `ubx`. **Gate fails**: held-out −44 (12 / 56), EGFR −29, RUNX2 −11; experimental +5 (all four CCNE1 cases). `ubc` the same. Not adopted; results in `research.md`.

Adam, 2026-10-19: "I would check the generator first", and "ok" to rebasing the specs/057
branch and running the arm against the new defaults.

## The defect (traced in specs/077)

Under `DS_CARRIER_LOOP_OPEN=fed_single` (default since specs/077), the SCF(Skp2):Cks1
ligase reads its outside producers, not its own release. Free CKS1B is the ligase's
supply, and it is drained by depletion edges from the generator's ubiquitin pass
(`_emit_substrate_depletion_edges`, LNG `5e663d9`):

- R-HSA-187575 "Ubiquitination of phospho-p27/p21". The catalyst is its own input,
  cyclin E/A:CDK2:p-p27/p21:SCF(Skp2):CKS1B. The output is the same complex
  with 3xubiquitin added.
- R-HSA-1363331 "Ubiquitination of p130 (RBL2) by SCF (Skp2)". The output includes the
  released ligase CUL1:SKP1:SKP2:CKS1B.

The pass takes every input leaf whose gene is also in the output as the substrate. So
CKS1B, SKP2, CUL1 and SKP1 (the ligase) are treated as being ubiquitinated, along with
p27/p21 and p130. In Mitotic G1 and S phase this decides 4 experimental CCNE1 cases
(specs/077 research). Before specs/077, the ligase's own release masked the drain.

## The rule (generator, unchanged from specs/057)

LNG branch `ub-substrate-078`: the specs/057 commits (LNG PR #113, closed) cherry-picked
cleanly onto LNG main `5e663d9`, the canonical build's generator. With both flags off it
is the generator that built canonical.

- `LNG_UB_SUBSTRATE_CHANGED=1`: a ubiquitination's substrate is an input leaf whose gene
  is in the output **and whose own stId is not**. A partner that leaves unchanged is
  carried, not depleted. This drops both CKS1B drains above.
- `LNG_UB_SUBSTRATE_UBX=1` (specs/057 amendment 1, read only with CHANGED): a leaf is
  still the substrate if it is a direct component of an output complex that has Ub as a
  direct component. This keeps the R-HSA-187575 drain (CKS1B sits inside the
  3xubiquitin complex) and drops the R-HSA-1363331 one.

## Arms

| arm | build | solver | against | gated |
|---|---|---|---|---|
| `ubc` | variant `--env LNG_UB_SUBSTRATE_CHANGED=1` | main `3e96e7c` defaults | `results/3e96e7c/dflt` (canonical) | **yes** |
| `ubx` | variant `--env LNG_UB_SUBSTRATE_CHANGED=1 --env LNG_UB_SUBSTRATE_UBX=1` | the same | the same | no (reported) |

`dflt` is main with no overrides on canonical `20261017-1246_5e663d9`. It reproduces
specs/077's `ocs` on every prediction (24,100 curator, 849 experimental).

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade --pair-on-valid`, because
a different network can change which cases are valid. Both figures are reported, and
cases are checked one by one.

Before either arm, the depletion edges each build removes are counted by diffing the
variant networks against canonical and committed.

**Gate (specs/077's):** curator held-out net ≥ 0, experimental net ≥ −5, and no
pathway below −10 (RAF set aside). If it passes, the flag is recommended as the
generator default and a rebuilt canonical is the adoption step. Adam decides.

## Predictions

1. `ubc`: of the 4 CCNE1 experimental cases specs/077 traced (Mitotic G1 OE at
   CCNA1/CCNA2, S phase KO at R-HSA-187568 and R-HSA-187920), at least 2 are fixed.
2. `ubc`: CBL, EPS15 and SH3GL1 lose again, as in specs/057, because that
   mechanism does not involve carrier releases. MET/EGFR curator net from them: −20 or worse.
3. No prediction is made for the held-out net. specs/057 measured −39 (CHANGED) and −15
   (UBX) before specs/077, and the history says the gate most likely fails.

## Not in scope

- The faithful build: its generator is on another branch (`diagram-subpathways`).
  Combining them is a later step if `ubc` is worth it.
- A rule for the CBL class (specs/057 R1, not registered there either).

## Structural count (`count_depletion.py`, after the builds, before any arm)

Builds `20261020-0717_9770d28_ubc` and `20261020-0756_9770d28_ubx` (LNG `9770d28`),
compared with canonical `20261017-1246_5e663d9` by entity id:

| build | depletion edges | removed | pathways | other edges |
|---|---|---|---|---|
| canonical | 727 | — | — | — |
| `ubc` | 420 | 307 | 32 | −356 `variant_pool` (depleter pools no longer needed), −56 `assembly` (EGFR) |
| `ubx` | 441 | 286 | 32 | −349 `variant_pool`, −56 `assembly` (EGFR) |

The 56 assembly edges are the knock-on specs/057 recorded: CBL-containing EGFR complexes
are no longer decomposed as roots. Nothing is added.

Mitotic G1 has three such depletion edges on canonical. Two drain free CKS1B (from
R-HSA-187565 and from R-HSA-1363326), and one drains p-RBL2 (from R-HSA-1363326), the
real p130 substrate.
- `ubc` removes all three. p-RBL2 is the substrate, but Ub-p-RBL2 keeps its stId, the
  case specs/057 amendment 1 fixed.
- `ubx` removes the two CKS1B drains and keeps the p-RBL2 one. **The spec above was wrong
  on this:** CKS1B is not a direct component of the ubiquitinated output R-HSA-187568,
  whose direct components are R-HSA-187565 and Ub. So `ubx` treats CKS1B as carried there too.

## Amendment 1 (2026-10-20, from the structural count only, before any arm)

`ubx` removes exactly the ligase drains this spec targets and keeps the real substrate's.
`ubc` also removes real substrates' drains. **The gated arm is `ubx`**, with the same
gate; `ubc` is reported. Predictions 1 and 2 apply to `ubx` as written for `ubc`.
