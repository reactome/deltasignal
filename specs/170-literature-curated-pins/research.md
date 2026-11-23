# Research: literature edges with curated pins (specs/170)

Run `lit170` at `be1fbe08` (registered at `be85b322`, scheduler `be1fbe08`, before any arm), every arm with
`--bench DS_SOLVE_TIMEOUT=3600 --bench DS_PIN_CURATED_ROOT=<builds>/20261120-lit-v3-without`. Each arm log's
`Protocol:` line prints the flag. Results: `builds/<arm>/results/be1fbe08/lit170`; scores:
`~/deltasignal-catalogs/analysis/170/score_170.json` (`score_170.py <builds> 91189f0b be1fbe08 <out>`).

## Verdict: the gate FAILS (net −554). The edges do not help under curated pins either.

| check | result |
|---|---|
| FR-002 control parity (`v3-without` lit170 vs lit168) | identical on every axis (0 class changes, 0 cases scored in one run only) |
| FR-003 pin check | v3-single: a literature edge unroots 266 curated roots without the rule; 0 differences under it. v3: 316; 0 |
| gate: v3-single-with vs control | **primary net −554** (korkut −470: 17 fixed / 487 broken; experimental −84: 19 / 103) |
| shuffles | not run (net ≤ 0, as registered) |

The prediction was half right. korkut failed, as predicted. The experimental loss did not shrink, and ERBB2
over-expression in PIP3 does not read UP again. Both parts of the prediction rested on the specs/168 reading that
unrooting caused the PIP3 breaks, and the trace below shows that reading was wrong.

## Why it is worse than specs/168 (−554 against −374)

The pairing is on the cases both arms score. Unrooting in specs/168 did not only break cases: it made some invalid,
because a gene with no pin is not perturbed. Under curated pins the gate arm scores 874 more korkut cases and 32 more
experimental cases. On those it is −84 and −12.

On the cases both runs score (korkut 3,173, experimental 513):

| axis | specs/168 | specs/170 |
|---|---|---|
| korkut | −305 (75 / 380) | −386 (15 / 401) |
| experimental | −69 (21 / 90) | −72 (19 / 91) |

- **The specs/168 korkut "fixes" were unpinned drug targets.** 58 of the 75 cases it fixed read exactly 1.0 in specs/168:
  the target had lost its pin, so nothing moved and NOT was right. Most are SRC-inhibitor drug sets (SR in 51 of
  them), which specs/138 found did not engage in these cells. With the pin back, 62 of the 75 read DOWN and are wrong
  again.
- **Experimental barely moves.** 102 of the 103 experimental breaks are in PIP3, and they are the same cases that
  broke in specs/168 (84 of them read DOWN there too).

Protocol effect (v3-single, lit168 → lit170, paired): korkut −81, experimental −3, curator held-out −10, tuning +60.

## The traced case: PIP3 over-expressions read 0 because four literature edges inhibit MTOR

Traced on trace servers at `be1fbe08` with the bench's own `pin_occurrences` under the flag. The control reproduced all
93 recorded PIP3 over-expression folds of the gate arm. The scripts and edited pathway copies are in
`~/deltasignal-catalogs/analysis/170/trace/`.

- **Pins are not the cause.** ERBB2 is pinned in both arms, and its AKT readouts still read about 0.
- **The edges.** Four literature `regulator`/`neg` edges, each from one paper (prior confidence 0.475), point into
  MTOR, which is a curated root:

  | edge | key |
  |---|---|
  | PIK3CA ⊣ MTOR | `01deb7ece2bc701d` |
  | PIK3CB ⊣ MTOR | `0dd658b0a015a8cc` |
  | PIK3CD ⊣ MTOR | `333fa7d9fa92ff07` |
  | PIK3CG ⊣ MTOR | `1dd84fb029b57dfc` |

  The PI3K subunits are curated roots too. Literature in-edges, including a literature reaction that feeds itself
  among them, make them computed. Any over-expression upstream of AKT raises them about 20×.
- **The collapse.** With the subunits up about 20×, MTOR falls to about 3e-9 (internal scale). mTORC2, mTORC2:PIP3 and
  p-T309,S474-AKT2 go to about 0, and with them every AKT-target readout.
- **What holds it down.** A second literature loop, MTOR ⊣ AKT1S1 ⊣ MTOR (`d5e0bbbf855fd1c2`, `2fae6a06d5615627`),
  keeps MTOR low.
- **A switch, not a numerical failure.** Every solve converges (AKT1 OE residual 5.6e-7; the LM fallback never runs).
  At rest the two networks agree on all 4,532 curated nodes. But the switch detector (specs/154) finds a 1,438-node
  component that is bistable at rest, and its other state is exactly this collapse (MTOR 1e-9, AKT1S1 100×).
  Over-expression tips the solve into that state.
- **Causal check.** Remove only the four MTOR edges from a copy of the pathway:
  - AKT1 OE: readout 111910 goes 0 → 100, 198373 goes 0 → 50, 9614997 goes 0.0087 → 100.
  - Of the 78 cases that read UP without the literature edges and not UP with them (13 genes), 75 are restored.
  - The other 3 (readout 200163 under AKT1, AKT2 and ERBB3 OE) are a separate mechanism, not traced.
  - Removing the matching four PI3K ⊣ MLST8 edges as well changes nothing further.
- **Caveat.** Even without those four edges, the component is still bistable at rest. Its other state is then high
  (111910 at 100×), so the restored 100s are the high rail, not a graded signal.

**Correction to specs/168 research.md (reading item 2).** Unrooting ERBB2 (PHF8 → ERBB2) did not cause the PIP3
over-expression breaks. They come from literature inhibition of MTOR that tips a bistable component into its low state.
The "49 of 90 have the perturbed gene's starting form unrooted" was a correlation, not the mechanism.

## Described pairs

| pair | korkut | experimental | curator held-out (already reached: fixed / broken) | tuning (already reached) |
|---|---|---|---|---|
| gate (v3-single vs control) | −470 | −84 | 33 / 103 | 105 / 286 |
| v3 (every copy) vs control | −606 | −12 | 39 / 86 | 98 / 130 |

- **Readouts only the literature reaches** (reported as uninformative under the 2026-11-20 rule): the curator
  change-call precision is 0.027, and korkut's is 0.015 (403 calls) in the gate and 0.011 (530) in v3.
- **v3 holds up better on experimental** (−12) because every-copy attachment does not tip the PIP3 switch for most
  over-expressions, as in specs/168.

## What this means for the literature bundles

1. The pin question is now settled. Even with every curated gene pinned exactly as in the curated network, the edges
   lose on korkut (newly reached readouts are almost all over-calls) and on experimental (PIP3).
2. The PIP3 loss is four edges and one loop. Edges that make a curated root computed can also create alternative states
   in a curated loop. A sign-checked edge that inhibits a curated activator (PI3K ⊣ MTOR runs against the curated
   PI3K → AKT → mTORC1 direction) should at least be flagged at bundle time.
3. Option B (specs/171) is still a protocol question on its own. This result does not argue for or against it.
