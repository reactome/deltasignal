# specs/035 — Conserved loop participants held at baseline (design D1)

**Status:** pre-registered 2026-09-26, before any code or arm.

## Why (specs/034 §6)

- **7,989 reaction-input edges** (9% of all) in 70 pathways are supplied ONLY by
  the reaction's own downstream: every producer of the input node is reached
  from the reaction.
- 87% have no other usable copy of the entity, so nothing outside the loop
  anchors its baseline. A loop whose AND inputs are all loop-carried has a gain
  of 1 at baseline (the knife-edge, specs/013 and 014), and a push can drain it
  to the all-zero state.
- **RAF/MAP kinase is traced** (specs/032, 033): the 128-node component still
  collapses to 0 under every RAS/RAF perturbation after pooling. Experimental
  is 17 of 49, where MP-BioPath gets 46.

**The nearest failed attempt, and the difference:**
- `DS_SCC_BREAK_CATALYST` froze every recycling catalyst at its entry value:
  TP53 +104, held-out −65, not adopted.
- It froze catalysts that also carry real signal.
- This rule is narrower. It touches only entities **released unchanged** by the
  loop: every producer of the node consumes a complex that CONTAINS the node's
  entity (or is a `dissociation` edge from one). That is a conserved pool,
  bound and released, like a cofactor.
- Entities transformed inside the loop (RAS:GTP, the RAF complexes) are not
  touched.

**Split of the 7,989:**

| class | edges | where |
|---|---|---|
| released unchanged (conserved) | 1,722 (22%) | 36 pathways: Class I MHC 499, DNA Damage Bypass 150, PTK6 150, FGFR2 140, ERBB2 119, ERBB4 112, HDR 81, RAF 67 |
| transformed in the loop | 6,267 (78%) | not touched |

RAF's conserved set:
- the RAS GAPs: NF1, RASA1–4, RASAL1–3, SYNGAP1, DAB2IP;
- the scaffold partners: F-actin, CNKSR1/2, IQGAP1, KSR2, TLN1, ARRB1/2,
  RAP1A/B, APBB1IP;
- the ions.

## The rule

`DS_CONSERVED_MODE` = `off` (default, byte-identical) or `inert`.

A node is **conserved** if all four hold:
1. it is an input of a reaction in a non-trivial strongly connected component;
2. every non-depletion producer of it is reached from that reaction;
3. every producer is a release: a reaction whose input or catalyst CONTAINS the
   node's entity (by the bundle's containment table, or the same stId), or a
   `dissociation` edge from such a complex;
4. it has at least one producer.

Under `inert`, conserved nodes are pinned at baseline exactly as cofactors are:
participants at fold 1.0 that never carry a perturbation. An explicit pinning
observation (confidence above the gate) wins, so a perturbed GAP such as NF1
is still perturbed. The solve reports `conserved_held`, and the benchmark
refuses an `inert` arm that held nothing.

## Pre-registration

**Arms:**
- `consctrl`: code defaults;
- `consinert`: `DS_CONSERVED_MODE=inert`;
- both on the canonical build `20260926-1221_590301c`, through
  `scripts/run_arm.sh`.

**Adopt only if all hold:**
- curator held-out net > +15 with p < 0.05;
- experimental net ≥ −15;
- no single pathway loses more than 10 held-out cases;
- the gain is not a single-pathway effect: at least 3 pathways moved, or the
  result is reported as such.

**Predictions:**
- RAF experimental rises (at least +10 of 49);
- RAF convergence improves;
- Class I MHC, DNA Damage Bypass, PTK6 and the FGFR/ERBB pathways move;
- pathways with no conserved nodes are byte-identical.

**Also reported:**
- tuning curator;
- per-pathway net on both axes;
- concentration;
- convergence counts;
- `conserved_held`.

## Result — NOT ADOPTED (2026-09-26)

**Setup:** canonical build `20260926-1221_590301c`; arms `consctrl` and
`consinert` at solver 1fcc1bc. The inert arm held 12,070 curator and 1,256
experimental node-solves. Convergence: curator 1,652 → 1,658, experimental
217 → 226 of 244.

| split | net | fixed / broke | p | pathways |
|---|---|---|---|---|
| curator held-out | **+42** | 110 / 68 | 0.002 | 13 (+9 / −3) |
| curator tuning | −265 | 100 / 365 | 1.7e-36 | PIP3 −308 |
| **experimental** | **−122** | 3 / 125 | 2.1e-33 | 3 (all negative) |

- Best: HDR +54, ERBB2 +26, ERBB4 +18, DNA Damage Bypass +17.
- Worst: PIP3 −308, Mitotic G2/M −42, **TGF-β −30**, **Class I MHC −18**.

**Gates:**
- held-out passes (+42, p = 0.002);
- **experimental fails** (−122);
- **no-pathway-below −10 fails** (TGF-β −30 and Class I MHC −18, both held-out).

**Not adopted.** `DS_CONSERVED_MODE` stays in the code as a record, default
`off`.

**Why PIP3 collapsed:** PIP3 itself (R-ALL-179838) is one of the two nodes held
in that pathway.
- That PIP3 node is produced only by reactions that take PIP3 in and give it
  back (it binds its effectors and is released), so it meets "released
  unchanged".
- But its level IS the pathway's signal, a second messenger.
- **"Released unchanged" does not separate a conserved scaffold (F-actin, a RAS
  GAP) from a signalling molecule that is bound and released.**
- It is also copy multiplicity again: the loop-only PIP3 copy is not the one
  PI3K produces.

**What survives:**
- The held-out gain is real: +42, 13 pathways, 66 perturbations.
- It sits in the pathways the inventory predicted: HDR, ERBB2/4 and DNA Damage
  Bypass, the loop-only-input pathways of specs/034 §6.
- A narrower rule (exclude small molecules and second messengers, or require a
  protein or complex entity) would be a NEW pre-registration. It is not a tuning
  of this one. Not run.

**Scored per pathway (experimental, ctrl → inert):**

| pathway | ctrl | inert |
|---|---|---|
| **RAF** | 14 / 49 | **14 / 49** |
| PIP3 | 168 / 200 | **49 / 200** |
| Mitotic G1 | 35 / 88 | 35 / 88 |
| ERBB2 | 42 / 49 | 41 / 49 |

- **The whole experimental loss is PIP3. RAF, the motivating case, did not
  move.** Holding its 29 conserved entities at baseline does not stop the
  collapse either.
- Three attempts have now failed to fix RAF: inert drugs (specs/032), set
  pools (specs/033) and conserved participants (here). Each targets a
  structural suspect found by reading the graph.
- **The next step is to trace the RAF solve's iterations directly** (which
  node first leaves baseline, in which direction, and why the loop does not
  settle) before proposing another rule.
