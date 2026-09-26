# 024 — Why genes lost their root form, and a minimal pull towards 1 in loops

## 1. The genes MP-BioPath could perturb and we cannot (measured 2026-09-25)

MP-BioPath maps a gene to its own protein, DNA and RNA nodes
(`geneToRootNodes`, mp-biopath/src/idMap.jl), and perturbs those that are ROOT
in its network, i.e. have no parents (pi.jl). Across the 877 benchmark
perturbations:

| MP-BioPath network | ours (root pins, specs/023) | perturbations |
|---|---|---|
| root form | root form | 765 |
| **root form** | **no root form** | **90** |
| gene absent | root form | 10 |
| gene absent | no root form | 8 |
| in network, no root | root form | 4 |

Why the 90 have no root form in ours:
- **29:** the gene appears only inside complexes.
- **20:** its own entity node is produced by a reaction.
- **27:** it is fed partly or wholly by our derived `dissociation` / `depletion`
  edges.
- **14:** the gene is not in our network (version skew, or resolution; not yet
  split).

**Traced, and it is the same mechanism twice: a catalytic cycle.**
- **ALKBH2:Fe2+** binds damaged DNA, demethylates it and is released
  unchanged. It is the output of its own reaction, so it is fed only by its
  own regeneration.
- **CBLB**, an E3 ligase, is released by *"Release of E3 from
  polyubiquitinated substrate"* and fed back the same way.

MP-BioPath's published networks had loops removed by hand, which cut these
cycles and left the free enzyme as a root. Our generator reproduces Reactome's
cycle, so nothing is a root. This is the catalytic-recycling loop type from the
2026-06-12 loop taxonomy. **It is how our generation differs from the published
networks, not a change in Reactome.**

**The rule (`DS_PIN_SCOPE=root_cycle`, bench side).** For a gene with no true
root:
- pin the resolved nodes that are fed only by their own downstream (a cycle
  regenerating them) or by derived dissociation/depletion edges;
- prefer the gene's own entity node;
- otherwise take the entry complexes, or the whole cycle's pool if they
  regenerate each other.

Genes with a true root are untouched. Measured offline, it recovers **75 of the
76** in-network cases and changes the pins of **0** of the other 780
perturbations.

## 2. A minimal pull towards 1 in loops (Adam, 2026-09-25)

Adam's answer on question 1:
- The pull should be *"extremely minimal for loops"*.
- A cycle **without** inhibitors should show *"feedback amplification of the
  values when they are above one and also be pulled away from one in the
  negative direction if input values are below one"*.
- A loop **with** negative interactions *"should act as a dampener, so pulling
  values towards 1 is less of a big deal as long as it is extremely minimal"*.

Why a minimal pull delivers exactly that for positive loops:
- With every in-loop edge at elasticity 1, a positive cycle's gain is exactly
  1. That is the specs/014 knife-edge: any leak rails the loop to 100x or 0, and
  which one depends on sweep order.
- Reading loop edges at fold^ε with ε slightly below 1 makes the loop gain g
  slightly below 1. The loop then *amplifies* a sustained input by a finite
  factor (≈ 1/(1 − g) in log terms) instead of railing. Up stays up and grows,
  and down stays down and deepens, which is the behaviour Adam describes.
- The pull's strength sets the amplification.
- For loops containing an inhibition, the same small pull only adds to damping
  the loop already has.

The existing knob is `DS_LOOP_ELASTICITY`: constant elasticity on activator
edges whose source and target share a component (width 0).

## Pre-registration (committed before any arm runs)

All arms run through `scripts/run_arm.sh` on build `20260925-1039_d4f4f64`,
against the current defaults (`results/5f7b9c7/new_defaults`).

- **Arm P, `--bench DS_PIN_SCOPE=root_cycle`:** a protocol refinement, the
  faithful reproduction of MP-BioPath's roots.
  - Paired results on shared cases must be identical (only newly valid cases
    change), otherwise the arm is not interpreted.
  - Reported: how many cases become valid, and their accuracy.
  - **Adopt as the default** if the newly valid cases score at least as well as
    the old every-case figure scored them as NORMAL. It is a protocol decision,
    so it goes to Adam with the numbers.
- **Arms L1 and L2, `--server DS_LOOP_ELASTICITY=0.99` and `0.95`.**
  - Decided on curator held-out: net, McNemar, pathways moved, perturbations;
    experimental reported.
  - **Adopt** only if held-out net > +15 with p < 0.05, the gain is not one
    pathway, and experimental is no worse than −15.
  - If both pass, adopt 0.99 unless 0.95 is better by more than 15 held-out,
    because Adam asked for minimal.
  - Report the positive-loop and has-inhibition loop classes separately (MH,
    `failure_structure.py`).

## Found while splitting the 14 "gene not in our network" cases: one catalog pathway is the wrong pathway

`bench/catalog_pathways.tsv` builds "Interleukin-2_family_signaling" from
**R-HSA-447115**, taken from MP-BioPath's `pathway_list.tsv`. In Reactome v97,
R-HSA-447115 is **"Interleukin-12 family signaling"**; Interleukin-2 family
signaling is **R-HSA-451927**, where IL2 sits in 18 reactions and JAK3 in 44.
So the benchmark has been scoring **260 curator cases** of IL-2 family biology
against the IL-12 family network. That accounts for 7 of the 14 genes missing
from our networks (IL2, IL21R, JAK3, LCK, PIK3CA, STAT5B, SYK).

The other two ids whose v97 name differs are renames of the same pathway:
- R-HSA-453279, Mitotic G1 phase and G1/S transition;
- R-HSA-388841, Regulation of T cell activation by CD28 family.

**Fix, not yet applied:** point the entry at R-HSA-451927 and rebuild the
catalog. That is a new build (`scripts/catalog.sh build`), so it is left for a
deliberate step. Of the remaining 7 missing genes, only COL1A1 and COL2A1
(GPVI, immunoregulatory, PDGF), EP300 (NER), GDI1 (Rho GTPase cycle) and FBXW7
(RUNX2) show 0 v97 reactions in their pathway. Those are version skew, which no
generator change can recover.

## Results (build `20260925-1039_d4f4f64`, arms at `5208139`, base `5f7b9c7/new_defaults`)

**Arm P, `DS_PIN_SCOPE=root_cycle`.**
- Shared cases are identical (0 / 0), as required.
- Pins: 1,246 over 854 perturbations, up from 1,023 over 779.

| newly scorable | n | correct | as scored today (NORMAL) | under the old broad pins |
|---|---|---|---|---|
| curator held-out | 726 | **81.7%** | 32.0% | 87.7% |
| curator all | 1,052 | 74.5% | 31.6% | 81.9% |
| experimental | 32 | 50.0% | 0% | 71.9% |

Every-case figures:
- curator held-out 82.89% → **84.79%**, macro-F1 0.7678 → **0.7979**;
- curator all 81.00% → 82.88%;
- experimental 64.55% → 66.43%.

It passes the pre-registered bar: the newly valid cases score far better than
NORMAL. As registered, adopting it as the default protocol is **Adam's
decision**, and it is not flipped here.

**Arms L1 and L2, `DS_LOOP_ELASTICITY` = 0.99 and 0.95. Not adopted.**

| ε | held-out net | pathways / perturbations | without Interferon-γ | tuning | experimental |
|---|---|---|---|---|---|
| 0.99 | +30 (41 / 11), p 4e-5 | 3 / 4 | −4 | −17 (TP53, RUNX) | −1 |
| 0.95 | +48 (63 / 15), p 4e-8 | 10 / 16 | +14 | −72 (TP53 −57, PIP3 −14) | −2 |

Both fail the concentration condition: Interferon-γ is +34 in each. The
tuning split loses, mostly TP53, whose giant component the pull affects most.
A constant elasticity on every in-loop activator edge compounds around a long
cycle (ε^L), so it is **not** the "extremely minimal" pull the design asks for
on large loops. A pull that is minimal per *loop* rather than per *edge*, e.g.
ε applied once at the loop's entry, is the obvious next design. Not attempted
here.
