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
