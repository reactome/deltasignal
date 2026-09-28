# problem3: derived independently (Opus), before reading the Fable derivation

## The principle

Amendments 1–4 each removed one way the pool's protein acted on the pool twice:
- catalyst ⊣ source state;
- the pool depleting its own supply;
- pool-fed inputs of its own steps.

They are one principle, which the losses show was applied only in pieces:

> **The pool's own protein acts on the pool only as mass.** Every R-form the
> pool can reach is part of the pool's mass balance, as a state, an
> intermediate or a *sink*. An edge from any R-form onto the pool's own
> transitions is not applied as an edge, because its effect is already the
> mass that form holds.

What fixed HDR and failed on Interferon α/β and γ was holding *inputs* at
baseline without adding their mass. The sequestration then vanished instead of
being moved into the balance.

## Q1: forms that leave and do not return are sink states

- A:SOCS, degradation-bound forms and dissociation sinks become **sink states**:
  R-forms reached from a state by an R-step (as defined), with no R-step back
  into the chain.
- The chain is then open. Supply σ (the supply fold × the baseline inflow)
  enters at the states fed from outside, and sinks drain.
- The steady state is x = σ·(−Q̃)⁻¹, where Q̃ is Q restricted to the non-sink
  states, with sink exits as killing rates.
- **Baseline:** at rest, total inflow equals total sink outflow (turnover), so
  no new parameter is needed:
  - The sink rates' baseline split among several sinks is by exit weight, as
    for transitions (enzyme-driven 1, intrinsic 1e-3).
  - φ₀ still sets the state split.
  - A pool with no sink stays closed, and reduces to amendment 3.
- A sink state reads its inflow flux relative to baseline, as an intermediate
  reads v/k_cat. Downstream, the degradation readout reads that flux.
- **Consequence:** the pool total is no longer fixed. It is supply ÷ outflow.
  Knocking out the only sink's driver (SOCS1 KO) sends the total to the cap:
  UP, bounded by 100.

## Q2: which regulating complexes are double counts

- **Rule:** any inhibitor or activator edge whose SOURCE is an R-form of the
  pool (state, intermediate or sink), acting on one of the pool's own steps,
  is not applied. Its effect is the mass it holds: SOCS-bound receptor cannot
  be activated because it is not in the base state.
- Edges from **non-R** nodes act as they do now. So SOCS itself, as the free
  protein, drives the sink step by mass action, and SOCS KO removes the leak.
  That is the genuine negative feedback, and it survives.
- This distinguishes nothing between α/β and γ, and does not need to: under
  this rule both give UP.
  - **α/β:** A:SOCS is a sink off the base state, and its inhibition of
    "Activation of JAK kinases" is dropped. pp ≈ s ≈ 1.25 under JAK1 80x.
    - That is UP, but only just: s is capped by the pathway's USP18-type
      inhibition of supply, not by the pool.
    - A **stated risk**: 1.25 is close to the 1.15 cutoff.
  - **γ:** the SOCS complex from p-JAK2 is a sink off p-JAK2, and its
    inhibition of the second transition is dropped. pp ≈ s = 80 · (leak
    correction) → UP.
- **Amendment 4's loss of γ's +34 is unexplained**, and must be traced before
  this is pre-registered. By this reasoning, holding the inhibitor should have
  left γ UP. Something else in the held set, such as a supply or readout node,
  is suspected.

## Q3: HDR is a detection error

- A path whose step co-input is reachable by mass flow from the pool's own
  intermediates is a **machine**, not an enzyme cycle. Such paths are dropped,
  and a pool left without a path in either direction is not a pool. Detectable
  in the generator by reachability in the network, and counted.
- This is cleaner than amendment 4's hold, which kept the machine and patched
  its inputs.
- **Prediction:** HDR returns to the `off` behaviour, so its +51 is kept.

## Q4: detection

- Sinks: R-steps out of the chain with no return. The generator already has the
  R-graph; the sink states are the nodes one step outside the strongly
  connected set.
- Self-regulating edges: edges whose source is any pool form (state,
  intermediate or sink) onto a pool step node. Listable at generation time.
- Machine paths: as in Q3.

## Q5: numbers (φ₀ = 0.1)

| case | prediction |
|---|---|
| α/β JAK1 80x | pp ≈ 1.25, UP (marginal) |
| α/β SOCS1 KO | total → cap, UP |
| γ IFNG 80x | pp ≈ 80 × (1 − baseline leak share, restored), UP |
| γ SOCS1 KO | UP |
| HDR RAD51 80x | = off, 100, UP |
| HDR any KO | = off |

## Q6: generalisation, checked before measuring

- Before any arm, list every pool with a sink state or a self-regulating edge,
  catalog-wide, and the curator and experimental cases whose readouts lie
  downstream of them.
- **Pre-register the gate on the cases OUTSIDE the three motivating pathways**
  (α/β, γ, HDR). If that subset is too small to test (say fewer than 100 cases
  in fewer than 3 pathways), say so in advance, and treat the result as
  descriptive.
- **Risk:** the open chain makes pool totals supply-sensitive again. That is
  right for turnover-driven proteins, but it reintroduces magnitudes that
  conservation had bounded. A pool with a large sink share can run to the cap.
  Report how many pools reach the cap in the arm.
