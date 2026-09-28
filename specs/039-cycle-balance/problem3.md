# specs/039 problem 3: where does a pool end?

For two blind derivations, as for `problem.md` and `problem2.md`. Read
`research.md` (the results, amendments 1–4 and the traced losses) first.

## What exists

The pool of protein R is built as follows:
- It is a strongly connected set of R-forms joined by R-steps (reactions with
  exactly one R-containing input and one R-containing output).
- Its **states** are the core-only forms, one per modification signature.
- Its **intermediates** are the enzyme complexes on the paths between states.
  They hold zero baseline share and read v/k_cat.
- It is solved as π Q = 0 with drives evaluated with the source held at
  baseline.
- Supply comes from outside the pool; carriers are the enzyme's free forms.

Measured, and not adopted on either axis. The losses trace to one class:

1. **Interferon α/β, −56.** SOCS binds the receptor's base state (A + SOCS →
   A:SOCS). A:SOCS divide-inhibits the pool's own forward step. A:SOCS does not
   return to the pool; it leads to degradation. So it is not in the strongly
   connected set, and it is not an intermediate. The inhibitor tracks A
   exactly, and the forward flux is pinned at 1.000 whatever happens upstream.
2. **HDR, −63.** The pool's exit step has co-inputs made from the pool's own
   intermediates (p-CHEK1 via ATR). With the pool bounded, the zero root is the
   only attractor. Holding those inputs at baseline (amendment 4, post hoc)
   fixed HDR (+51).
3. The same hold did **not** fix α/β, which was left at −56. It also **removed**
   Interferon γ's +34 gain. There, SOCS is built from the p-JAK2 state and
   inhibits the SECOND transition, and conservation had redistributed mass
   correctly.

## Questions

1. **R-forms that leave the chain and do not return.** A:SOCS (sequestration
   then degradation), dissociation sinks, and degradation are examples. Are they
   states, sink states (outflow), or outside the pool?
   - If they are sink states, what is the steady state? Without a return, π
     leaks, and supply must balance the outflow.
   - What do they read?
2. **A complex built from a pool form that regulates the pool's own step.**
   When is its effect a double count, as with the pre-registered catalyst ⊣
   source-state edge, and when is it genuine negative feedback that must act?
   - The α/β and γ cases give opposite answers under the same hold. State a
     rule that distinguishes them, or show that none exists at this level.
3. **Inputs made from intermediates** (the HDR case). Are they a detection
   error (the pool is a machine, not an enzyme cycle), or something the solve
   should handle?
4. Say which of these are detectable from Reactome by the generator, and how.
5. Work the three cases numerically: α/β JAK1 80x (expected UP), γ IFNG 80x
   (expected UP), and HDR RAD51 80x (expected UP). Also state what the rule
   predicts for a SOCS1 KO in each.
6. **What can go wrong.** The fix must not be tuned to these three pathways.
   Say what it predicts elsewhere, and how to tell in advance whether it
   generalises.

Keep it to what a pre-registration can state, and give numbers.
