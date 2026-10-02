# specs/042: derived independently (Opus), before reading the Fable derivation

## 1. The rule: extend an adopted mechanism instead of adding a new one

The generator already derives `depletion` edges (catalyst → substrate) for two
enzyme classes: phosphatases (Pi as an output) and ubiquitin ligases (Ub as an
input).
- These are **removal** reactions: they take the substrate's active form away
  and do not give it back.
- They are load-bearing (removing them is held-out −90) and bounded (specs/011:
  a fold between 0.1 and 10).

**Rule (R1): every non-returning removal of X gives its drivers a depletion
edge onto X.**
- **Drivers:** the reaction's catalysts, and its non-X inputs that are proteins
  or complexes (a binder, a trap, a destruction complex).
- The two classes above are special cases. Nothing else changes: the same
  divide form, the same 0.1–10 bound, the same solver code.

**Mass balance.** At steady state, X = supply / outflow. If removal is X's
dominant exit, a remover with drive u gives X ≈ s / u, which is exactly what a
divide-form depletion edge computes.
- "Removal is dominant at steady state" is the standard turnover assumption:
  everything synthesised is eventually degraded.
- So **no new parameter** is needed. The weights would only matter if
  non-removal exits shared the outflow, and R1 deliberately does not let them
  act (see §2).

## 2. What counts as non-returning, and why the blanket rules failed

An exit reaction r (X is an *input* of r, not a catalyst) **removes** X if
both hold:
- **(a)** no output of r contains X's protein, i.e. degradation; or every
  output containing X's protein is a **dead end** in the pathway: nothing
  consumes it, or it is consumed only by reactions that are themselves
  removals (recursively, depth ≤ 4).
- **(b)** r is not a step of a specs/039 pool. Interconversion is conserved
  there and already handled.

**Exits that carry X onward** (β-catenin → nucleus → TCF; a complex that goes
on to signal) are **not** removals. They do not get a depletion edge.

That is the difference from the failures:
- **The blanket "general consumption" rule** made every consumption deplete the
  co-input, including signalling steps whose product carries the signal. Every
  signalling cascade then fought itself, giving over-coupling, the dominant
  error class (48% false change).
- **Naive substrate depletion** added catalyst → substrate on every catalysed
  reaction. That is the same mistake for catalysts.
- **specs/019 sink bridges** connected sinks to *other* nodes, adding positive
  routes. R1 adds only negative edges, from removers to what they remove.
- **specs/035** held conserved loop inputs at baseline: it removed signal. R1
  adds routes.

## 3. Composition with existing rules

- **specs/039 pools:** pool forms are excluded as targets. A pool's own exits
  are its transitions. A degradation exit from a pool state could feed the
  pool's supply, but that is deferred.
- **specs/022/040 self-contained inhibitors:** the depleter must be the free
  remover, never a complex that contains X. So the destruction complex
  *bound to β-catenin* is not the depleter; the free APC:AXIN:GSK3 complex is.
  That avoids a self-tracking inhibitor.
- **specs/011 bound:** unchanged, at most 10x either way.

## 4. Detection and census (generator)

- **Detection:** reuse the pool code's R-graph (one reference protein through
  its R-steps). Mark exits that meet (a) and are not (b). Emit a `depletion`
  edge from each driver to the X node, and tag it `derived_removal`, so it can
  be switched off as an arm (`LNG_REMOVAL_DEPLETION=1`, default off).
- **Census before any arm:**
  - removal exits, new edges and affected pathways;
  - overlap with the existing phosphatase and ligase depletion edges, which
    must be a superset of them (a check that the rule reproduces adopted
    behaviour);
  - **agreement with MP-BioPath's hand brakes** (specs/034 §12d). Those are
    independent expert judgements, used only as a validation set: precision
    and recall of R1's edges against their 219 additions. This is a check,
    not a fit.

## 5. The cases, worked

- **WNT, APC.** Free β-catenin's exit into the destruction complex leads to
  p-β-catenin → Ub → degradation, which qualifies under (a).
  - New edge: free destruction complex ⊣ free β-catenin.
  - APC KO: the remover's drive is 0, so β-catenin is de-repressed to the
    10x bound: **UP**.
  - APC 80x: β-catenin goes down to the 0.1 bound: **DOWN**.
  - Both currently inverted.
- **p27 (CDKN1B).** p27 + cyclin:CDK2 → p27:cyclin:CDK2, an inactive complex
  that leads to p27 degradation and dead ends; this qualifies if it is a dead
  end in S phase.
  - Edge: p27 ⊣ cyclin:CDK2. p27 KO: CDK2 **UP** (bounded).
- **RB1.** RB1 + E2F → RB1:E2F, a dead end (repression), so it qualifies.
  - Edge: RB1 ⊣ E2F. RB1 KO: E2F targets **UP**.
- **HDR branch competition (HR vs SSA).** Both branches carry the intermediate
  onward, so neither is a removal. **R1 does not cover it**, and the HDR cases
  stay wrong. Branch competition is all-exits mass balance, the family that
  failed as a blanket rule. It is left out on purpose.

## 6. What can go wrong

- **Incompletely curated dead ends.** Many Reactome complexes have no curated
  downstream because curation stopped, not because they are biologically
  inert. R1 would turn every such stop into a depletion, which is over-coupling
  again.
  - **Mitigation (pre-registered as a variant):** R1-strict counts a dead end
    only if it is degradation, or a complex whose name or GO annotation marks it
    inactive or inhibitory, or a curated negative regulation of something.
  - Measure R1 and R1-strict.
- **False change rises.** Report the false-change rate per arm. Gate: the
  false-change count must not rise by more than the fixed-case count.
- **Prediction outside the motivating pathways:** binders and traps
  (FST/Activin, IL13RA2, LEFTY, SOCS6/KIT, CDK6/RUNX1, BCL2/BAK) gain routes,
  which matches the ~190 no-route cases in §12d. If R1 helps there, on
  pathways it was not derived from, that is the evidence it generalises.

## 7. Gates

The same as specs/040:
- curator held-out > +15, p < 0.05;
- experimental ≥ 0 and not significantly negative;
- no pathway loses > 10;
- ≥ 2 pathways gain;
- net ≥ 0 outside WNT, HDR, S phase and G1/S.

Plus the false-change check, and agreement with the specs/041 expected-biology
review where its rows apply.
