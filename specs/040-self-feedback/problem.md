# specs/040 problem: a node's own downstream multiplied back into it

Given to two derivations independently (blind), as for specs/039.

## Context

DeltaSignal propagates fold-changes (1 = normal, cap 100) through a
Reactome-derived logic network:
- a reaction is the AND of its inputs (product, hill_sat), divided by its
  inhibitors;
- an entity is the OR-mean of its producers;
- a set-valued catalyst or regulator is one **set-pool** node that multiplies
  its members (`DS_SET_POOL_MODE=product`, specs/033/038).

Loops are solved per strongly connected component by fixed-point iteration.
Curated interconversion cycles are solved as conserved pools (specs/039, the
default). Its principles are that a pool cannot count its own output as supply,
and that the pool's protein acts on its own transitions only through its
states.

## The defect, traced (RAF/MAP kinase cascade, R-HSA-5673001)

p21 RAS:GTP responds correctly (KRAS 80x → 40x). The activated RAF dimer
complex one step below reads about 0 under every perturbation, and everything
downstream is dead. RAF gets 18 of 49 experimental cases right. Four
amplifiers stack at the RAF activation step:

1. **Recycled set-pool members.**
   - "Phosphorylation of RAF" is AND(ATP, RAS:GTP:RAF dimer, catalyst set-pool
     "RAF activating kinases").
   - Three of the pool's nine members are p-MEK dimers, produced only by
     "Dissociation of RAS:RAF complex", downstream of the activated RAF dimer
     itself. That is faithful curation (MEK phosphorylates RAF).
   - Under `product` the loop gain becomes fold³.
2. **A capped bundle's leaves fed by the step's own downstream.**
   - The next step, "MAP2Ks and MAPKs bind to the activated RAF complex", is a
     variant-capped bundle (generator cap, 12,000 > 512 variants).
   - Its 28 AND inputs include 22 scaffold leaves (F-actin, TLN1, ARRB1/2, KSR2,
     Ca2+ …), each fed ONLY by the dissociation downstream.
   - A dip of 0.47 is raised to the 22nd power each pass.
3. **A self-built inhibitor.** That step is divide-inhibited by RAS:GTP:activated
   RAF1 dimer:PEBP1. The inhibitor is built by AND from RAS:GTP and RAF1, the
   same species the step consumes, so it tracks the step's input and cancels
   it. specs/022's self-inhibitor rule misses it: that rule matches whole-node
   stId containment, and here the two share leaves, not a node.
4. **Drug-bound complexes** inhibit the MEK/MAPK phosphorylation steps
   (specs/032). This is **out of scope**: "drugs inert" was measured and
   refuted (experimental −5).

Every value is 1.0 at rest; the collapse happens only under perturbation. The
iteration shows the activated dimer rising (4 → 12 → 28 → 60 → 89) while the
next step falls (0.47 → 0.09) until the loop settles at zero.

Holding each of 1–3 at baseline (and 4) makes RAS/RAF 80x read UP and KO read
DOWN, which is up to ~19 experimental and ~24 curator RAF cases. No single
hold does it alone.

**Two earlier results that must not be ignored:**
- Holding capped-bundle structure (`LNG_CAP_POOLS`, specs/036/037) cost
  held-out −35 to −162 catalog-wide.
- De-duplicating a catalyst that is also a substrate (specs/012) cost −15.

A rule that fixes RAF by any means that reproduces those will lose elsewhere.

## Questions

1. **State one rule** (or show that none exists) that covers 1–3. Candidate
   wording: "a node's own downstream products must not multiply back into it".
   Make it precise:
   - what counts as "own downstream";
   - what "multiply back" is (a multiplicative factor in an AND, a member of a
     product pool, the numerator-tracking part of a divide inhibitor);
   - what the node reads instead (baseline? its value at entry to the
     component? the part not explained by the feedback?);
   - how it relates to specs/039's principles, specs/018 (derived-edge closures
     read at entry), specs/022 (self-contained inhibitors) and the loop
     knife-edge (specs/013/014).
2. **When is downstream feedback genuine biology** that must act? MEK
   phosphorylating RAF is real positive feedback. Is the rule removing real
   feedback, or only its runaway multiplication? What is the magnitude of the
   feedback after the rule, and is it bounded?
3. **Is it detectable** by the solver from the network alone (reachability in
   the component), or does it need the generator? Keep it label-independent.
4. **Work the RAF case numerically** under the rule: KRAS 80x, BRAF 80x, NF1
   KO, KRAS KO, at the activated RAF dimer and the p-MAPK dimers. Show which
   amplifiers it neutralises.
5. **Say what it predicts OUTSIDE RAF**, and how to check before measuring that
   it does not reproduce the cap-pools or de-duplication losses:
   - which other components have recycled pool members, bundle leaves fed only
     by their downstream, or leaf-sharing inhibitors;
   - a census method.
6. **What can go wrong:** removing load-bearing feedback, and label dependence.
   Which pre-registrable gates, and which census, should come before an arm?

Constraints: no parameter fitted to the evaluation; baseline exact; the default
otherwise byte-identical when the rule is off. Give numbers.
