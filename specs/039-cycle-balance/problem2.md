# specs/039 problem 2: modification cycles curated in several steps

Given to two derivations independently (blind), as for `problem.md`.

## What exists

specs/039 solves a **pool**: the forms of one protein joined by curated
transition reactions. It is solved as a Markov chain, π Q = 0, with:
- rates k_ij · u_ij, where u is the transition reaction's drive (the reaction
  with its source form held at baseline, over baseline);
- baseline π₀ by detailed balance, with each step away from the base form
  weighted φ₀/(1 − φ₀), and φ₀ = 0.1;
- form value = baseline × s × π/π₀, where s is the supply fold.

A transition must convert one form directly into another: A + small molecules
→ B. The rebuilt catalog has 25 such pools (12 of them RHO GTPases).

## The problem

Reactome usually curates modification in steps:

    S + E → S:E           (binding)
    S:E → S*:E            (modification; ATP → ADP etc.)
    S*:E → S* + E         (release)
    S* + P → S*:P → S:P → S + P    (the reverse, via a phosphatase/GAP P)

or with the enzyme acting on a pre-formed complex: RAS:GTP + NF1 →
RAS:GTP:NF1 → RAS:GDP + NF1 (+ Pi). The single-step detection cannot see
these, so:
- NF1 knockout in RAF/MAP reads "no change" in RAS:GTP, where both biology and
  the curators say UP;
- 164 of 203 curated A ⇄ B entity pairs are enzyme binding steps E ⇄ E:S,
  excluded because they are not interconversions of one protein.

Complication: the complex S:E contains both S's protein and E's protein.
E + S → S:E is a step for S (S → S:E) and also for E (E → S:E). Enzymes are
conserved too.

## Questions

1. Which reaction sequences form a pool, stated so a generator can detect them
   from Reactome (inputs, outputs, catalysts, reference entities, modified
   residues)? Which protein's pool is it, when a complex holds several?
2. How is a reaction that is a step for two proteins treated? (It is one node
   with one value.)
3. What is the baseline π₀ over multi-step forms (free S, S:E, S*:E, S*, S*:P,
   S:P)? φ₀ is defined for "modified vs unmodified". How are the transient
   enzyme complexes weighted, and with what justification? No fitting against
   the evaluation.
4. What drives each step? For example, the binding step's drive = E's fold
   under the existing AND semantics with S held at baseline. Is that right?
5. What do the enzyme's own free form and its complexes read, and is its
   conservation needed at the first order?
6. Check the rule on the fixtures:
   - RAF: NF1 KO → RAS:GTP UP; SOS1 80x → RAS:GTP UP; KRAS gene 80x → both
     forms UP.
   - A kinase/phosphatase cycle through E:S complexes: kinase KO → S* ≈ 0;
     phosphatase KO → S* up to its bound (1/φ₀).
   - Baseline exactly baseline.
7. What can go wrong: detection false positives, and coverage.

Keep it to what a pre-registration can state. Give numbers for the fixtures.
