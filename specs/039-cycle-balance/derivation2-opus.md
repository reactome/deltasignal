# problem2: derived independently (Opus), before reading the Fable derivation

## 1. Detection: the pool is the protein whose own state cycles

Work per protein P (a reference entity of an EWAS leaf).

- **Form graph of P.** The nodes are the non-small entities that contain P. There is an
  edge X → Y through reaction r when X is r's **only** input containing P and Y is
  r's **only** output containing P.
  - Other participants are **partners**: an enzyme or cofactor on the input side,
    a released partner on the output side.
  - A reaction with two P-containing inputs (dimerisation) is not an edge, and nor
    is one with no P-containing output (degradation).
- **Candidate pool.** A strongly connected set of P's forms with at least 2 forms.
- **Whose pool it is.**
  - The *core* is the set of proteins present in **every** form of the candidate.
  - The candidate is a pool only if **at least two forms contain exactly the
    core**. Those are the "free" forms, such as S and p-S, or RAS:GDP and RAS:GTP
    (small molecules allowed).
  - Forms that also carry other proteins are **intermediates**, such as S:E or
    RAS:GTP:NF1.
  - For an enzyme E, the cycle E → E:S → … → E has only one free form (E), so it is
    **not** a pool. So the binding steps belong to the substrate's pool, not the
    enzyme's.
  - This one rule removes both the "which protein?" ambiguity and most double
    membership.
- **Orientation.** The base form is the free form chosen as in amendment 1:
  residues first, then the donor, then components.
- **Single-step pools** (amendment 1) are the special case with no intermediates.

## 2. A reaction that is a step for two pools

After §1 this is rare: both proteins would have to have their own state cycle
through it (reciprocal modification). Such reactions are dropped with their
pools, and counted, as amendment 1 does for multi-use nodes. The drop is
reported per pathway.

## 3. Baseline π₀ and baseline rates

Let r = φ₀/(1 − φ₀) = 1/9.

- **Free forms:** as before, r^(modification distance from the base free form)
  over the free-form graph. Two free forms have weights 1 and r.
- **Intermediates:** r × the weight of the nearest free form upstream along the
  directed arc. Enzyme–substrate complexes are transient and substoichiometric
  at rest, the same "minority at rest" judgement as φ₀, so no new parameter
  is introduced.
- **Normalisation.** For the ring S → S:E → S*:E → S* → S*:P → S:P → S the
  weights are
  1, r, r², r, r², r (Z ≈ 1.358), so π₀(S) = 0.736 and π₀(S*) = 0.0818.

**The baseline rates must make π₀ stationary on a graph that is not
reversible.**
- The current `k_ij = sqrt(π0_j/π0_i)` does this only when every step has a
  reverse. On a directed ring it gives π_i ∝ 1/k_out,i ≠ π₀. **This is a bug for
  multi-step pools, and the code must change.**
- Rule: choose a positive baseline **circulation** f on the transitions, so that
  inflow equals outflow at every form. Then k_ij = f_ij / π0_i.
- Construction (deterministic, label-independent through sorted stable ids): for
  each transition i → j, add one unit of flow along it and along the shortest
  directed path back from j to i.
- Parallel reactions for one pair share the edge's flow by weight, catalysed 1
  and intrinsic 1e-3 (amendment 2).
- A two-form reversible pool gives f_AB = f_BA and reproduces the current rates
  up to scale, so the amendment 1 fixtures are unchanged.

## 4. Drives

These are unchanged: each step's reaction is evaluated with its source form at
baseline.
- **Binding** (S + E → S:E): drive = E's fold. The cofactor is inert.
- **Chemistry inside the complex:** usually uncatalysed, so the drive is 1, or it
  follows a regulator if one is curated.
- **Release:** drive 1.

So the enzyme acts through the binding rate, which is mass action in E. That is
first order in E: enzyme in excess, low occupancy.

## 5. The enzyme itself

- Free E, and E's other uses, are computed by the existing semantics. E is
  produced by the release reactions, whose nodes carry flux, and by any other
  producer.
- E's conservation is **not** modelled. Under substrate overexpression the
  release flux rises, so free E reads UP. This is the pre-existing behaviour, not
  a new one (specs/035 held such inputs at baseline and lost).
- The loop E → bind → … → release → E stays iterated. It is no longer a gain-1
  knife-edge, because the pool bounds the intermediate.

## 6. Fixtures (φ₀ = 0.1)

**Ring S → S:E → S*:E → S* → S*:P → S:P → S**, unidirectional, one binding step
per enzyme. The stationary π_i ∝ f_i / k_out,i with k ∝ drive.
- **Baseline:** π = π₀, every fold 1.
- **Kinase KO** (the S + E binding drive is 0): mass leaves the S*-side arcs and
  cannot return, so it all collects at S. S fold = 1/0.736 = 1.36; S* → 0
  (the regulariser only).
- **Phosphatase KO:** mass collects at S*, fold 1/0.0818 = **12.2** (the bound,
  about 1/φ₀).
- **Kinase 80x:** the S → S:E rate is ×80, so S is drained. With the circulation
  f = 1 on the ring, the stationary share is proportional to 1/k_out:
  - S drops by 80×, to π(S) ≈ 0.736/80, normalised.
  - The other forms keep their relative weights, rising by the renormalisation:
    Z' = 0.736/80 + 0.264 = 0.273, so each non-S form is ×3.66. S* reads 3.66
    (UP), and S reads 0.034.
  - This is weaker than the single-step 8.99. The bound on any enzyme-up
    perturbation is 1/(1 − π₀(S)) = 3.79, because S's mass is shared among the
    five other forms. Corrected in writing, before the Fable derivation was read.
  - Ring-wide, the modified side takes 1 − (intermediates on the base side). A
    **stated consequence**, since the intermediates dilute the gain.

**RAF**, assuming NF1's arc is RAS:GTP + NF1 → RAS:GTP:NF1 → RAS:GDP + NF1 and
SOS1 exchange is catalysed on the free form:
- **NF1 KO:** RAS:GTP's only exit (if NF1 is the only GAP) is blocked. RAS:GTP
  fold ≈ 1/π₀ ≈ 10: **UP**.
- **SOS1 80x:** RAS:GTP is UP, bounded as above.
- **KRAS gene 80x:** s = 80, and every form is 80x.

## 7. Failure modes

- **Signalling complexes mistaken for intermediates.** RAS:GTP:RAF is an
  intermediate only if it returns to the cycle. Otherwise it lies outside the
  strongly connected set and is computed as now.
- **The core is ambiguous** where every form holds several proteins (SMAD2/3:SMAD4
  ⇄ p-SMAD2/3:SMAD4). The core is then the complex, which is fine, but the
  pre-registration should say so.
- **Compression by intermediates** (the kinase 80x fixture). This can cost cases
  where curators expect strong UP. Set it against the single-step rule's 1/φ₀.
- **Multiple GAPs:** one GAP KO leaves the others' arcs, so the fold is partial,
  which is correct.
- **Coverage** is unknown until counted. The generator needs the P-containment of
  every form (EWAS reference entities, which it already queries).
