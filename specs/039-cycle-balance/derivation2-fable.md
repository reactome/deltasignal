# specs/039 problem 2 — derivation 2 (Fable, blind)

2026-09-28, from `problem2.md`, `research.md`, `cycles.jl`, the generator at
1358d03, the `pools039b` RAF bundle and Reactome R97. The other derivation was
not read.

## The rule in one paragraph

Track **one protein** (a reference entity, R) through every reaction that carries
it. The cycle a curator drew in several steps is a strongly connected component
of that R-graph. Its **states** are the least-bound form of each modification
signature (S, S\*; RAS:GDP, RAS:GTP); the enzyme complexes on the way (S:E,
RAS:GTP:GAP) are **intermediates**. A path state → intermediates → state is one
transition of the pre-registered chain, with drive = the product of its steps'
drives. Intermediates carry **no baseline share** and read the flux through
them ([ES] = v/k_cat). The enzyme is not conserved. This reduces exactly to the
pre-registered two-form rule, so φ₀ = 0.1 keeps its meaning and no parameter
enters.

## 1. Detection (Q1)

For each reaction in the pathway's event closure and each reference entity R
reached from an input or output by `hasComponent|hasMember|hasCandidate*` →
EWAS → `referenceEntity`:

- **R-step:** exactly one input PE and one output PE contain R, both at
  stoichiometry 1, and they differ. Binding, catalysis inside a complex,
  release, translocation and single-step modification qualify; synthesis,
  degradation (no R output), dimerisation and two-R-input reactions do not.
- **Signature** σ_R(PE) = modified residues on R's EWAS leaves ∪ SimpleEntity
  components of the complex whose *direct* children include R's leaf (GTP/GDP
  attach to RAS in `p21 RAS:GTP` and not to the GAP in `p21 RAS:GTP:RAS GAPs`,
  which is nested; checked in R97).
- **Pool candidate:** an SCC of the R-graph with ≥ 2 forms and ≥ 2 distinct
  signatures. An SCC with one signature is a **carrier** loop (the enzyme's
  E → E:S → E, a scaffold, a cyclin): not a pool, but see §5.
- **States:** in each signature class, the forms with the fewest non-small-
  molecule slots (flattened). Everything else is an intermediate. Base state:
  the pre-registered rule (residues, donor, components) over states.
- **Transitions:** simple directed paths state → intermediates\* → *other*
  state inside the SCC, ≤ 6 steps (longer: dropped, counted). Bind/unbind
  paths returning to their own state are dropped. A path is **enzyme-driven**
  if any step has a catalyst or joining input containing a protein other than R
  (self-catalysis, e.g. RAS's intrinsic hydrolysis with `catalyst = RAS:GTP`,
  is not).
- **Node level** as now: every step must exist as a reaction node between the
  form nodes; each uuid path is a parallel transition. Shipped: forms (role
  state | intermediate, is_base), paths (from, to, step order, reaction node,
  enzyme-driven), carriers (enzyme free-form node, release nodes to ignore).

**Which protein's pool** (a complex holds several): the one whose signature
changes across the SCC. In RAS:GTP:GAP only RAS's changes.

## 2. A step for two proteins (Q2)

E + S → S:E is an R-step for both. For S it is a pool step (drive = E's fold);
for E a carrier step, and carriers write nothing but their free form (§5), so
the reaction node and S:E are written once, by S's pool. If both proteins are
modified in the same SCC (MAPK12 phosphorylates PTPN3 inside PTPN3:p-MAPK12,
R-HSA-8868118) a node is claimed by two pools: remove it from both, recompute,
count. R97 stId level: 13 such forms in 92 pathways (7 in WNT, 4 in insulin).

## 3. Baseline over multi-step forms (Q3)

π₀ over **states** is the pre-registered rule: ρ^m normalised, ρ = φ₀/(1−φ₀),
m = signature distance (0.9 / 0.1 for two signatures; forms of one signature in
two compartments split it equally). **Intermediates hold zero baseline share**
(β₀ → 0). Justification, not fitting: an enzyme–substrate complex is a transient
whose occupancy is set by the flux through it; any baseline share β₀ makes
k_cat = J₀/β₀ finite and caps the cycle's flux. Six-form ring
S → S:E → S\*:E → S\* → S\*:P → S:P → S, kinase 80×, S\* reads:

| β₀ (bound share per bound form, relative to its free form) | 1/9 | 1e-2 | 1e-3 | → 0 (lumped) |
|---|---|---|---|---|
| S\* | 3.67 | 7.77 | 8.85 | **8.99** = two-form rule |
| phosphatase KO bound on S\* | 12.2 | 10.2 | 10.02 | **10 = 1/φ₀** |

Only the limit reproduces the pre-registered numbers and bound. Baseline rates:
k_ij = J₀_ij/π₀_i, with J₀ the stationary flow of the random walk on the state
graph whose exit weights are 1 (enzyme-driven) or ε_int = 1e-3 (a non-enzyme
exit beside an enzyme-driven one *from the same state*: per exit, not per pair,
so RAS's intrinsic hydrolysis is weighed against the GAP *path*). For a pair
this is the current √ rule up to a constant; for any graph π₀Q₀ = 0 holds by
construction.

## 4. Drives (Q4)

Transition drive u = Π over steps of (the step evaluated with its source form
held at baseline) / (its baseline): hill_sat product, 0 absorbing. The binding
step's drive is E's fold with S at baseline: right at first order (mass action
in E, E in excess). Self-catalysed and release steps read 1 unless inhibited.
Rate = k(u + 1e-9); πQ = 0; state = baseline × s × modifier × π/π₀; a path's
flux fold = u·π_from/π₀_from, read by its reaction nodes and intermediates (an
intermediate on several paths reads their sum over the baseline sum). Supply,
pins and the unapplied catalyst ⊣ source-form edge are unchanged.

## 5. The enzyme (Q5)

Carrier free form E: producers that are pool release steps are removed; it
reads baseline × s_E × modifier. Today that release edge closes a gain-1 loop
of its own (NF1 → bind → RAS:GTP:GAP → release → NF1). Its complexes read S's
pool. Conservation of E is **not** needed at first order: every fixture
direction holds without it. Cost: substrate 80× reads its complexes and the
catalytic flux at 80× (no fixed E_total, no V_max), and sequestration of E by
a competing substrate is invisible.

## 6. Fixtures (Q6), φ₀ = 0.1

**RAF, hand-traced on `pools039b`.** RAS pool = {RAS:GDP (G, base), RAS:GTP (T),
RAS:GTP:RAS GAPs (TX, intermediate)}. G → T: `5672965` GEF ×4 copies,
`9649735` intrinsic ×4 (ε_int). T → G: `9649736` intrinsic, self-catalysed ×4
(ε_int); paths `5658435` bind → TX → `5658231` release ×48 (12 GAP variants × 4
RAS; 12 carry NF1 via SPRED1/2/3:NF1). Exits: RAF, PI3K, BRAP binding. Supply:
`9649733`. Carriers: the 12 GAP nodes and the `RAS GEFs` set pool.

| perturbation | RAS:GDP | RAS:GTP | RAS:GTP:GAP |
|---|---|---|---|
| none | 1 | 1 | 1 |
| NF1 KO (12/48 paths at 0) | 0.968 | **1.290 UP** | 0.968 |
| all GAPs KO | 0.001 | 9.99 | 0 |
| GAPs 80× | 1.110 | 0.014 | 1.110 |
| SOS1 80× (GEF node g = 80; 100 if the set pool caps) | 0.112 (0.092) | **8.99 (9.17) UP** | 8.99 |
| GEF KO | 1.111 | 0.001 | 0.001 |
| KRAS 80× (s = 80 via the maturation chain) | **80** | **80** | 80 |

NF1 KO is UP at the 1.15 cutoff only because SPRED variants give NF1 3 of 12
slots; with equal weight over the 10 curated GAPs it would read 1.103, NORM.

**Kinase/phosphatase through E:S complexes** (six-form ring, lumped): kinase KO
→ S\* = 0, S = 1.111, complexes 0; phosphatase KO → S\* = 10 (= 1/φ₀), S = 0;
kinase 80× → S\* 8.99, S 0.112, complexes 8.99 (flux); both enzymes 80× →
S, S\* = 1, complexes 80. **Baseline:** π₀Q₀ = 0, every drive and flux 1:
exactly baseline.

## 7. Coverage (R97, stId level, 92 pathways)

**78 pools in 37 pathways**: 32 multi-step (with intermediates) in 24
pathways, the class the current rule cannot see; 46 single-step (its 55 pairs
after stoichiometry and dedup); 5 with > 2 signatures. Expect node level to
halve this, as 55 → 25 did. Experimental axis (`bench/analysis/_common.py:
EXP_PATHWAYS`; cases from the MP-BioPath supplementary workbook,
`--supplementary-workbook`): 7 of 9 pathways have pools (12): TP53 3 (2
multi), WNT 3 (2), Prophase 2, HRR 1 (multi), G1/S 1, S phase 1, PIP3/AKT 1;
ERBB2 and Checkpoints none. The experimental ≥ +15 gate stays a stretch; RAF
is curator-axis only.

## 8. What can go wrong (Q7)

- **Set-member enzyme KO reads (N−1)/N.** UP at 1.15 needs N ≤ 6 (N = 6:
  1.176; 7: 1.148; 10: 1.099). NF1 passes by variant multiplicity, not biology.
- **Set-closed cycles are invisible at stId level.** MAPK in RAF: DUSPs act on
  the set `p-T,Y MAPK monomers and dimers` → set `MAPK monomers and dimers`; no
  reaction returns to `MAPKs`. Node-level set expansion may close it; count.
- **Cross-pathway halves** (PTPN3 dephosphorylates MAPK12 in RAF; the kinase is
  elsewhere): no SCC, no pool.
- **Stable complexes read as transients.** p-SMAD2,3:SMAD4 is an intermediate
  of SMAD4's ubiquitin pool and reads flux: direction right, magnitude not.
  HRR's 13-form SCC (11 intermediates) is a repair machine, not an enzyme cycle.
- **Signature by nesting.** A flat complex [RAS, GTP, GAP] would give the GAP a
  nucleotide signature; R97 nests these, but the generator must check.
- **Shared forms** (13) drop both pools; paths > 6, pinned intermediates and
  supply into an intermediate fall back, counted.
- **Product of step drives** over-multiplies when two steps of one path are
  both perturbed (rate-limiting physics is nearer min); double perturbations only.
- φ₀ ≈ 1 sites, autocatalysis, influence scores: as pre-registered.
