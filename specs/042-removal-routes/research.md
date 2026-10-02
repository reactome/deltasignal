# specs/042 research

## Neo4j check before any rule (Adam's rule, 2026-10-02)

Before writing a rule, confirm in Reactome that the topology is curated as
observed, and is not a logic-network-generator artefact.

### WNT / APC, the main motivating M1 case: a GENERATOR ARTEFACT

**In Reactome** (Signaling by WNT, R-HSA-195721), free CTNNB1 [cytosol]
(R-HSA-448839) is **one entity in a closed binding/release cycle:**
- an input of R-HSA-195304, "Association of beta-catenin with the destruction
  complex";
- an input of R-HSA-201669, "Beta-catenin translocates to the nucleus";
- an output **only** of R-HSA-201685, "Beta-catenin is released from the
  destruction complex";
- no synthesis of β-catenin is curated in the pathway.

**In our build** (`20260928-1110_06ccb63`) it is split into three uuid copies:

| copy | produced by | feeds |
|---|---|---|
| `9aa158df` | nothing (in-degree 0) | association only |
| `094a7dc9` | release only | nuclear translocation |
| `2c39329c` | a dissociation | nothing (dead end) |

- **The curated cycle is severed.** APC KO → 0 follows from the split (no
  complex, so no release), not from a missing removal route.
- **Rejoined,** it is a binding/release cycle in which conservation gives the
  biological answer: no destruction complex leaves all β-catenin free, which
  is UP. That is specs/039's pool logic, whose detection currently excludes
  one-signature (binding) loops as "carriers".
- **Consequence:** derivation-opus.md's rule R1 (a depletion edge per
  non-returning removal) was built on this case as stated in problem.md. Its
  premise does not hold for WNT. The remaining M1 cases (p27, RB1, HDR branch
  competition, CDKN2A) are being checked the same way before any rule is
  chosen.

### Every M1 case checked against Reactome (Fable, `derivation-fable.md` §0)

| cases | verdict | mechanism |
|---|---|---|
| WNT APC/AMER1 (8) | **generator artefact** | CTNNB1 [cytosol] and AMER1 [cytosol] are each one entity in Reactome (bind 195304, release 201685, nuclear import 201669, proteasome exit 2130282). Our build makes 3 copies each. |
| TP53 CDKN2A (10) | **generator artefact** | TP53 Tetramer 3209194 is one entity. Our build has 4 copies; the p14ARF and USP7 products dead-end in dissociations. |
| HDR BLM/RTEL1/PALB2/RAD51 (10), S CDKN1B (1), G1 RB1 (1), G1 CDKN1B (1) | **curated as-is** | Shared states with 2–3 consuming fates in Reactome. RB1:E2F and p27:CDK2 are dead ends in Reactome itself. |

**What the verified cases call for:**
- **A generator fix first: rejoin.** One node per stId per pathway where that
  stId is both produced and consumed.
- **Rejoining alone does not fix WNT.** The rejoined loop has no supply, so it
  goes to 0. With a closed specs/039 pool it reads 1.11 (NORMAL), because
  APC's lever is degradation, which a closed pool cannot see.
- **The curated cases need a solver rule.** Fable proposes three parts, each
  separately switchable and baseline-exact:
  - **R1, open-pool turnover:** total = supply / mean removal drive over the
    curated removal routes, with the specs/011 clamp.
  - **R2, fate competition** at a state with ≥ 2 transforming fates. The
    perturbed fate is never divided, which is the precise difference from naive
    depletion and general consumption.
  - **R3, dead-end sequestration as an exit.** A separate arm.
- **Census (canonical build):**
  - root-split stIds: 499 in 50 pathways;
  - dead-end produced copies: 614 stIds;
  - protein removal reactions: 97, on 243 species in 40 pathways;
  - branch states with ≥ 2 transforming fates: 1,302 in 86 pathways;
  - dead-end bindings: 1,257 in 89 pathways.
- **Worked cases:**
  - WNT: APC KO → β-catenin 10 (UP); APC 80x → 0.1 (DOWN).
  - HDR: PALB2 80x → SSA 0.025; KO → 2.0.
  - CDKN1B 80x → p-RB1 0.037.
  - RB1 KO → E2F1 pool 10.

**Not yet reconciled with derivation-opus.md** (R1, committed earlier). The
order is fixed: **the generator rejoin is designed first**, as its own spec
change, because it is an artefact. A solver rule for the curated cases is
pre-registered only after the rejoin is measured.

### Correction (2026-10-02, diagram check, see specs/043)

The WNT CTNNB1 split is **faithful to Reactome's diagrams**, not a generator
artefact.
- Association with the destruction complex and release / nuclear import are
  drawn in separate sub-pathway diagrams, on separate glyphs.
- So the WNT M1 cases return to the solver/semantics column. Free cytosolic
  β-catenin is curated only as a release product, and APC's real lever
  (degradation) is not represented as acting on it.
- The TP53 tetramer split **is** an artefact: a diagram draws one shared glyph.
