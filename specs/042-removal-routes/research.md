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
