# specs/185 research: the site-map correction re-measured, and what the traces found instead

Run 2026-11-27. Build `20261126-2322_d6d026b_top98`, predictions `results/27b16759/stim185` (54 solves, 51 converged,
363 rows). Output outside the repo under `analysis/184/stim_185.*`. **Post hoc on FIT data.**

## Verdict

- **Control passes:** all 349 rows of the 45 unchanged readouts equal specs/184's file (0 differences).
- **The correction is right and changes little.** G1 is +0.012 [−0.012, +0.036] (specs/184: +0.015 [−0.005, +0.036]);
  enrichment 1.97 [1.18, 3.53]; direction 78 of 84.
- **HER3 pY1289 is now read:** 120 nodes, 67× under NRG1 (measured UP in 2 of 3 lines), unchanged under the other
  ligands.
- **p70S6K is still silent.** With the T412 form included the readout still reads exactly 1 under every ligand. My
  prediction that it would be called UP was wrong. The reason is below, and it is not the site map.

| | specs/184 | specs/185 |
|---|---|---|
| pairs | 700 | 718 |
| measured changes | 131 | 138 |
| model: right of calls | 76 of 308 | 78 of 311 |
| reachability: right of calls | 94 of 438 | 99 of 450 |
| model balanced accuracy | 0.412 | 0.408 |
| reachability balanced accuracy | 0.397 | 0.395 |

By engagement, δ doubled and the late window move as little (δ doubled: +0.046 [−0.001, +0.085]).

## Why p70S6K does not move: made copies and used copies of one entity do not meet

Reactome has one entity "Active mTORC1 complex" (R-HSA-165678). The Signal Transduction network has 16 copies of it:

| Copies | Made by | Used by |
|---|---|---|
| 4 | Formation of active mTORC1 complex (RHEB:GTP route) | catalyse S6K, 4EBP1 and the other mTORC1 substrates |
| 4 | Activated AKT phosphorylates AKT1S1 (PRAS40) | nothing |
| 8 | nothing (starting copies) | inputs of the AMPK and one other reaction |

The four copies AKT makes have exactly the member content of copies that are used. So AKT's activation of mTORC1
ends in a node nothing reads. The only route from EGF to p70S6K T412 is 31 edges long, through PLCγ, calcium,
calmodulin, a phosphodiesterase, AMP and AMPK, and it is negative; insulin has no route at all. (Reactome also
curates no regulation from AKT-phosphorylated TSC2 to RHEB: the product only binds 14-3-3. That part is curation.)

The same pattern cuts two other textbook routes in this network:

- **RAS → RAF.** 93 of 148 RAS:GTP copies are made by exchange reactions and used by nothing (the KRAS and NRAS
  copies); each has an identical-content copy that RAF, PI3K and the GAPs do use.
- **RAF cascade → free phospho-ERK.** "Dissociation of RAS:RAF complex" releases the set "p-T,Y MAPKs" as one node
  with no out-edge; the next reaction's member copies are fed from elsewhere. Growth factors reach 5 of 53
  phospho-ERK nodes, all inside scaffold complexes. Free phospho-ERK, its dimers, nuclear ERK and p90RSK are
  unreachable from any ligand (issue #300).

**How common (all 28 networks of the build, from the catalog tables):**

- 5,625 entities have variant copies. For 270, a copy that a reaction makes is used by nothing while another copy of
  the same entity is used; for 260 of those the unused copy has exactly the content of a used copy; for 181 a used
  copy is made by nothing.
- 43 sets released whole by a reaction feed nothing while their members are consumed elsewhere.
- Largest counts: Immune System 73, Signal Transduction 60, Metabolism of proteins 24, Gene expression 21.

This is not counted by specs/046's "cut hand-offs" check, which is why that reads 0.

## IGF1

IGF1's own break is different (issue #301): the IRS complex on the IGF1 receptor is linked to "p-Y-IRS bind PI3K"
only by `precedingEvent`; no reaction releases IRS.

## What follows

The unconnected copies are a construction defect at the centre of signalling, found on data where the answer is
textbook. They may also account for part of the unrouted changes on Lun (1,001 of 1,437) and korkut; that is not
measured here. The next spec registers a fix and measures it on every axis under the usual gates. Earlier merges of
copies re-welded cycles (specs/053: −493), so the fix has to be one-way and narrow, and judged case by case.

specs/112's corrected re-score was not run: it needs the Release 97 `hpn112` build served. p70S6K is a direct mTOR
target, so the corrected sub-challenge score may differ from the recorded 0.693; the re-score is still owed.
