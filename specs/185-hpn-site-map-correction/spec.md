# specs/185: two residue-numbering errors in the HPN-DREAM site map, corrected and re-measured

**Status:** registered 2026-11-27, after the audit below and before any re-run.
**Requested by:** Adam, 2026-11-27 ("ok", on correcting the gaps specs/184 traced and re-measuring).
**Post hoc.** HPN-DREAM is FIT and specs/184's figures have been seen. This corrects inputs; it is not a test.

## The audit

Every one of the 47 readouts in `specs/112-hpn-dream/site_map.tsv` was compared with every phosphorylated EWAS
Reactome 98 curates for its genes (coordinate and reference isoform). 36 readouts match a curated form at the mapped
coordinate. A coordinate was changed only where the antibody's stated residue is provably the same residue under
another numbering, checked against the UniProt sequences:

| Readout | Was | Now | Evidence |
|---|---|---|---|
| p70S6K pT389 | RPS6KB1 T389 | T389 and T412 | The 15-residue window around T389 of isoform 2 (P23443-2) sits at T412 of the canonical sequence (P23443). Reactome curates the mTORC1 target on `P23443-1` at T412 and a second form on `P23443-2` at T389. bench/analysis/127 already maps this antibody to T412. |
| HER3 pY1298 | ERBB3 Y1298 | Y1289 | Residue 1298 of P21860 is a proline; the tyrosine in that stretch is Y1289, which Reactome curates. |

Checked and left unchanged, because no curated form is the same residue: BAD pS112 (mouse numbering; human S75, not
curated; Reactome has S91 and S99), FOXO3 pS318/S321 (Reactome has S315 and S325, different residues), 4EBP1 pS65,
p27 pT198, mTOR pS2448, NDRG1, PEA15, Rictor, YB-1, ACACA S79 (no curated phospho form), JNK (curated without a
coordinate).

## Re-measurement

- specs/184's arm is re-run unchanged at this spec's commit (`run_arm.sh --axis hpn_stim_predict`, build
  `20261126-2322_d6d026b_top98`) and re-scored with `bench/analysis/184/stim_184.py`, unchanged.
- **Control:** every prediction row for the 45 other readouts must equal specs/184's file (`results/2e399999/stim184`):
  same rows, same calls, `pred_ui` within 1e-6 relative. A difference stops the analysis until it is traced.
- Reported beside specs/184: all of its registered figures, and the two readouts' cases one by one.
- G1 (model − reachability, 24 conditions) is read again and labelled post hoc whatever it shows.

Written before the run: p70S6K gains a second node and is called UP where a ligand reaches mTORC1; HER3 becomes a
scored readout, called UP under NRG1. Recall rises; I do not expect the gate to change side.

## Not corrected here

specs/184's other two gaps are network questions, traced in research.md and filed as issues, with no change made:
IGF1's handoff from IRS on the receptor to PI3K, and the ERK → RSK branch no ligand reaches.

## specs/112

The same map scored sub-challenge 1A, where p70S6K is the canonical mTOR target. The recorded 0.693 stays the
scored-once result. A corrected re-score needs the Release 97 `hpn112` build served again; whether it is run is
recorded in research.md.
