# specs/115: why the frozen model does poorly on the phospho datasets (analysis)

**Status:** analysis only, 2026-10-27. No rule, no setting and no score is changed. korkut2015 and decryptM are FIT.

Adam, 2026-10-27: "The big question is why we did so poorly ... we need to be clear on what deficiency we have with our
method (DeltaSignal + Reactome networks) ... also how are we doing compared to others on these other datasets."

**Questions.**
1. Where do the errors fall? Every scored case on the top-level build (specs/114) is put in one class: changed with no
   route, changed and right, changed but called NOT, changed in the wrong direction, unchanged and right, unchanged but
   called changed (`bench/analysis/error_anatomy_115.py`).
2. Do predictions carry information near the drug's target and lose it further out? korkut cases are split by the
   shortest hop distance from any pinned target node to a readout node in the same network
   (`bench/analysis/korkut_hops.py`, graph only, the bench's own pin and readout resolution).
3. What have other methods scored on korkut2015, decryptM and LINCS P100, and were they fitted to that data?
   (Literature, read-only.)
