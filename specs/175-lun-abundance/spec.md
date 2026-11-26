# specs/175: the Lun 2019 total readouts under the adopted abundance rule

**Status:** registered 2026-11-26, before any new prediction was made.
**Requested by:** Adam, 2026-11-26 ("go", on option 1 of specs/174's candidates).

## Why

specs/174 traced total β-catenin (0 of 44 right change calls on Lun half A) to the readout, not the model. The Lun
predictions (specs/127, 2026-11-03) read a total antibody as every node containing the protein. specs/134 adopted the
abundance rule for total antibodies on the korkut and CPPA benches on 2026-11-04 (`DS_TOTAL_READOUT=abundance`, now
amendment 3). The Lun bench was never brought into line. This applies an already-adopted rule; nothing is fitted.

## Change

`bench/lun_bench.py` reads `DS_TOTAL_READOUT` exactly as `korkut_bench.py` does:
- total readouts resolve to protein EWAS only;
- each is cut to its abundance copies (`total_readout_nodes`);
- the default is `abundance` (amendment 3);
- `forms` reproduces specs/127.

Reachability uses the same readout nodes. Phospho readouts are untouched. `run_arm.sh` passes the knob for
`lun_predict`.

## Run

One arm on the top114 build (`20261027-1346_8b6527e_top114`), the specs/127 protocol, code at this spec's
implementation commit. Neo4j Release 97 serves the gene lookup, as it did for specs/127.

**Control, inside the same arm.** Every phospho-readout row must equal the specs/127 file (sha256 `0c0ec033…`):
same rows, same `predicted` and `reach_predicted`, and `pred_ui` within 1e-6 relative. A difference stops the analysis
until it is traced. This checks that the solver at HEAD reproduces 0d16eb9 on these networks. The total rows are the
only intended change.

## Registered reading (half A only)

Total readouts measured on Lun: β-catenin only (cyclin B1 is in Table S2 but not in the binned file's targets).
Scored with specs/173's measured call on the P ∩ A pairs, against the specs/127 predictions. specs/161's rule applies,
because the purpose is removing false change calls.

- **Paired gate on the moved β-catenin cases.** Fixed > broken, sign test p < 0.05.
- **Sensitivity guard.** Over all P ∩ A pairs, recall of measured changes may fall by at most 0.01.
- **Reported:**
  - three-class balanced accuracy, all pairs and β-catenin alone (new rule, old rule, reachability under each);
  - R2 enrichment and lift (specs/173 definitions).

The implementation makes `abundance` the Lun bench default, as it already is for korkut and CPPA: specs/134 is the
adopted rule for total antibodies. If the gate fails, that is recorded, and Adam decides whether the Lun default
goes back to `forms`.

**Half B is not scored here.** It stays sealed for one scoring of the candidate Lun protocol, once the other specs/174
candidates (route length, hub proteins) have each been decided on half A under their own registrations. The Lun
predictions for half B will use whatever readout rule this spec settles.
