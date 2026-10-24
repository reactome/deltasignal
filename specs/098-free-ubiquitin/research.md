# specs/098 research

## Arm `ubi`: void

`results/2d0c227/ubi` ran on a branch cut from a local `main` that had not been pulled after
PR #163, so it **lacked specs/095's default** (`DS_SELF_INHIBITOR_COPIES=on`), while its baseline
`results/b61b365/dflt095` has it. Its differences mix the two rules: the 9 PTEN overexpression
"breaks" it shows are specs/095's absence, not this rule (`pten_ub.jl` reproduces 83.5× with free
ubiquitin propagating). Not scored. The branch was merged with main (`specs/095` default on) and
the arm re-run as `ubi2`.

Lesson recorded: before an arm, check that the arm's commit contains the baseline's defaults
(`git merge-base --is-ancestor <baseline sha> HEAD`).
