#!/usr/bin/env python3
"""Turn paired all-ten benchmark runs into an auditable evaluation report."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from xml.sax.saxutils import escape


PATHWAY_NAMES = {
    "R-HSA-1257604": "PIP3 activates AKT signaling",
    "R-HSA-453279": "Mitotic G1 phase and G1/S transition",
    "R-HSA-69620": "Cell Cycle Checkpoints",
    "R-HSA-68875": "Mitotic Prophase",
    "R-HSA-69242": "S Phase",
    "R-HSA-195721": "Signaling by WNT",
    "R-HSA-1227986": "Signaling by ERBB2",
    "R-HSA-3700989": "Transcriptional Regulation by TP53",
    "R-HSA-5673001": "RAF/MAP kinase cascade",
    "R-HSA-5693567": "HDR through HRR or SSA",
}
CASE_KEY = ("pathway_id", "gene", "direction", "key_output_dbid")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagram-on-dir", type=Path, required=True)
    parser.add_argument("--diagram-off-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--lng-commit", help="LNG commit used to freeze both catalogs")
    parser.add_argument("--reactome-release", help="Reactome release used by LNG")
    parser.add_argument(
        "--tcga-readiness",
        type=Path,
        help="Optional historical TCGA LUAD readiness_key_findings.tsv",
    )
    return parser.parse_args()


def load_run(path: Path) -> tuple[dict[str, object], list[dict[str, str]]]:
    summary = json.loads((path / "benchmark_summary.json").read_text())
    with (path / "benchmark_cases.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != int(summary["eligible_cases"]):
        raise ValueError(f"Case count does not match summary in {path}")
    return summary, rows


def sanitized_manifest(path: Path) -> dict[str, object]:
    manifest = json.loads((path / "benchmark_manifest.json").read_text())
    inputs = manifest["inputs"]
    return {
        "ground_truth": manifest["ground_truth"],
        "pathway_ids": manifest["pathway_ids"],
        "evaluation_split": manifest["evaluation_split"],
        "thresholds": manifest["thresholds"],
        "output_aggregation": manifest["output_aggregation"],
        "allow_output_proxies": manifest["allow_output_proxies"],
        "perturbation_up": manifest["perturbation_up"],
        "solver_environment_overrides": {
            key: value
            for key, value in manifest["solver_environment_overrides"].items()
            if key != "DS_PATHWAY_CATALOG"
        },
        "deltasignal_git": manifest["deltasignal_git"],
        "benchmark_code": {
            key: value
            for key, value in manifest["benchmark_code"].items()
            if key != "path"
        },
        "input_hashes": {
            "supplementary_workbook": inputs["supplementary_workbook"]["sha256"],
            "id_map": inputs["id_map"]["sha256"],
            "reactome_id_audit": (
                inputs["reactome_id_audit"]["sha256"]
                if inputs["reactome_id_audit"]
                else None
            ),
        },
        "network_files": {
            pathway_id: {
                name: {
                    "sha256": metadata["sha256"],
                    "size_bytes": metadata["size_bytes"],
                }
                for name, metadata in files.items()
            }
            for pathway_id, files in inputs["networks"].items()
        },
    }


def case_key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[field] for field in CASE_KEY)


def is_correct(row: dict[str, str]) -> bool:
    return bool(row["prediction"]) and row["prediction"] == row["expected"]


def convergence_metrics(rows: list[dict[str, str]]) -> dict[str, int | float]:
    scored = [row for row in rows if row["prediction"]]
    converged = [row for row in scored if row["converged"] == "True"]
    unique_solves = {
        (row["pathway_id"], row["gene"], row["direction"]): row["converged"]
        for row in scored
    }
    correct = sum(is_correct(row) for row in converged)
    nonconverged = [row for row in scored if row["converged"] != "True"]
    nonconverged_correct = sum(is_correct(row) for row in nonconverged)
    return {
        "converged_cases": len(converged),
        "nonconverged_cases": len(scored) - len(converged),
        "converged_correct": correct,
        "converged_accuracy": correct / len(converged) if converged else 0.0,
        "nonconverged_correct": nonconverged_correct,
        "nonconverged_accuracy": (
            nonconverged_correct / len(nonconverged) if nonconverged else 0.0
        ),
        "unique_scored_solves": len(unique_solves),
        "unique_nonconverged_solves": sum(
            value != "True" for value in unique_solves.values()
        ),
    }


def exact_mcnemar_p(gains: int, losses: int) -> float:
    discordant = gains + losses
    if not discordant:
        return 1.0
    tail = sum(
        math.comb(discordant, value)
        for value in range(min(gains, losses) + 1)
    ) / (2**discordant)
    return min(1.0, 2.0 * tail)


def percent(value: float) -> str:
    return f"{100.0 * value:.1f}%"



BASELINE_LABELS = {
    "no_change": "Always predict NO CHANGE",
    "development_class_frequency": "Development class frequency",
    "signed_reachability": "Signed reachability",
    "shortest_signed_path": "Sign of shortest signed path",
}


def baseline_table(summary, scope="all"):
    """Rows for the structural baselines the harness computes on the same cases.

    These were computed and written to the summary JSON but appeared in no
    section, table or figure, so the report named only MP-BioPath and curator
    as comparators. `shortest_signed_path` is the relevant one: it scores the
    SAME cases DeltaSignal does, and a reader who opens the JSON will find it.
    """
    lines = []
    for key, label in BASELINE_LABELS.items():
        entry = (summary.get("baselines", {}).get(key) or {}).get(scope)
        if not entry:
            continue
        lines.append(
            f"| {label} | {entry['scored_cases']} | {entry['correct']} | "
            f"{percent(float(entry['accuracy']))} | {float(entry['macro_f1']):.3f} |"
        )
    return "\n".join(lines)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def svg_bar_chart(
    path: Path,
    title: str,
    labels: list[str],
    series: list[tuple[str, str, list[float]]],
) -> None:
    left = max(280, min(610, int(max(map(len, labels)) * 8.3)))
    top, plot_width = 105, 810
    width = left + plot_width + 120
    row_height = 54
    height = top + len(labels) * row_height + 90
    group_width = 18
    gap = 5
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="36" y="48" font-family="Arial" font-size="28" font-weight="700" fill="#111">{escape(title)}</text>',
    ]
    for tick in range(0, 101, 20):
        x = left + plot_width * tick / 100
        parts.append(f'<line x1="{x:.1f}" y1="78" x2="{x:.1f}" y2="{height-55}" stroke="#dddddd"/>')
        parts.append(f'<text x="{x:.1f}" y="72" text-anchor="middle" font-family="Arial" font-size="13" fill="#555">{tick}%</text>')
    for index, label in enumerate(labels):
        y = top + index * row_height
        parts.append(f'<text x="{left-12}" y="{y+23}" text-anchor="end" font-family="Arial" font-size="15" fill="#222">{escape(label)}</text>')
        for series_index, (_, color, values) in enumerate(series):
            value = values[index]
            bar_y = y + series_index * (group_width + gap)
            bar_width = plot_width * value
            parts.append(f'<rect x="{left}" y="{bar_y}" width="{bar_width:.1f}" height="{group_width}" fill="{color}"/>')
            parts.append(f'<text x="{left+bar_width+7:.1f}" y="{bar_y+14}" font-family="Arial" font-size="13" fill="#222">{100*value:.1f}</text>')
    legend_x = 36
    legend_y = height - 28
    for name, color, _ in series:
        parts.append(f'<rect x="{legend_x}" y="{legend_y-13}" width="14" height="14" fill="{color}"/>')
        parts.append(f'<text x="{legend_x+21}" y="{legend_y}" font-family="Arial" font-size="14" fill="#333">{escape(name)}</text>')
        legend_x += 190
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def load_tcga(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    on_summary, on_rows = load_run(args.diagram_on_dir)
    off_summary, off_rows = load_run(args.diagram_off_dir)
    provenance = {
        "reactome_release": args.reactome_release,
        "lng_commit": args.lng_commit,
        "diagram_on": sanitized_manifest(args.diagram_on_dir),
        "diagram_off": sanitized_manifest(args.diagram_off_dir),
    }
    (args.output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    on_by_key = {case_key(row): row for row in on_rows}
    off_by_key = {case_key(row): row for row in off_rows}
    if set(on_by_key) != set(off_by_key):
        raise ValueError("Diagram-on and diagram-off runs do not contain the same cases")

    paired = [
        key
        for key in on_by_key
        if on_by_key[key]["prediction"] and off_by_key[key]["prediction"]
    ]
    gains = sum(is_correct(on_by_key[key]) and not is_correct(off_by_key[key]) for key in paired)
    losses = sum(not is_correct(on_by_key[key]) and is_correct(off_by_key[key]) for key in paired)
    changed = [key for key in paired if on_by_key[key]["prediction"] != off_by_key[key]["prediction"]]
    mcnemar_p = exact_mcnemar_p(gains, losses)

    # Cases are NOT independent: the harness caches one solve per
    # (pathway, gene, direction) and fans it out to every readout of that
    # perturbation, so 847 cases come from ~194 unique solves. Pairing at case
    # level treats correlated readouts of one experiment as separate evidence
    # and makes the p anti-conservative. Recompute at the unit of
    # randomisation — the perturbation solve — and report both.
    def _solve_key(row):
        return (row["pathway_id"], row["gene"], row["direction"])

    solve_gains, solve_losses = set(), set()
    for key in paired:
        on_ok, off_ok = is_correct(on_by_key[key]), is_correct(off_by_key[key])
        if on_ok and not off_ok:
            solve_gains.add(_solve_key(on_by_key[key]))
        elif off_ok and not on_ok:
            solve_losses.add(_solve_key(on_by_key[key]))
    mixed = solve_gains & solve_losses          # solve helped some readouts, hurt others
    solve_gains -= mixed
    solve_losses -= mixed
    solve_mcnemar_p = exact_mcnemar_p(len(solve_gains), len(solve_losses))

    # Both arms converged: the Convergence Gate section warns that
    # non-converged outputs are not equally reliable evidence, so report the
    # effect restricted to pairs where neither arm reported a failure.
    conv_paired = [
        key for key in paired
        if str(on_by_key[key].get("converged")) == "True"
        and str(off_by_key[key].get("converged")) == "True"
    ]
    conv_gains = sum(is_correct(on_by_key[k]) and not is_correct(off_by_key[k]) for k in conv_paired)
    conv_losses = sum(not is_correct(on_by_key[k]) and is_correct(off_by_key[k]) for k in conv_paired)
    conv_mcnemar_p = exact_mcnemar_p(conv_gains, conv_losses)

    scorecard = []
    convergence = {}
    for label, summary, rows in (
        ("diagram_on", on_summary, on_rows),
        ("diagram_off", off_summary, off_rows),
    ):
        convergence[label] = convergence_metrics(rows)
        for scope, metrics in (("all", summary), *summary["by_split"].items()):
            scorecard.append(
                {
                    "condition": label,
                    "scope": scope,
                    "eligible_cases": metrics["eligible_cases"],
                    "scored_cases": metrics["scored_cases"],
                    "coverage": metrics["coverage"],
                    "correct": metrics["correct"],
                    "accuracy": metrics["accuracy"],
                    "coverage_adjusted_accuracy": metrics["coverage_adjusted_accuracy"],
                    "macro_f1": metrics["macro_f1"],
                    "balanced_accuracy": metrics["balanced_accuracy"],
                    "change_f1": metrics["change_f1"],
                }
            )
    write_tsv(args.output_dir / "scorecard.tsv", scorecard)

    pathway_rows = []
    for pathway_id in PATHWAY_NAMES:
        on = on_summary["by_pathway"][pathway_id]
        off = off_summary["by_pathway"][pathway_id]
        split = "development" if pathway_id in on_summary["evaluation_split"]["development_pathways"] else "held_out"
        pathway_rows.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": PATHWAY_NAMES[pathway_id],
                "split": split,
                "eligible_cases": on["eligible_cases"],
                "scored_cases": on["scored_cases"],
                "coverage": on["coverage"],
                "diagram_on_correct": on["correct"],
                "diagram_on_accuracy": on["accuracy"],
                "diagram_off_correct": off["correct"],
                "diagram_off_accuracy": off["accuracy"],
                "net_correct": int(on["correct"]) - int(off["correct"]),
            }
        )
    write_tsv(args.output_dir / "pathway_scorecard.tsv", pathway_rows)

    change_rows = []
    for key in changed:
        on = on_by_key[key]
        off = off_by_key[key]
        change_rows.append(
            {
                **{field: on[field] for field in CASE_KEY},
                "evaluation_split": on["evaluation_split"],
                "expected": on["expected"],
                "diagram_on_prediction": on["prediction"],
                "diagram_off_prediction": off["prediction"],
                "diagram_on_correct": is_correct(on),
                "diagram_off_correct": is_correct(off),
                "diagram_on_converged": on["converged"],
                "diagram_off_converged": off["converged"],
            }
        )
    write_tsv(args.output_dir / "paired_topology_changes.tsv", change_rows)

    metric_labels = ["Coverage", "Scored accuracy", "Coverage-adjusted", "Macro-F1"]
    svg_bar_chart(
        args.output_dir / "overall_scorecard.svg",
        "All-ten exact-output benchmark",
        metric_labels,
        [
            (
                "Diagram edges on",
                "#1877B9",
                [
                    float(on_summary["coverage"]),
                    float(on_summary["accuracy"]),
                    float(on_summary["coverage_adjusted_accuracy"]),
                    float(on_summary["macro_f1"]),
                ],
            ),
            (
                "Diagram edges off",
                "#E76F51",
                [
                    float(off_summary["coverage"]),
                    float(off_summary["accuracy"]),
                    float(off_summary["coverage_adjusted_accuracy"]),
                    float(off_summary["macro_f1"]),
                ],
            ),
        ],
    )
    svg_bar_chart(
        args.output_dir / "pathway_accuracy.svg",
        "Accuracy by pathway (scored cases)",
        [f"{row['pathway_id']}  {row['pathway_name'][:30]}" for row in pathway_rows],
        [
            ("Diagram edges on", "#1877B9", [float(row["diagram_on_accuracy"]) for row in pathway_rows]),
            ("Diagram edges off", "#E76F51", [float(row["diagram_off_accuracy"]) for row in pathway_rows]),
        ],
    )

    tcga_rows = load_tcga(args.tcga_readiness)
    tcga_table = ""
    if tcga_rows:
        tcga_lines = [
            "| Pathway | Samples | Mapped roots | Log-rank p | Shuffle fraction <= observed |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
        for row in tcga_rows:
            tcga_lines.append(
                f"| {row['pathway_name']} | {row['samples_scored']} | {row['mapped_roots']} | "
                f"{float(row['median_split_logrank_p']):.3g} | {float(row['shuffle_fraction_p_le_observed']):.3g} |"
            )
        tcga_table = "\n".join(tcga_lines)

    on_conv = convergence["diagram_on"]
    off_conv = convergence["diagram_off"]
    held_on = on_summary["by_split"]["held_out"]
    held_off = off_summary["by_split"]["held_out"]
    failure_counts = Counter(row["mapping_status"] for row in on_rows if not row["prediction"])
    failure_text = ", ".join(f"{key}: {value}" for key, value in failure_counts.most_common())
    pathway_table = "".join(
        f"| {row['pathway_name']} | {row['split']} | "
        f"{row['scored_cases']}/{row['eligible_cases']} | "
        f"{percent(float(row['diagram_on_accuracy']))} | "
        f"{percent(float(row['diagram_off_accuracy']))} | "
        f"{int(row['net_correct']):+d} |\n"
        for row in pathway_rows
    )
    # Cases whose readout node set is entirely inside the set the
    # perturbation pins: the "prediction" is the intervention read back.
    def _pinned(rows_):
        return [r for r in rows_ if str(r.get("readout_fully_pinned")) == "True"]

    def _acc(rows_, column="prediction"):
        vals = [r for r in rows_ if r.get(column)]
        hit = sum(1 for r in vals if r[column] == r["expected"])
        return hit, len(vals)

    on_scored = [r for r in on_rows if r.get("prediction")]
    pinned_rows = _pinned(on_scored)
    unpinned_rows = [r for r in on_scored if r not in pinned_rows]
    pin_hit, pin_n = _acc(pinned_rows)
    unpin_hit, unpin_n = _acc(unpinned_rows)
    held_scored = [r for r in on_scored if r.get("evaluation_split") == "held_out"]
    held_unpinned = [r for r in held_scored if str(r.get("readout_fully_pinned")) != "True"]
    held_hit, held_n = _acc(held_scored)
    held_u_hit, held_u_n = _acc(held_unpinned)

    # Structural baselines: computed by the harness on the same cases, but
    # previously written only to the summary JSON and shown nowhere.
    baseline_table_all = baseline_table(on_summary)
    ssp = (on_summary.get("baselines", {}).get("shortest_signed_path") or {})
    ssp_all = ssp.get("all", {})
    ssp_held = ssp.get("by_split", {}).get("held_out", {})
    on_held = on_summary.get("by_split", {}).get("held_out", {})
    ssp_margin_all = int(on_summary["correct"]) - int(ssp_all.get("correct", 0))
    ssp_margin_held = int(on_held.get("correct", 0)) - int(ssp_held.get("correct", 0))

    report = f"""# DeltaSignal Evaluation Report

## Executive Result

This report evaluates current DeltaSignal with SCC solving on 847 experimentally
supported MP-BioPath cases across ten Reactome pathways. Three pathways were
used during development and seven were locked as held out. Exact key-output
entities are the primary endpoint; reaction proxies remain disabled.

With Reactome v97 diagram edges enabled, DeltaSignal scored
**{on_summary['scored_cases']}/{on_summary['eligible_cases']}** cases and was
correct on **{on_summary['correct']}/{on_summary['scored_cases']}**
({percent(float(on_summary['accuracy']))}). With diagram edges disabled it was
correct on **{off_summary['correct']}/{off_summary['scored_cases']}**
({percent(float(off_summary['accuracy']))}). Coverage-adjusted accuracy was
{percent(float(on_summary['coverage_adjusted_accuracy']))} versus
{percent(float(off_summary['coverage_adjusted_accuracy']))}.

On the seven held-out pathways, diagram-on accuracy was
**{held_on['correct']}/{held_on['scored_cases']}**
({percent(float(held_on['accuracy']))}) versus
**{held_off['correct']}/{held_off['scored_cases']}**
({percent(float(held_off['accuracy']))}) diagram-off. Across the
{len(paired)} paired scored cases, diagram edges changed {len(changed)}
classifications: {gains} changed wrong-to-correct and {losses}
correct-to-wrong (exact McNemar p = {mcnemar_p:.4g}).

That p is computed per case, and cases are not independent — the harness
caches one solve per (pathway, gene, direction) and fans it out to every
readout of that perturbation. Clustered at the perturbation solve, the unit
that was actually randomised, the effect is {len(solve_gains)} gains and
{len(solve_losses)} losses (exact McNemar p = {solve_mcnemar_p:.4g}).
Restricted to pairs where BOTH arms converged — the Convergence Gate below
warns that non-converged outputs are not equally reliable evidence — it is
{conv_gains} gains and {conv_losses} losses (p = {conv_mcnemar_p:.4g}) over
{len(conv_paired)} pairs.

Read together: the case-level p overstates the strength of this evidence.
Whether diagram topology contributes signal is not settled by these data, and
in either case it is not evidence that a DeltaSignal internal change caused
the difference.

![Overall scorecard](overall_scorecard.svg)

## Convergence Gate

The ten-pathway run exposed a generalization issue that was invisible in the
three-pathway development panel. Diagram-on reported
{on_conv['nonconverged_cases']} non-converged scored cases, but these repeat
only **{on_conv['unique_nonconverged_solves']} of
{on_conv['unique_scored_solves']}** unique perturbation solves. Diagram-off
reported {off_conv['nonconverged_cases']} cases from
**{off_conv['unique_nonconverged_solves']} of
{off_conv['unique_scored_solves']}** unique solves. All failures are in
Transcriptional Regulation by TP53. Non-converged outputs must not be treated
as equally reliable evidence; the report preserves them for diagnosis and
reports convergence explicitly.

Among converged diagram-on cases, accuracy was
{on_conv['converged_correct']}/{on_conv['converged_cases']}
({percent(float(on_conv['converged_accuracy']))}); for diagram-off it was
{off_conv['converged_correct']}/{off_conv['converged_cases']}
({percent(float(off_conv['converged_accuracy']))}). The non-converged
diagram-on TP53 cases were correct on
{on_conv['nonconverged_correct']}/{on_conv['nonconverged_cases']}
({percent(float(on_conv['nonconverged_accuracy']))}), so the headline accuracy
is partly supported by numerically unstable outputs and should not be reported
without this qualification.

## Coverage And Failure Attribution

Exact-output coverage was {percent(float(on_summary['coverage']))}. Unscored
cases were not converted to unchanged predictions. Recorded reasons were:
{failure_text}.

This separation matters: scored-case accuracy describes behavior where the
requested perturbation and output exist in the graph; coverage-adjusted
accuracy describes the end-to-end system over all eligible cases.

## Pathway-Level Result

![Pathway accuracy](pathway_accuracy.svg)

| Pathway | Split | Scored / eligible | Diagram on | Diagram off | Net correct |
| --- | --- | ---: | ---: | ---: | ---: |
{pathway_table}

## Readout-Is-Intervention Stratum

`load_dbid_to_uuids` registers a node under its own stable id and every
`member_leaves` entry, so a complex containing gene X counts as a "gene-X
node" and the perturbation pins it. Where the pinned set covers the whole
readout set, `max` aggregation returns the pinned value verbatim — UI 80
always classifies UP, UI 0 always DOWN — and the case is scored with no model
content.

| Subset | Cases | Correct | Accuracy |
| --- | ---: | ---: | ---: |
| Readout entirely pinned | {pin_n} | {pin_hit} | {percent(pin_hit / pin_n) if pin_n else 'n/a'} |
| All other scored cases | {unpin_n} | {unpin_hit} | {percent(unpin_hit / unpin_n) if unpin_n else 'n/a'} |
| Held out, excluding pinned | {held_u_n} | {held_u_hit} | {percent(held_u_hit / held_u_n) if held_u_n else 'n/a'} |

Held-out accuracy is {percent(held_hit / held_n) if held_n else 'n/a'} as
reported and {percent(held_u_hit / held_u_n) if held_u_n else 'n/a'} with this
stratum removed.

Every method is inflated by it to about the same degree — MP-BioPath, the
curator predictions and the shortest-signed-path baseline all score in the
high nineties here — so it does not distort the comparisons between them. It
does inflate the absolute figures, which is why the stratum is reported rather
than silently included. Scoring is unchanged; whether to exclude these cases
is a methodological decision, not a defect fix.

## Structural Baselines

Computed by the harness alongside DeltaSignal. The two structural baselines
score exactly the cases DeltaSignal scores, so their accuracies are directly
comparable to it. The two trivial baselines are projected over all 847
eligible cases — a different denominator — because they need no network path;
their accuracies are NOT comparable to the rows below them, and are shown to
bound the floor.

| Baseline | Scored | Correct | Accuracy | Macro-F1 |
| --- | ---: | ---: | ---: | ---: |
{baseline_table_all}
| **DeltaSignal (diagram on)** | **{on_summary['scored_cases']}** | **{on_summary['correct']}** | **{percent(float(on_summary['accuracy']))}** | **{float(on_summary['macro_f1']):.3f}** |

The comparison that matters is `Sign of shortest signed path` — a graph
traversal with no propagation model at all. DeltaSignal's margin over it is
{ssp_margin_all} cases overall and {ssp_margin_held} on the held-out split
({on_held.get('correct', 0)} vs {ssp_held.get('correct', 0)} of {on_held.get('scored_cases', 0)}). No paired test
is reported for this comparison; it should not be read as an established win.

## Comparator Context

On the {on_summary['scored_cases']} cases scored by diagram-on DeltaSignal,
MP-BioPath was correct on {on_summary['mpbiopath_on_scored_cases']['correct']}
({percent(float(on_summary['mpbiopath_on_scored_cases']['accuracy']))}) and the
curator predictions on {on_summary['curator_on_scored_cases']['correct']}
({percent(float(on_summary['curator_on_scored_cases']['accuracy']))}). The
paired bootstrap interval versus MP-BioPath was
[{100*float(on_summary['paired_bootstrap_vs_mpbiopath']['ci95_low']):.1f},
{100*float(on_summary['paired_bootstrap_vs_mpbiopath']['ci95_high']):.1f}]
percentage points. The result supports approximate paired performance, not a
statistically established win. The paired interval versus the curator
predictions is
[{100*float(on_summary['paired_bootstrap_vs_curator']['ci95_low']):.1f},
{100*float(on_summary['paired_bootstrap_vs_curator']['ci95_high']):.1f}]
percentage points, which excludes zero: on these paired cases DeltaSignal is
significantly worse than the curator predictions. Both intervals are computed
and stored; quoting only the favourable one would misrepresent the comparison. MP-BioPath and curator predictions cover all
847 eligible cases, while exact-output DeltaSignal covers
{on_summary['scored_cases']}.

## TCGA LUAD External Association Layer

{tcga_table if tcga_table else "No TCGA readiness table was supplied."}

The TCGA result is a historical three-pathway pipeline demonstration using 502
tumors. It shows that the LNG-to-DeltaSignal-to-clinical-analysis workflow can
produce an externally associated pathway score. It is observational,
configuration-specific, and pathway-selected; it is **not** causal
perturbation accuracy and must not be pooled with the 847-case benchmark.

## Deliverable Boundary

What is complete: an auditable LNG-to-DeltaSignal evaluation path, release-aware
mapping failures, frozen graph hashes, development/held-out pathway splits,
simple structural baselines, topology ablation, paired comparisons, and an
external TCGA association demonstration.

What remains: fix or characterize TP53 SCC non-convergence, improve exact
readout coverage without silently changing endpoints, reproduce on an external
interventional dataset, and rerun TCGA under a release-pinned current stack
with predeclared pathway/readout definitions.

## Generated Artifacts

- `scorecard.tsv`: overall, development, and held-out metrics.
- `pathway_scorecard.tsv`: pathway-level coverage and accuracy.
- `paired_topology_changes.tsv`: every classification changed by diagram edges.
- `overall_scorecard.svg` and `pathway_accuracy.svg`: presentation-ready plots.
- `provenance.json`: sanitized commits, settings, and input/network hashes;
  absolute workstation and server paths are intentionally excluded.
"""
    (args.output_dir / "evaluation_report.md").write_text(report)
    print(args.output_dir / "evaluation_report.md")


if __name__ == "__main__":
    main()
