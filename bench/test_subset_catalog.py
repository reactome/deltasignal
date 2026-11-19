"""specs/166: a build of a few pathways is a subset catalog, so the liveness guards are skipped."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402

LIST = Path(__file__).resolve().parent / "catalog_pathways.tsv"


def write(tmp, doc):
    p = tmp / "BUILD.json"
    p.write_text(json.dumps(doc))
    return p


def test_ten_listed_pathways_is_a_subset(tmp_path):
    assert bm.subset_catalog(write(tmp_path, {"pathways_requested": [f"R-HSA-{i}" for i in range(10)]}), LIST)


def test_full_catalog_count_is_not_a_subset(tmp_path):
    full = sum(1 for line in open(LIST) if line.startswith("R-HSA-"))
    assert not bm.subset_catalog(write(tmp_path, {"pathways_requested": full}), LIST)


def test_unreadable_or_silent_build_is_not_a_subset(tmp_path):
    assert not bm.subset_catalog(tmp_path / "missing.json", LIST)
    assert not bm.subset_catalog(write(tmp_path, {}), LIST)
