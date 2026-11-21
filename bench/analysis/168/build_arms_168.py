#!/usr/bin/env python3
"""specs/168: build arms from a FullTextMining bundle with the specs/166 machinery.

Arms per bundle: `with` (as delivered), `without` (literature rows and nodes removed) and, when asked,
`shuf<N>` (the specs/166 degree- and sign-preserving rewiring). Each arm is a catalog build directory
`<prefix>-<arm>`; an existing directory is refused, never overwritten.

Usage: build_arms_168.py SRC_DIR BUILDS_ROOT PREFIX ARM [ARM ...]     (ARM: with | without | shufNN)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "166"))
import build_arms_166 as b  # noqa: E402


def build(src_root, builds, prefix, arm):
    pathways = sorted(p.name for p in src_root.iterdir() if (p / "logic_network.csv").exists())
    arm_dir = builds / f"{prefix}-{arm}"
    if arm_dir.exists():
        raise SystemExit(f"{arm_dir} exists; refusing to overwrite")
    arm_dir.mkdir(parents=True)
    handled = []
    for p in pathways:
        if arm == "with":
            b.copy_pathway(src_root / p, arm_dir / p)
        elif arm == "without":
            handled.append(b.make_without(src_root / p, arm_dir / p))
        elif arm.startswith("shuf"):
            handled.append(b.make_shuffle(src_root / p, arm_dir / p, int(arm[4:]))[0])
        else:
            raise SystemExit(f"unknown arm {arm}")
    b.build_json(arm_dir, f"specs/168 arm {arm} of {src_root.name}", pathways)
    return arm_dir, len(pathways), sum(handled)


def main():
    src, builds, prefix = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    for arm in sys.argv[4:]:
        d, n, h = build(src, builds, prefix, arm)
        print(f"{d.name}: {n} pathways" + (f"; literature rows handled {h}" if h else ""))


if __name__ == "__main__":
    main()
