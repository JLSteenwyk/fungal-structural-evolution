#!/usr/bin/env python3
"""Fail when Markdown links point to missing local files.

External URLs and same-document anchors are intentionally outside this check.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path


LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")


def local_target(target: str) -> str | None:
    """Return a local file component, excluding URLs and anchor-only links."""
    target = target.strip().split(maxsplit=1)[0].strip("<>")
    if not target or target.startswith("#") or "://" in target or target.startswith("mailto:"):
        return None
    return target.split("#", 1)[0] or None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("documents", nargs="+", type=Path, help="Markdown documents to check")
    args = parser.parse_args()
    missing: list[tuple[Path, str]] = []
    checked = 0
    for document in args.documents:
        if not document.is_file():
            missing.append((document, "document itself is missing"))
            continue
        for raw_target in LINK.findall(document.read_text()):
            target = local_target(raw_target)
            if target is None:
                continue
            checked += 1
            if not (document.parent / target).exists():
                missing.append((document, raw_target))
    print(f"checked_local_links={checked}")
    print(f"missing_local_links={len(missing)}")
    for document, target in missing:
        print(f"{document}: {target}")
    if missing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
