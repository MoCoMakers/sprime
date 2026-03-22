#!/usr/bin/env python3
"""
Generate README-PyPI.md from README.md for the PyPI project description.

PyPI does not render GitHub-style fenced math or LaTeX; this script replaces the
fenced math block with plain text that reads well on pypi.org.

Machine-readable sync line (line 1 of README-PyPI.md):
  [//]: # pypi-readme-sync FORMULA_REVISION:N LAST_VALIDATED_WITH_PYPI_RELEASE:vX.Y.Z

- FORMULA_REVISION: bump in this script when changing PYPI_MATH_REPLACEMENT below.
- LAST_VALIDATED_WITH_PYPI_RELEASE: last release whose PyPI long description was
  validated; update with: python scripts/sync_readme_pypi.py --set-release-tag v0.3.0

Usage:
  python scripts/sync_readme_pypi.py              # write README-PyPI.md
  python scripts/sync_readme_pypi.py --check      # exit 1 if README-PyPI.md is stale
  python scripts/sync_readme_pypi.py --set-release-tag v0.3.0
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
README_MD = REPO_ROOT / "README.md"
README_PYPI = REPO_ROOT / "README-PyPI.md"

# Bump when the PyPI-specific math/plain-text template changes (not when README LaTeX changes).
FORMULA_REVISION = 1

# Single fenced math block in README.md (GitHub renders this).
_MATH_FENCE = re.compile(r"```math\n.*?```\n?", re.DOTALL)

# Replacement for PyPI (no LaTeX). Link to GitHub for rendered math.
PYPI_MATH_REPLACEMENT = """**Equivalent logarithmic form** — LaTeX renders on [GitHub README](https://github.com/MoCoMakers/sprime/blob/main/README.md#about-s-s-prime). On PyPI (plain text):

`S' = ln( (Zero_asymptote - Inf_asymptote) / EC50 + sqrt( ((Zero_asymptote - Inf_asymptote) / EC50)^2 + 1 ) )`

"""

_SYNC_LINE_RE = re.compile(
    r"^\[//\]: # pypi-readme-sync FORMULA_REVISION:(\d+) LAST_VALIDATED_WITH_PYPI_RELEASE:(\S+)\s*$",
    re.MULTILINE,
)


def _parse_sync_meta(text: str) -> Tuple[Optional[int], Optional[str]]:
    m = _SYNC_LINE_RE.search(text)
    if not m:
        return None, None
    return int(m.group(1)), m.group(2)


def _build_sync_line(release_tag: str) -> str:
    return (
        f"[//]: # pypi-readme-sync FORMULA_REVISION:{FORMULA_REVISION} "
        f"LAST_VALIDATED_WITH_PYPI_RELEASE:{release_tag}\n\n"
    )


def _default_release_tag() -> str:
    return "v0.2.1"


def _read_release_tag_from_existing() -> str:
    if not README_PYPI.is_file():
        return _default_release_tag()
    body = README_PYPI.read_text(encoding="utf-8")
    _, tag = _parse_sync_meta(body)
    return tag if tag else _default_release_tag()


def generate_content(release_tag: str) -> str:
    source = README_MD.read_text(encoding="utf-8")
    if not _MATH_FENCE.search(source):
        print(
            "error: README.md has no ```math ... ``` block; sync script needs updating.",
            file=sys.stderr,
        )
        sys.exit(2)
    replaced = _MATH_FENCE.sub(PYPI_MATH_REPLACEMENT, source, count=1)
    if "```math" in replaced:
        print("error: multiple ```math blocks or replacement failed.", file=sys.stderr)
        sys.exit(2)
    return _build_sync_line(release_tag) + replaced


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit 1 if README-PyPI.md does not match generator output.",
    )
    parser.add_argument(
        "--set-release-tag",
        metavar="TAG",
        help="Set LAST_VALIDATED_WITH_PYPI_RELEASE (e.g. v0.3.0) and regenerate.",
    )
    args = parser.parse_args()

    if not README_MD.is_file():
        print(f"error: missing {README_MD}", file=sys.stderr)
        sys.exit(2)

    if args.set_release_tag:
        tag = args.set_release_tag.strip()
        if not tag.startswith("v"):
            print("warning: release tags usually look like v0.3.0", file=sys.stderr)
    else:
        tag = _read_release_tag_from_existing()

    content = generate_content(tag)

    if args.check:
        if not README_PYPI.is_file():
            print(
                f"error: {README_PYPI} missing; run: python scripts/sync_readme_pypi.py",
                file=sys.stderr,
            )
            sys.exit(1)
        existing = README_PYPI.read_text(encoding="utf-8")
        if existing.replace("\r\n", "\n") != content.replace("\r\n", "\n"):
            print(
                "README-PyPI.md is out of sync with README.md (or release tag metadata).\n"
                "Run: python scripts/sync_readme_pypi.py\n"
                "If you released to PyPI, optionally:\n"
                "  python scripts/sync_readme_pypi.py --set-release-tag vX.Y.Z",
                file=sys.stderr,
            )
            sys.exit(1)
        return

    README_PYPI.write_text(content, encoding="utf-8", newline="\n")
    print(f"Wrote {README_PYPI} (LAST_VALIDATED_WITH_PYPI_RELEASE={tag})")


if __name__ == "__main__":
    main()
