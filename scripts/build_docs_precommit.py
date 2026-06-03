"""Run pdoc and stage pdoc_html for commit. Used by pre-commit hook."""

import subprocess
import sys

subprocess.run(
    [
        sys.executable,
        "-m",
        "pdoc",
        "--html",
        "--output-dir",
        "pdoc_html",
        "--force",
        "sprime",
    ],
    check=True,
)
subprocess.run(["git", "add", "pdoc_html"], check=True)
