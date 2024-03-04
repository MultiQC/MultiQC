"""
Update CSP.txt whitelist of JavaScript files embedded in MultiQC reports
"""

import base64
import hashlib
import os
from pathlib import Path

workspace_path = Path(os.environ.get("GITHUB_WORKSPACE", ""))
whitelist_path = workspace_path / "CSP.txt"


def get_hash(script: str) -> str:
    return "sha256-" + base64.b64encode(hashlib.sha256(script.encode("utf-8")).digest()).decode()


# Read current whitelist
current_hash_by_snippet = dict()
with whitelist_path.open() as f:
    for line in f:
        line = line.strip()
        if line.startswith("'sha256-"):
            js_hash = line.split("#")[0].strip().strip("'")
            snippet = line.split("#")[1].strip()
            current_hash_by_snippet[snippet] = js_hash

# Find all JS scripts in the codebase
dirs_to_scan = [
    workspace_path / "multiqc" / "templates",
    workspace_path / "multiqc" / "modules",
]

new_hash_by_snippet = dict()
for dir_to_scan in dirs_to_scan:
    for root, dirs, files in os.walk(dir_to_scan):
        for file in files:
            if file.endswith(".js"):
                path = Path(root) / file
                if path.with_suffix(".min.js").exists():
                    path = path.with_suffix(".min.js")
                with path.open() as f:
                    contents = f.read()
                    js_hash = get_hash(contents)
                    new_hash_by_snippet[path.relative_to(workspace_path / "multiqc")] = js_hash

# Write updated whitelist
with whitelist_path.open("w") as f:
    whitelist = "\n    ".join(f"'{js_hash}' # {snippet}" for snippet, js_hash in new_hash_by_snippet.items())
    f.write(
        f"""\
# If you are hosting MultiQC reports >= v1.8 on a web application with CSP
# (Content Security Policy), you will need the following scripts whitelisted:

script-src 'self'
    {whitelist}
    'unsafe-eval'
;"""
    )
