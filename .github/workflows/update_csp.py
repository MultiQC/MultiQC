import argparse
import base64
import hashlib
import logging
from io import open

from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="Checks inline scripts in report.html against a whitelist of hashes (Content Security Policy)"
)

parser.add_argument("--report", help="Report file (.html)", required=True)
parser.add_argument("--whitelist", help="Whitelist to update", required=True)
args = parser.parse_args()


def is_executable(script):
    executable_types = {
        "text/javascript",
        "application/javascript",
        "module",
        "text/ecmascript",
        "application/ecmascript",
    }
    script_type = script.attrs.get("type", "text/javascript")
    return script_type in executable_types


def get_hash(script):
    return "sha256-" + base64.b64encode(hashlib.sha256(script.encode("utf-8")).digest()).decode()


current_hash_by_snippet = dict()
with open(args.whitelist, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if line.startswith("'sha256-"):
            js_hash = line.split("#")[0].strip().strip("'")
            snippet = line.split("#")[1].strip()
            current_hash_by_snippet[snippet] = js_hash

html_report = open(args.report, "r", encoding="utf-8").read()
soup = BeautifulSoup(html_report, features="html.parser")
html_scripts = [script.get_text() for script in soup.select("script") if is_executable(script)]
new_hash_by_snippet = dict()
for script in html_scripts:
    js_hash = get_hash(script)
    snippet = script.replace("\n", "")[0:80].strip()
    new_hash_by_snippet[snippet] = js_hash

# Write updated whitelist
with open(args.whitelist, "w", encoding="utf-8") as f:
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
