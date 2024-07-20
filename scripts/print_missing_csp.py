"""
Print missing hashes for a Content Security Policy (CSP) whitelist (see https://github.com/MultiQC/MultiQC/pull/911)
for the scripts inlined in a MultiQC report. Usage:

multiqc test_data/data/modules --filename full_report.html
pip install beautifulsoup4
python scripts/print_missing_csp.py --report full_report.html [--current csp.txt]

If the current list is provided, will print only missing hashes.
"""

import argparse
import base64
import hashlib
import re

from bs4 import BeautifulSoup


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


def get_hash(script: str) -> str:
    return "sha256-" + base64.b64encode(hashlib.sha256(script.encode("utf-8")).digest()).decode()


def main():
    parser = argparse.ArgumentParser(
        description="Checks inline scripts in report.html against a whitelist of hashes (Content Security Policy)"
    )
    parser.add_argument("--report", help="Report file (.html)", required=True)
    parser.add_argument("--current", help="Current whitelist", required=False)
    args = parser.parse_args()

    whitelist: str = ""
    if args.current:
        with open(args.current) as f:
            whitelist_with_comments: str = f.read()
        whitelist = re.sub(r"#.*", "", whitelist_with_comments)

    with open(args.report) as f:
        soup = BeautifulSoup(f.read(), features="html.parser")
    scripts = [script.get_text() for script in soup.select("script") if is_executable(script)]

    for script in scripts:
        if get_hash(script) not in whitelist:
            snippet = script.replace("\n", "")[0:80]
            print("  '{}' # {}".format(get_hash(script), snippet))


if __name__ == "__main__":
    main()
