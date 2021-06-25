from __future__ import print_function
from io import open
import argparse
import logging
import re
import base64
import hashlib
from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="Checks inline scripts in report.html against a " "whitelist of hashes (Content Security Policy)"
)

parser.add_argument("--report", help="Report file (.html)", required=True)
parser.add_argument("--whitelist", help="Whitelist to compare against", required=True)
args = parser.parse_args()


def is_executable(script):
    executable_types = {
        "text/javascript",
        "application/javascript",
        "module" "text/ecmascript",
        "application/ecmascript",
    }
    script_type = script.attrs.get("type", "text/javascript")
    return script_type in executable_types


def get_hash(script):
    return "sha256-" + base64.b64encode(hashlib.sha256(script.encode("utf-8")).digest()).decode()


whitelist_with_comments = open(args.whitelist, "r", encoding="utf-8").read()
whitelist = re.sub(r"#.*", "", whitelist_with_comments)

html_report = open(args.report, "r", encoding="utf-8").read()
soup = BeautifulSoup(html_report, features="html.parser")

scripts = [script.get_text() for script in soup.select("script") if is_executable(script)]
missing_scripts = [script for script in scripts if get_hash(script) not in whitelist]

if missing_scripts:
    logger.warning("The following scripts are missing from {}".format(args.whitelist))
    for script in missing_scripts:
        hash = get_hash(script)
        snippet = script.replace("\n", "")[0:80]
        logger.warning("  '{}' # {}".format(hash, snippet))
    exit(1)
