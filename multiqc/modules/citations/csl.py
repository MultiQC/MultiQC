"""Parse CSL-JSON tool citations into the shared Citation model.

The canonical citations input is
[CSL-JSON](https://citeproc-js.readthedocs.io/en/latest/csl-json/): a JSON
array of CSL items (a JSON object keyed by id, or a single item, are also
accepted). Each item carries the bibliographic fields plus an optional `custom`
object holding the tool name and the version that actually ran:

    {
      "id": "fastqc",
      "type": "software",
      "title": "FastQC: A Quality Control Tool ...",
      "author": [{"family": "Andrews", "given": "S"}],
      "issued": {"date-parts": [[2010]]},
      "container-title": "Bioinformatics",
      "DOI": "10.1093/bioinformatics/btw354",
      "URL": "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
      "custom": {"tool": "fastqc", "version": "0.12.1"}
    }

This is the CSL adapter: it loads the JSON, detects its shape, and emits
`Citation` objects via `Citation.from_csl`.
"""

import json
from typing import List

from .citation import Citation


def parse_csl(content: str, fn: str = "<citations>") -> List[Citation]:
    """Parse CSL-JSON file content into Citations.

    Accepts a JSON array of CSL items (canonical), a JSON object keyed by id, or
    a single CSL item object. Raises ValueError with the file path on malformed
    JSON rather than silently producing an empty report.
    """
    try:
        data = json.loads(content)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Could not parse CSL-JSON citations file '{fn}': {exc}") from exc

    if isinstance(data, list):
        items = data
    elif isinstance(data, dict):
        # A single CSL item carries identifying keys at the top level; anything
        # else is treated as a mapping of {id: item}.
        if {"id", "custom", "type", "title"} & set(data.keys()):
            items = [data]
        else:
            items = list(data.values())
    else:
        raise ValueError(f"Unexpected top-level JSON type in citations file '{fn}': {type(data).__name__}")

    return [Citation.from_csl(item) for item in items]
