"""Parse riker `error` outputs (error-mismatch.txt, error-indel.txt, error-overlap.txt).

Each file has one row per (stratifier, covariate) combination, e.g. error rate by base 
quality ('bq') or by read cycle. The row where both are 'all' holds the overall metrics for
the sample, and is one reported in general stats. 

"""

import logging
from typing import Dict, Optional
from .util import read_tsv, to_float, to_int
log = logging.getLogger(__name__)


def parse_reports(module):
    parsed_samples = set()
    parsed_samples |= _parse_mismatch(module)
    # TODO: parse error-indel.txt and error-overlap.txt
    return parsed_samples


def _is_overall_row(row) -> bool:
    """The `all / all` row holds the overall (unstratified) metrics for a sample."""
    return row.get("stratifier") == "all" and row.get("covariate") == "all"


def _parse_mismatch(module) -> set:
    data_by_sample: Dict[str, Dict[str, Optional[float]]] = {}

    for f in module.find_log_files("riker/error_mismatch", filehandles=True):
        for row in read_tsv(f["f"], source=f["fn"]):
            sample = row.get("sample")
            #skip the stratified rows, only the overall row is used for the General Stats table. 
            if not sample or not _is_overall_row(row):
                continue
            s_name = module.clean_s_name(sample, f)
            try:
                data_by_sample[s_name] = {
                    "total_bases": to_int(row["total_bases"]),
                    "error_bases": to_int(row["error_bases"]),
                    "frac_error": to_float(row["frac_error"]),
                    "q_score": to_float(row["q_score"]),
                }
            except (KeyError, TypeError, ValueError) as e:
                log.warning(f"riker: skipping row in {f['fn']} for sample {sample}: {e}")
                continue
            module.add_data_source(f, s_name, section="error_mismatch")

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    headers = {
        "frac_error": {
            "title": "Mismatch rate",
            "description": "Percentage of aligned bases that mismatch the reference (riker error)",
            "min": 0,
            "suffix": "%",
            # riker reports a fraction and shows it as a percentage
            "modify": lambda x: x * 100 if x is not None else None,
            "format": "{:,.3f}",
            "scale": "OrRd",
        },
        "q_score": {
            "title": "Mismatch Q",
            "description": "Mismatch error rate as a Phred-scaled quality score (higher is better)",
            "min": 0,
            "format": "{:,.1f}",
            "scale": "RdYlGn",
            # same information as the mismatch rate, so it is hidden by default.
            "hidden": True,
        },
    }
    module.general_stats_addcols(data_by_sample, headers, namespace="error")
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_error_mismatch")
    return set(data_by_sample.keys())