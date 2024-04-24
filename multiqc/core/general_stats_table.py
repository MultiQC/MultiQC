import re

from multiqc.plots import table
from multiqc.utils import config, report


def _general_stats_table() -> None:
    """
    Construct HTML for the general stats table.
    """

    # Remove empty data sections from the General Stats table
    empty_keys = [i for i, d in enumerate(report.general_stats_data[:]) if len(d) == 0]
    empty_keys.sort(reverse=True)
    for i in empty_keys:
        del report.general_stats_data[i]
        del report.general_stats_headers[i]

    # Add general-stats IDs to table row headers
    for idx, h in enumerate(report.general_stats_headers):
        for k in h.keys():
            if "rid" not in h[k]:
                h[k]["rid"] = re.sub(r"\W+", "_", k).strip().strip("_")
            ns_html = re.sub(r"\W+", "_", h[k]["namespace"]).strip().strip("_").lower()
            report.general_stats_headers[idx][k]["rid"] = report.save_htmlid(
                f"mqc-generalstats-{ns_html}-{h[k]['rid']}"
            )

    all_hidden = True
    for headers in report.general_stats_headers:
        for h in headers.values():
            if not h.get("hidden", False):
                all_hidden = False
                break

    # Generate the General Statistics HTML & write to file
    if len(report.general_stats_data) > 0 and not config.skip_generalstats and not all_hidden:
        pconfig = {
            "id": "general_stats_table",
            "table_title": "General Statistics",
            "save_file": True,
            "raw_data_fn": "multiqc_general_stats",
        }
        p = table.plot(report.general_stats_data, report.general_stats_headers, pconfig)
        report.general_stats_html = p.add_to_report(clean_html_id=False)
    else:
        config.skip_generalstats = True
