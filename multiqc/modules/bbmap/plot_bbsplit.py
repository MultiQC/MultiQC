"""MultiQC submodule to plot BBSplit alignment distribution"""

from multiqc.plots import bargraph

# BBSplit stats file column indices
BBSPLIT_COLS = [
    "pct_unambiguous",
    "unambiguous_mb",
    "pct_ambiguous",
    "ambiguous_mb",
    "unambiguous_reads",
    "ambiguous_reads",
    "assigned_reads",
    "assigned_bases",
]


def plot_bbsplit(samples, file_type, plot_title, plot_params):
    """Create a stacked bar chart showing read distribution across reference genomes"""

    # Prepare data structure for bargraph
    data = {}
    cats = {}

    for s_name, sample in samples.items():
        data[s_name] = {}

        for ref_name, values in sample["data"].items():
            assigned_reads = values[BBSPLIT_COLS.index("assigned_reads")]
            data[s_name][ref_name] = assigned_reads

            # Define category (reference genome) if not already done
            if ref_name not in cats:
                cats[ref_name] = {"name": ref_name}

    # Configure the plot
    pconfig = {
        "id": "bbmap-" + file_type + "_plot",
        "title": "BBTools: " + plot_title,
        "ylab": "Number of Reads",
        "cpswitch_counts_label": "Number of Reads",
        "cpswitch_percent_label": "Percentage of Reads",
    }
    pconfig.update(plot_params)

    return bargraph.plot(data, cats, pconfig)
