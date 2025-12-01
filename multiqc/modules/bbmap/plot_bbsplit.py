"""MultiQC submodule to plot BBSplit alignment distribution"""

from multiqc.plots import bargraph

# BBSplit stats file column indices
# Columns: %unambiguousReads, unambiguousMB, %ambiguousReads, ambiguousMB,
#          unambiguousReads, ambiguousReads, assignedReads, assignedBases
BBSPLIT_COLS = {
    "pct_unambiguous": 0,
    "unambiguous_mb": 1,
    "pct_ambiguous": 2,
    "ambiguous_mb": 3,
    "unambiguous_reads": 4,
    "ambiguous_reads": 5,
    "assigned_reads": 6,
    "assigned_bases": 7,
}


def plot_bbsplit(samples, file_type, plot_title, plot_params):
    """Create a stacked bar chart showing read distribution across reference genomes"""

    # Prepare data structure for bargraph
    data = {}
    cats = {}

    for s_name, sample in samples.items():
        data[s_name] = {}

        for ref_name, values in sample["data"].items():
            # Use BBSPLIT_COLS dictionary for array access
            assigned_reads = values[BBSPLIT_COLS["assigned_reads"]]
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
        "cpswitch_c_active": True,  # Start with counts view (True = counts active)
    }
    pconfig.update(plot_params)

    return bargraph.plot(data, cats, pconfig)
