"""MultiQC submodule to plot BBSplit alignment distribution"""

from multiqc.plots import bargraph


def plot_bbsplit(samples, file_type, plot_title, plot_params):
    """Create a stacked bar chart showing read distribution across reference genomes"""
    
    # Prepare data structure for bargraph
    data = {}
    cats = {}
    
    for s_name, sample in samples.items():
        data[s_name] = {}
        
        for ref_name, values in sample["data"].items():
            # Values array: [%unambiguousReads, unambiguousMB, %ambiguousReads, ambiguousMB,
            #                unambiguousReads, ambiguousReads, assignedReads, assignedBases]
            assigned_reads = values[6]  # assignedReads column
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
        "cpswitch_c_active": False,  # Start with percentage view
    }
    for key, value in plot_params.items():
        setattr(pconfig, key, value)
    
    return bargraph.plot(data, cats, pconfig)
