#!/usr/bin/env python

""" MultiQC module to parse output from odgi stats """

from __future__ import print_function
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

import yaml

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """graph
    The MultiQC module to parse and plot odgi stats output.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="ODGI",
            anchor="odgi",
            href="https://github.com/pangenome/odgi",
            info="is an optimized dynamic graph/genome implementation.",
        )

        self.odgi_stats_map = dict()

        # Find and load any odgi stats file
        log_files = self.find_log_files("odgi", filehandles=True)

        # Parse odgi stats data
        for f in log_files:
            self.parse_odgi_stats_report(f)
            self.add_data_source(f)

        self.odgi_stats_map = self.ignore_samples(self.odgi_stats_map)

        if not self.odgi_stats_map:
            raise UserWarning

        if len(self.odgi_stats_map) == 1:
            log.info("Found {} report".format(len(self.odgi_stats_map)))
        else:
            log.info("Found {} reports".format(len(self.odgi_stats_map)))

        # Write parsed report data to a file
        self.write_data_file(self.odgi_stats_map, "multiqc_odgi_stats")

        # Add the general section containing the general odgi stats
        self.plot_general_odgi_stats()

        # Plot the odgi stats metrics as a lineplot
        self.plot_odgi_metrics()

    def parse_odgi_stats_report(self, f):
        """
        Load the odgi stats YAML file.
        The stats reports assume the following YAML structure:

        ---
        length: 206263
        nodes: 3751
        edges: 5195
        paths: 13
        num_weakly_connected_components: 1
        weakly_connected_components:
          - component:
              id: 0
              nodes: 3751
              is_acyclic: 'no'
        num_nodes_self_loops:
          total: 0
          unique: 0
        A: 57554
        C: 43275
        G: 41944
        T: 63490
        mean_links_length:
          - length:
              path: all_paths
              in_node_space: 1.64973
              in_nucleotide_space: 7.34035
              num_links_considered: 202793
              num_gap_links_not_penalized: 147940
        sum_of_path_node_distances:
          - distance:
              path: all_paths
              in_node_space: 5.53383
              in_nucleotide_space: 2.1454
              nodes: 202806
              nucleotides: 3757597
              num_penalties: 231
              num_penalties_different_orientation: 0
        """
        data = {}
        try:
            data = yaml.load(f["f"], Loader=yaml.SafeLoader)
        except Exception as e:
            log.warning("Could not parse YAML for '{}': \n  {}".format(f, e))
        self.odgi_stats_map.update(self.compress_stats_data(self.extract_sample_name(f["fn"]), data))

    def plot_general_odgi_stats(self):
        """
        Plot the odgi stats by adding them to the general statistics table.
        """
        data = dict()
        headers = OrderedDict()
        headers["Length"] = {"title": "Length", "description": "Graph length in nucleotides.", "scale": "BuPu"}
        headers["Nodes"] = {
            "title": "Nodes",
            "description": "Number of nodes in the graph.",
            "scale": "OrRd",
            'format': '{:,.0f}', 
        }
        headers["Edges"] = {"title": "Edges", "description": "Number of edges in the graph.", "scale": "PuBu"}
        headers["Paths"] = {"title": "Paths", "description": "Number of paths in the graph.", "scale": "Greens"}
        headers["Components"] = {
            "title": "Components",
            "description": "Number of weakly connected components in the " "graph.",
            "scale": "Oranges",
        }
        headers["A"] = {"title": "A", "description": "Number of adenine bases in the graph.", "scale": "Spectral"}
        headers["C"] = {"title": "C", "description": "Number of cytosine bases in the graph.", "scale": "Greys"}
        headers["T"] = {"title": "T", "description": "Number of thymine bases in the graph.", "scale": "Blues"}
        headers["G"] = {"title": "G", "description": "Number of guanine bases in the graph.", "scale": "BrBG"}
        headers["N"] = {"title": "N", "description": "Number of `N` basis in the graph.", "scale": "BrBG"}
        headers["Total"] = {
            "title": "Nodes With Self Loops - Total",
            "description": "Total number of nodes having self loops in the " "graph.",
            "scale": "Set1",
            "hidden": True,
        }
        headers["Unique"] = {
            "title": "Nodes With Self Loops - Unique",
            "description": "Number of unique nodes having self loops in the " "graph.",
            "scale": "Set2",
            "hidden": True,
        }
        for fn in self.odgi_stats_map.keys():
            file_stats = self.odgi_stats_map[fn]
            data.update({fn: file_stats["General stats {}".format(fn)]})
        self.general_stats_addcols(data, headers)

    def plot_odgi_metrics(self):
        """
        Plot odgi metrics as a lineplot supporting raw counts and log10 counts.
        """
        mean_links_length_nodes_space = dict()
        mean_links_length_nucleotide_space = dict()
        sum_of_path_nodes_distances_nodes_space = dict()
        sum_of_path_nodes_distances_nucleotide_space = dict()
        odgi_stats_file_names = sorted(self.odgi_stats_map.keys())
        seq_smooth = odgi_stats_file_names[-2:]
        sorted_filenames = seq_smooth + odgi_stats_file_names[:-2]
        for sample_name in sorted_filenames:
            file_stats = self.odgi_stats_map[sample_name]
            mean_links_length_nodes_space.update(
                {sample_name: float(file_stats["Mean_links_length {}".format(sample_name)]["In_node_space"])}
            )
            mean_links_length_nucleotide_space.update(
                {sample_name: float(file_stats["Mean_links_length {}".format(sample_name)]["In_nucleotide_space"])}
            )
            sum_of_path_nodes_distances_nodes_space.update(
                {sample_name: float(file_stats["Sum_of_path_nodes_distances {}".format(sample_name)]["In_node_space"])}
            )
            sum_of_path_nodes_distances_nucleotide_space.update(
                {
                    sample_name: float(
                        file_stats["Sum_of_path_nodes_distances {}".format(sample_name)]["In_nucleotide_space"]
                    )
                }
            )

        metrics_lineplot_config = {
            "id": "odgi_metrics_lineplot",
            "ylab": "Count",
            "categories": True,
            "yDecimals": False,
            "logswitch": True,
            "title": "ODGI: ODGI stats metrics - Evaluating sorting goodness - How linear is a graph?",
        }

        self.add_section(
            name="ODGI metrics",
            anchor="odgi_stats",
            description="The ODGI metrics section",
            helptext="<b>sum-of-path-node-distances</b>:<br>For each path we iterate from node to node and count the node distance of nodes on the pangenome "
            "level normalized by the path length. If a node is reversed, we count the node distance twice.<br><br>"
            "<b>sum-of-path-nucleotide-distances</b>:<br>"
            "For each path we iterate from node to node and count the nucleotide distance of nodes on the pangenome level normalized by the path "
            "length. If the node is reversed, we count the nucleotide distance twice.<br><br>"
            "<b>mean-links-length</b>:<br>In contrast to the <em>sum-of-path-distances</em>, for each path we iterate from node to node and count the node distance <b>mean_links_length_nodes_space</b> or nucleotide distance <b>mean_links_length_nucleotide_space</b> of nodes within the same path only! "
            "We then normalized by the path length.",
            plot=linegraph.plot(
                {
                    "mean_links_length_nodes_space": mean_links_length_nodes_space,
                    "mean_links_length_nucleotide_space": mean_links_length_nucleotide_space,
                    "sum_of_path_nodes_distances_nodes_space": sum_of_path_nodes_distances_nodes_space,
                    "sum_of_path_nodes_distances_nucleotide_space": sum_of_path_nodes_distances_nucleotide_space,
                },
                pconfig=metrics_lineplot_config,
            ),
        )

    def extract_sample_name(self, file_name):
        """
        Extracts the sample name from a given file name.
        Expects and returns one of seqwish, smooth or cons@*.
        If none of the above keywords are present, the full filename except for the ending is returned.
        """
        # TODO if non of the above, just take the full name
        if "cons@" in file_name:
            file_name = file_name.split(".")
            consensus_identifier = list((e for e in file_name if "cons@" in e))
            return consensus_identifier[0]
        elif "smooth" in file_name:
            return "smooth"
        elif "seqwish" in file_name:
            return "seqwish"
        else:
            return ".".join(file_name.split(".")[:-3])

    def compress_stats_data(self, sample_name, data) -> dict:
        """
        Compress odgi stats into a single dictionary to visualize.
        """
        mean_links_length = data["mean_links_length"]
        # we have to find the entry with path: 'all_paths', because odgi stats could emit a list of path names
        for l in mean_links_length:
            if l["length"]["path"] == "all_paths":
                length = l["length"]
        sum_of_path_nodes_distances = data["sum_of_path_node_distances"]
        # we have to find the entry with path: 'all_paths', because odgi stats could emit a list of path names
        for d in sum_of_path_nodes_distances:
            if d["distance"]["path"] == "all_paths":
                distance = d["distance"]
        n = data.get("N", 0)
        num_nodes_self_loops = data["num_nodes_self_loops"]
        return {
            sample_name: {
                "General stats {}".format(sample_name): {
                    "Length": float(data["length"]),
                    "Nodes": float(data["nodes"]),
                    "Edges": float(data["edges"]),
                    "Paths": float(data["paths"]),
                    "Components": float(data["num_weakly_connected_components"]),
                    "A": float(data["A"]),
                    "C": float(data["C"]),
                    "T": float(data["T"]),
                    "G": float(data["G"]),
                    "N": float(n),
                    "Total": float(num_nodes_self_loops["total"]),
                    "Unique": float(num_nodes_self_loops["unique"]),
                },
                "Mean_links_length {}".format(sample_name): {
                    "Path": length["path"],
                    "In_node_space": length["in_node_space"],
                    "In_nucleotide_space": length["in_nucleotide_space"],
                    "Num_links_considered": length["num_links_considered"],
                    "Num_gap_links_not_penalized": length["num_gap_links_not_penalized"],
                },
                "Sum_of_path_nodes_distances {}".format(sample_name): {
                    "Path": distance["path"],
                    "In_node_space": distance["in_node_space"],
                    "In_nucleotide_space": distance["in_nucleotide_space"],
                    "Nodes": distance["nodes"],
                    "Nucleotides": distance["nucleotides"],
                    "Num_penalties": distance["num_penalties"],
                    "Num_penalties_different_orientation": distance["num_penalties_different_orientation"],
                },
            }
        }
