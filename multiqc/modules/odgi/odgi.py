#!/usr/bin/env python

""" MultiQC module to parse output from odgi stats """

from __future__ import print_function
import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

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
            info="is an optimized dynamic graph/genome implementation, for efficient analysis and manipulation of pangenome graphs structured in the variation graph model.",
            # Can't find a DOI // doi=
        )

        # Parse odgi stats data
        self.odgi_stats_map = dict()
        for f in self.find_log_files("odgi", filehandles=True):
            self.parse_odgi_stats_report(f)

        self.odgi_stats_map = self.ignore_samples(self.odgi_stats_map)

        # No samples found
        if len(self.odgi_stats_map) == 0:
            raise UserWarning

        log.info("Found {} report{}".format(len(self.odgi_stats_map), "s" if len(self.odgi_stats_map) > 1 else ""))

        # Write parsed report data to a file
        self.write_data_file(self.odgi_stats_map, "multiqc_odgi_stats")

        # Plot detailed odgi stats in an extra section
        self.odgi_stats_table()

        # Plot the odgi stats metrics as a lineplot
        self.plot_sum_of_path_nodes_distances()
        self.plot_mean_links_length()

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
        # Load the YAML file
        try:
            data = yaml.load(f["f"], Loader=yaml.SafeLoader)
        except Exception as e:
            log.warning("Could not parse YAML for '{}': \n  {}".format(f, e))
            return

        # Flatten the data structure
        data_flat = self.compress_stats_data(data)

        # Calculate additional metrics
        try:
            data_flat["pct_gc"] = ((data_flat["G"] + data_flat["C"]) / data_flat["length"]) * 100.0
            data_flat["pct_n"] = (data_flat["N"] / data_flat["length"]) * 100.0
        except ZeroDivisionError:
            data_flat["pct_gc"] = 0.0
            data_flat["pct_n"] = 0.0

        # Clean up the sample name
        s_name = self.extract_sample_name(f["fn"])
        s_name = self.clean_s_name(s_name, f)

        # Add to the data sources file
        self.add_data_source(f, s_name)

        # Add the data to the main odgi_stats_map dict
        if s_name in self.odgi_stats_map:
            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        self.odgi_stats_map[s_name] = data_flat

    def odgi_stats_table(self):
        """
        Detailed odgi stats in this extra table.
        """
        headers = OrderedDict()
        headers["length"] = {
            "title": "Length",
            "description": "Graph length in nucleotides.",
            "scale": "BuPu",
            "format": "{:,.0f}",
        }
        headers["nodes"] = {
            "title": "Nodes",
            "description": "Number of nodes in the graph.",
            "scale": "OrRd",
            "format": "{:,.0f}",
        }
        headers["edges"] = {
            "title": "Edges",
            "description": "Number of edges in the graph.",
            "scale": "PuBu",
            "format": "{:,.0f}",
        }
        headers["paths"] = {
            "title": "Paths",
            "description": "Number of paths in the graph.",
            "scale": "Greens",
            "format": "{:,.0f}",
        }
        headers["components"] = {
            "title": "Components",
            "description": "Number of weakly connected components in the graph.",
            "scale": "Oranges",
            "format": "{:,.0f}",
        }
        headers["A"] = {
            "title": "A",
            "description": "Number of adenine bases in the graph.",
            "scale": "Spectral",
            "format": "{:,.0f}",
            "shared_key": "nucleotides",
        }
        headers["C"] = {
            "title": "C",
            "description": "Number of cytosine bases in the graph.",
            "scale": "Greys",
            "format": "{:,.0f}",
            "shared_key": "nucleotides",
        }
        headers["T"] = {
            "title": "T",
            "description": "Number of thymine bases in the graph.",
            "scale": "Blues",
            "format": "{:,.0f}",
            "shared_key": "nucleotides",
        }
        headers["G"] = {
            "title": "G",
            "description": "Number of guanine bases in the graph.",
            "scale": "RdPu",
            "format": "{:,.0f}",
            "shared_key": "nucleotides",
        }
        headers["N"] = {
            "title": "N",
            "description": "Number of `N` basis in the graph.",
            "scale": "Set3",
            "format": "{:,.0f}",
            "shared_key": "nucleotides",
        }
        headers["total"] = {
            "title": "Self Loops Nodes",
            "description": "Total number of nodes having self loops in the graph.",
            "scale": "Set1",
            "hidden": True,
            "format": "{:,.0f}",
        }
        headers["unique"] = {
            "title": "Unique Self Loops Nodes",
            "description": "Number of unique nodes having self loops in the graph.",
            "scale": "Set2",
            "hidden": True,
            "format": "{:,.0f}",
        }
        headers["pct_gc"] = {
            "title": "% GC",
            "description": "Percent of G/C bases in the graph.",
            "scale": "Spectral",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "hidden": True,
        }
        headers["pct_n"] = {
            "title": "% N",
            "description": "Percent of N bases in the graph.",
            "scale": "Reds",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "hidden": True,
        }
        # Some of the headers are quite general and can clash with other modules.
        # Prepend odgi_ to keep them unique
        prefix_headers = OrderedDict()
        prefix_data = {}
        for h, v in headers.items():
            prefix_headers[f"odgi_{h}"] = v
        for s_name, d in self.odgi_stats_map.items():
            prefix_data[s_name] = {}
            for h in headers:
                prefix_data[s_name][f"odgi_{h}"] = d[h]
        tconfig = {
            "id": "odgi_table",
            "namespace": "ODGI",
            "table_title": "ODGI Stats",
        }
        self.add_section(
            name="Detailed ODGI stats table.",
            anchor="extended_odgi_stats",
            plot=table.plot(prefix_data, prefix_headers, tconfig),
        )

    def plot_sum_of_path_nodes_distances(self):
        """
        Plot odgi path node distances bar plot.
        """

        cats = [OrderedDict(), OrderedDict()]
        cats[0]["sum_of_path_nodes_distances_in_node_space"] = {"name": "Sum of distances: node space"}
        cats[1]["sum_of_path_nodes_distances_in_nucleotide_space"] = {"name": "Sum of distances: nucleotide space"}
        pconfig = {
            "id": "odgi_sum_of_path_nodes_distances_plot",
            "title": "ODGI: Sum of path node distances",
            "ylab": "Distance",
            "yDecimals": False,
            "cpswitch": False,
            "logswitch": True,
            "logswitch_active": True,
            "data_labels": [
                {"name": "Node space", "ylab": "Distance"},
                {"name": "Nucleotide space", "ylab": "Distance"},
            ],
        }
        self.add_section(
            name="Sum of path node distances",
            anchor="odgi_sum_of_path_nodes_distances",
            description="""
                For each path we iterate from node to node and count the node / nucleotide distance of nodes on the pangenome
                level normalized by the path length. If a node is reversed, we count the node distance twice.
            """,
            helptext="""
                This value allows you to evaluate the sorting goodness - how linear the graph is.
            """,
            plot=bargraph.plot([self.odgi_stats_map, self.odgi_stats_map], cats, pconfig),
        )

    def plot_mean_links_length(self):
        """
        Plot odgi path node distances bar plot.
        """
        cats = [OrderedDict(), OrderedDict()]
        cats[0]["mean_links_length_in_node_space"] = {"name": "Mean links length: node space"}
        cats[1]["mean_links_length_in_nucleotide_space"] = {"name": "Mean links length: nucleotide space"}
        pconfig = {
            "id": "odgi_mean_links_length_plot",
            "title": "ODGI: Mean links length",
            "ylab": "Distance",
            "yDecimals": False,
            "cpswitch": False,
            "logswitch": True,
            "logswitch_active": True,
            "data_labels": [
                {"name": "Node space", "ylab": "Distance"},
                {"name": "Nucleotide space", "ylab": "Distance"},
            ],
        }

        self.add_section(
            name="Mean links length",
            anchor="odgi_mean_links_length",
            description="""
                For each path we iterate from node to node and count the node / nucleotide distance `mean_links_length`
                of nodes _within the same path only._ We then normalized by the path length.
            """,
            helptext="""
                This value allows you to evaluate the sorting goodness - how linear the graph is.
            """,
            plot=bargraph.plot([self.odgi_stats_map, self.odgi_stats_map], cats, pconfig),
        )

    def extract_sample_name(self, file_name):
        """
        Extracts the sample name from a given file name.
        Expects and returns one of seqwish, smooth or cons@*.
        If none of the above keywords are present, the full filename except for the ending is returned.
        """
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

    def compress_stats_data(self, data) -> dict:
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
        return {
            "length": float(data["length"]),
            "nodes": float(data["nodes"]),
            "edges": float(data["edges"]),
            "paths": float(data["paths"]),
            "components": float(data["num_weakly_connected_components"]),
            "A": float(data["A"]),
            "C": float(data["C"]),
            "T": float(data["T"]),
            "G": float(data["G"]),
            "N": float(data.get("N", 0)),
            "total": float(data["num_nodes_self_loops"]["total"]),
            "unique": float(data["num_nodes_self_loops"]["unique"]),
            "mean_links_length_path": length["path"],
            "mean_links_length_in_node_space": length["in_node_space"],
            "mean_links_length_in_nucleotide_space": length["in_nucleotide_space"],
            "mean_links_length_num_links_considered": length["num_links_considered"],
            "mean_links_length_num_gap_links_not_penalized": length["num_gap_links_not_penalized"],
            "sum_of_path_nodes_distances_path": distance["path"],
            "sum_of_path_nodes_distances_in_node_space": distance["in_node_space"],
            "sum_of_path_nodes_distances_in_nucleotide_space": distance["in_nucleotide_space"],
            "sum_of_path_nodes_distances_nodes": distance["nodes"],
            "sum_of_path_nodes_distances_nucleotides": distance["nucleotides"],
            "sum_of_path_nodes_distances_num_penalties": distance["num_penalties"],
            "sum_of_path_nodes_distances_num_penalties_different_orientation": distance[
                "num_penalties_different_orientation"
            ],
        }
