#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from gfastats --stats """

import logging

# Initialise the logger
log = logging.getLogger(__name__)


class StatsReportMixin:
    """Mixin class, loaded by main gfastats MutliqcModule class."""

    def parse_gfastats_stats(self):
        """Find gfastats stats logs and parse summary statistics"""

        self.gfastats_stats = dict()
        for f in self.find_log_files("gfastats/stats"):
            parsed_data = dict()
            for line in f["f"].splitlines()[1:]:
                sections = line.split("\t")
                field = sections[0]
                if not field == "Base composition (ACGT)":
                    value = float(sections[1])
                    parsed_data[field] = value

            if len(parsed_data) > 0:
                if f["s_name"] in self.gfastats_stats:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
                self.add_data_source(f, section="stats")
                self.gfastats_stats[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.gfastats_stats = self.ignore_samples(self.gfastats_stats)

        if len(self.gfastats_stats) > 0:

            # Write parsed report data to a file
            self.write_data_file(self.gfastats_stats, "multiqc_gfastats_stats")

            self.general_stats_addcols(self.gfastats_stats)

        # Return the number of logs that were found
        return len(self.gfastats_stats)
