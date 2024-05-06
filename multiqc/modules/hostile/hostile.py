#!/usr/bin/env python

import json
import logging
import os
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Hostile Module"""

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Hostile",
            anchor="hostile",
            href="https://github.com/bede/hostile",
            info="is a short and long host reads removal tool",
            doi="10.1093/bioinformatics/btad728",
        )

        data_by_sample = dict()
        for f in self.find_log_files("hostile", filehandles=True):
            try:
                parsed_entries = json.load(f["f"])
            except json.JSONDecodeError:
                log.warning(f"Could not parse JSON file {f['f']}")
                continue
            else:
                for entry in parsed_entries:
                    s_name = self.clean_s_name(entry["fastq1_in_name"])

                    if s_name in data_by_sample:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    data_by_sample[s_name] = entry

                    self.add_data_source(f, s_name=s_name)

                    version = entry.get("version")
                    if version:
                        self.add_software_version(version, sample=s_name)

        # Ignore samples
        data_by_sample = self.ignore_samples(data_by_sample)

        # Check if no matching log files were found
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound

        # Number of files found
        log.info(f"Found {len(data_by_sample)} reports")

        # Write hostile data.
        self.write_data_file(data_by_sample, "multiqc_hostile")

        # Plot
        self.hostile_plot(data_by_sample)

    def hostile_plot(self, data_by_sample):
        data = {}

        databases = set()
        for s_name, entry in data_by_sample.items():
            databases.add(os.path.basename(entry["index"]))
            if entry["reads_out"] + entry["reads_removed"] != entry["reads_in"]:
                raise ValueError(
                    f"Reads out + reads removed != reads in for sample {s_name}: "
                    f"{entry['reads_out']} + {entry['reads_removed']} != {entry['reads_in']}"
                )
            data[s_name] = {"Clean reads": entry["reads_out"], "Removed reads": entry["reads_removed"]}

        databases_message = ""
        if len(databases) == 1:
            databases_message = f"Database index: {list(databases)[0]}"
        elif len(databases) > 1:
            log.warning(f"Multiple database indices found in data: {', '.join(list(sorted(databases)))}")
            databases_message = (
                f"<div class='alert alert-warning'>Warning: multiple database indices found in data: "
                f"{', '.join(['<code>' + d + '</code>' for d in sorted(databases)])}. "
                f"Comparison between samples cleaned with different databases might be incorrect</div>"
            )

        cats = ["Clean reads", "Removed reads"]

        pconfig = {
            "title": "Hostile: Reads Filtered",
            "id": "he_reads_plots",
            "ylab": "# Reads",
        }

        self.add_section(
            name="Read Filtering",
            anchor="hostile-reads",
            description=(
                f"The number of reads after filtering (cleaned reads) vs. the number of removed host reads. "
                f"The numbers sum up to the total number of input reads. "
                f"{databases_message}"
            ),
            plot=bargraph.plot(data, cats, pconfig),
        )
