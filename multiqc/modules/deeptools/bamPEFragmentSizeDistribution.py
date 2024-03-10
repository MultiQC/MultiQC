""" MultiQC submodule to parse output from deepTools bamPEFragmentSize for read length distribution """

import logging

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class bamPEFragmentSizeDistributionMixin:
    def parse_bamPEFragmentSizeDistribution(self):
        """Find bamPEFragmentSize output. Supports the --outRawFragmentLengths option"""
        self.deeptools_bamPEFragmentSizeDistribution = dict()
        for f in self.find_log_files("deeptools/bamPEFragmentSizeDistribution", filehandles=False):
            parsed_data = self.parseBamPEFDistributionFile(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_bamPEFragmentSizeDistribution:
                    log.warning(f"Replacing duplicate sample {k}.")
                self.deeptools_bamPEFragmentSizeDistribution[k] = v
            if len(parsed_data) > 0:
                self.add_data_source(f, section="bamPEFragmentSizeDistribution")

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        self.deeptools_bamPEFragmentSizeDistribution = self.ignore_samples(self.deeptools_bamPEFragmentSizeDistribution)

        if len(self.deeptools_bamPEFragmentSizeDistribution) > 0:
            # Write data to file
            self.write_data_file(self.deeptools_bamPEFragmentSizeDistribution, "deeptools_frag_size_dist")

            config = {
                "id": "fragment_size_distribution_plot",
                "title": "deeptools: Fragment Size Distribution Plot",
                "ylab": "Occurrence",
                "xlab": "Fragment Size",
                "smooth_points": 50,
                "xmax": 1000,
                "tt_label": "<b>{point.x} bp</b>: {point.y}th occurrence",
            }

            self.add_section(
                name="Fragment size distribution",
                anchor="fragment_size_distribution",
                description="Distribution of paired-end fragment sizes",
                plot=linegraph.plot(self.deeptools_bamPEFragmentSizeDistribution, config),
            )

        return len(self.deeptools_bamPEFragmentSizeDistribution)

    def parseBamPEFDistributionFile(self, f):
        d = dict()
        lastsample = []
        for line in f["f"].splitlines():
            cols = line.rstrip().split("\t")
            if cols[0] == "#bamPEFragmentSize":
                continue
            elif cols[0] == "Size":
                continue
            else:
                s_name = self.clean_s_name(cols[2].rstrip().split("/")[-1], f)
                if s_name != lastsample:
                    d[s_name] = dict()
                    lastsample = s_name
                d[s_name].update({self._int(cols[0]): self._int(cols[1])})

        return d
