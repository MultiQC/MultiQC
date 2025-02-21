"""MultiQC submodule to parse output from deepTools bamPEFragmentSize for summary table"""

import logging

from multiqc.plots import linegraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class bamPEFragmentSizeTableMixin:
    def parse_bamPEFragmentSize(self):
        """Find bamPEFragmentSize output. Supports the --table option"""
        self.deeptools_bamPEFragmentSize = dict()
        for f in self.find_log_files("deeptools/bamPEFragmentSizeTable"):
            parsed_data = self.parseBamPEFile(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_bamPEFragmentSize:
                    log.warning(f"Replacing duplicate sample {k}.")
                self.deeptools_bamPEFragmentSize[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section="bamPEFragmentSize")

        self.deeptools_bamPEFragmentSize = self.ignore_samples(self.deeptools_bamPEFragmentSize)

        if len(self.deeptools_bamPEFragmentSize) == 0:
            return 0

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write data to file
        self.write_data_file(self.deeptools_bamPEFragmentSize, "deeptools_frag_size_table")

        headersSE = {
            "Reads Sampled": {
                "title": "# Sampled",
                "description": "Number of reads sampled",
                "format": "{:,.0f}",
            },
            "Read Len. Min.": {
                "title": "Min",
                "description": "Minimum read length",
                "format": "{:,.0f}",
                "shared_key": "read_length",
            },
            "Read Len. 1st. Qu.": {
                "title": "1st Quartile",
                "description": "1st quartile read length",
                "format": "{:,.0f}",
                "shared_key": "read_length",
            },
            "Read Len. Mean": {
                "title": "Mean",
                "description": "Mean read length",
                "shared_key": "read_length",
            },
            "Read Len. Median": {
                "title": "Median",
                "description": "Median read length",
                "format": "{:,.0f}",
                "shared_key": "read_length",
            },
            "Read Len. 3rd Qu.": {
                "title": "3rd Quartile",
                "description": "3rd quartile read length",
                "format": "{:,.0f}",
                "shared_key": "read_length",
            },
            "Read Len. Max": {
                "title": "Max",
                "description": "Maximum read length",
                "format": "{:,.0f}",
                "shared_key": "read_length",
            },
            "Read Len. Std.": {
                "title": "Std. Dev.",
                "description": "read length standard deviation",
                "shared_key": "read_length",
            },
            "Read Med. Abs. Dev.": {
                "title": "MAD",
                "description": "read length median absolute deviation",
                "shared_key": "read_length",
            },
        }
        config = {
            "namespace": "deepTools bamPEFragmentSize",
            "id": "deeptools_readlengths_table",
            "title": "deepTools: Read length metrics",
        }
        self.add_section(
            name="Read length metrics",
            anchor="deeptools_readlengths",
            plot=table.plot(self.deeptools_bamPEFragmentSize, headersSE, config),
        )

        headersPE = {
            "Frag. Sampled": {
                "title": "# Sampled",
                "description": "Number of fragments sampled",
                "format": "{:,.0f}",
            },
            "Frag. Len. Min.": {
                "title": "Min",
                "description": "Minimum fragment length",
                "format": "{:,.0f}",
                "shared_key": "frag_length",
            },
            "Frag. Len. 1st. Qu.": {
                "title": "1st Quartile",
                "description": "1st quartile fragment length",
                "format": "{:,.0f}",
                "shared_key": "frag_length",
            },
            "Frag. Len. Mean": {
                "title": "Mean",
                "description": "Mean fragment length",
                "format": "{:,.0f}",
                "shared_key": "frag_length",
            },
            "Frag. Len. Median": {
                "title": "Median",
                "description": "Median fragment length",
                "format": "{:,.0f}",
                "shared_key": "frag_length",
            },
            "Frag. Len. 3rd Qu.": {
                "title": "3rd Quartile",
                "description": "3rd quartile fragment length",
                "format": "{:,.0f}",
                "shared_key": "frag_length",
            },
            "Frag. Len. Max": {
                "title": "Max",
                "description": "Maximum fragment length",
                "format": "{:,.0f}",
                "shared_key": "frag_length",
            },
            "Frag. Len. Std.": {
                "title": "Std. Dev.",
                "description": "Fragment length standard deviation",
                "shared_key": "frag_length",
            },
            "Frag. Med. Abs. Dev.": {
                "title": "MAD",
                "description": "Fragment length median absolute deviation",
                "shared_key": "frag_length",
            },
        }

        # Are there any PE datasets?
        PE = False
        for k, v in self.deeptools_bamPEFragmentSize.items():
            if "Frag. Len. Min." in v:
                PE = True
                break
        if PE:
            config = {
                "namespace": "deepTools bamPEFragmentSize",
                "id": "deeptools_fragmentlengths_table",
                "title": "deepTools: Fragment length metrics",
            }
            self.add_section(
                name="Fragment length metrics",
                anchor="deeptools_fragmentlengths",
                plot=table.plot(self.deeptools_bamPEFragmentSize, headersPE, config),
            )

        # Read length plot
        config = {
            "data_labels": [
                {
                    "name": "Read length distribution",
                    "title": "Read length distribution",
                    "ylab": "Read length (bases)",
                },
                {
                    "name": "Fragment length distribution",
                    "title": "Fragment length distribution",
                    "ylab": "Fragment length (bases)",
                },
            ],
            "id": "deeptools_readlengthsPlot",
            "title": "deepTools: Read/Fragment length distribution",
            "ylab": "Read length (bases)",
            "xlab": "Percentile",
        }
        SE = dict()
        PE = dict()
        for k, v in self.deeptools_bamPEFragmentSize.items():
            SE[k] = {
                0: v["Read Len. Min."],
                10: v["Read Len. 10%"],
                20: v["Read Len. 20%"],
                25: v["Read Len. 1st. Qu."],
                30: v["Read Len. 30%"],
                40: v["Read Len. 40%"],
                50: v["Read Len. Median"],
                60: v["Read Len. 60%"],
                70: v["Read Len. 70%"],
                75: v["Read Len. 3rd Qu."],
                80: v["Read Len. 80%"],
                90: v["Read Len. 90%"],
                99: v["Read Len. 99%"],
                100: v["Read Len. Max"],
            }
            if "Frag. Len. Min." not in v:
                continue
            PE[k] = {
                0: v["Frag. Len. Min."],
                10: v["Frag. Len. 10%"],
                20: v["Frag. Len. 20%"],
                25: v["Frag. Len. 1st. Qu."],
                30: v["Frag. Len. 30%"],
                40: v["Frag. Len. 40%"],
                50: v["Frag. Len. Median"],
                60: v["Frag. Len. 60%"],
                70: v["Frag. Len. 70%"],
                75: v["Frag. Len. 3rd Qu."],
                80: v["Frag. Len. 80%"],
                90: v["Frag. Len. 90%"],
                99: v["Frag. Len. 99%"],
                100: v["Frag. Len. Max"],
            }
        self.add_section(
            name="Read/fragment length distribution",
            anchor="deeptools_fragmentlengths_dist",
            plot=linegraph.plot([SE, PE], config),
        )

        return len(self.deeptools_bamPEFragmentSize)

    def parseBamPEFile(self, f):
        d = {}
        headers = None
        for line in f["f"].splitlines():
            cols = line.rstrip().split("\t")
            if headers is None:
                headers = cols
            else:
                s_name = None
                for idx, h in enumerate(headers):
                    if idx == 0:
                        s_name = self.clean_s_name(cols[0], f)
                        if s_name in d:
                            log.debug(f"Replacing duplicate sample {s_name}.")
                        d[s_name] = dict()
                    else:
                        if idx < 19 and cols[1] == "0":
                            # Don't store fragment metrics for SE datasets, they're just 0.
                            continue
                        try:
                            # Most values are ac
                            d[s_name][h] = self._int(cols[idx])
                        except ValueError:
                            d[s_name][h] = float(cols[idx])

        return d
