import json
import logging
from json import JSONDecodeError

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import heatmap

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This module takes the JSON output of the HOPS postprocessing R script (version >= 0.34) to recreate the
    possible positives heatmap, with the heat intensity representing the number of 'ancient DNA characteristics'
    categories (small edit distance, damage, both edit distance and aDNA damage) that a particular taxon has.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="HOPS",
            anchor="hops",
            href="https://github.com/rhuebler/HOPS/",
            info="Ancient DNA characteristics screening tool of output from the metagenomic aligner MALT.",
            doi="10.1186/s13059-019-1903-0",
        )

        # Find and load any HOPS post-processing JSONs
        self.hops_data = dict()

        for f in self.find_log_files("hops", filehandles=True):
            try:
                self.parseJSON(f)
            except KeyError:
                logging.warning(f"Error loading file {f['fn']}")

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        self.hops_data = self.ignore_samples(self.hops_data)

        if len(self.hops_data) == 0:
            raise ModuleNoSamplesFound

        # Write data to file
        self.write_data_file(self.hops_data, "hops")

        log.info(f"Found {len(self.hops_data)} samples")

        # This type of data isn't 'summarise-able' for general stats, so
        # skipping straight to heatmap. We also won't write data file to the
        # multiqc_data directory because it would be exactly same as input JSON.
        self.hops_heatmap()

    def parseJSON(self, f):
        """Parse the JSON output from HOPS and save the summary statistics"""

        try:
            parsed_json = json.load(f["f"])
        except JSONDecodeError as e:
            log.debug(f"Could not parse HOPS JSON: '{f['fn']}'")
            log.debug(e)
            return None

        # Convert JSON to dict for easier manipulation
        for s in parsed_json:
            s_name = self.clean_s_name(s, f)
            if s_name in self.hops_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name=s_name)
            self.hops_data[s_name] = {}
            for t in parsed_json[s]:
                self.hops_data[s_name][t] = parsed_json[s][t]

    def hops_heatmap(self):
        """Heatmap showing all statuses for every sample"""
        samples = []
        for s in self.hops_data:
            samples.append(s)

        # As all samples always have same taxa, will take from the first sample
        taxa = []
        for t in self.hops_data[samples[0]]:
            taxa.append(t.replace("_", " "))

        # Get values from named list into a list of lists required for heatmap
        data = []
        for s, d in self.hops_data.items():
            row = []
            for val in d.values():
                # Values can be lists of 1 element, so flattening to single value
                if isinstance(val, list):
                    val = val[0]
                row.append(val)
            data.append(row)

        pconfig = {
            "id": "hops-heatmap",
            "title": "HOPS: Potential Candidates",
            "xlab": "Node",
            "ylab": "Sample",
            "square": False,
            "colstops": [
                [1, "#ededed"],
                [2, "#FFFFC5"],
                [3, "#F2B26C"],
                [4, "#AD2A2B"],
            ],
            "tt_decimals": 0,
            "legend": True,
            "display_values": False,
            "xcats_samples": False,
        }

        extra_warning = ""
        if len(self.hops_data) > 20:
            extra_warning = """
            <div class="alert alert-warning">
                Large numbers of samples can result in Y-axis labels
                overlapping. Drag the handle at the bottom of the plot down
                to expand and see all samples names.
            </div>
                """

        self.add_section(
            name="Potential Candidates",
            anchor="hops_heatmap",
            description="""
            Heatmap of candidate taxa for downstream aDNA analysis, with
            intensity representing additive categories of possible 'positive'
            hits.
            """
            + extra_warning,
            helptext="""
            HOPS assigns a category based on how many ancient DNA
            characteristics a given node (i.e. taxon) in a sample has.
            The colours indicate the following:

            * <span style="background-color: #ededed; padding:0.2rem 1rem;">**Grey**</span> - No characteristics detected
            * <span style="background-color: #FFFFC5; padding:0.2rem 1rem;">**Yellow**</span> - Small edit distance from reference
            * <span style="background-color: #F2B26C; padding:0.2rem 1rem;">**Orange**</span> - Typical aDNA damage pattern
            * <span style="background-color: #AD2a2B; padding:0.2rem 1rem;">**Red**</span> - Small edit distance _and_ aDNA damage pattern

            A red category typically indicates a good candidate for further investigation
            in downstream analysis.
            """,
            plot=heatmap.plot(data, xcats=taxa, ycats=samples, pconfig=pconfig),
        )
