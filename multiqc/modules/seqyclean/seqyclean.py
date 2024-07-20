import logging
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    SeqyClean is a comprehensive preprocessing software application for NGS reads, that removes noise from FastQ
    files to improve de-novo genome assembly and genome mapping.

    The module parses the `*SummaryStatistics.tsv` files that results from a SeqyClean cleaning.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="SeqyClean",
            anchor="seqyclean",
            href="https://github.com/ibest/seqyclean",
            info="Filters adapters, vectors, and contaminants while quality trimming.",
            doi="10.1145/3107411.3107446",
        )

        # Parse logs
        self.seqyclean_data = dict()
        for f in self.find_log_files("seqyclean"):
            rows = f["f"].splitlines()
            headers = rows[0].split("\t")
            cols = rows[1].split("\t")

            self.seqyclean_data[f["s_name"]] = dict()
            for header, col in zip(headers, cols):
                # Add verions info
                if header == "Version":
                    self.add_software_version(col, f["s_name"])

                # Attempt to convert into a float if we can
                try:
                    col = float(col)
                except (ValueError, TypeError):
                    pass
                self.seqyclean_data[f["s_name"]].update({header: col})

            self.add_data_source(f)

        if len(self.seqyclean_data) == 0:
            raise ModuleNoSamplesFound

        self.seqyclean_data = self.ignore_samples(self.seqyclean_data)

        if len(self.seqyclean_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.seqyclean_data)} logs")

        # Adding the bar plots
        self.add_section(
            name="Summary",
            anchor="seqyclean-summary",
            description="This plot shows the number of reads that were kept and discarded.",
            plot=self.seqyclean_summary(),
        )

        self.add_section(
            name="Annotations",
            anchor="seqyclean-annotation",
            description="This plot shows how reads were annotated.",
            plot=self.seqyclean_analysis(),
        )

        self.add_section(
            name="Discarded",
            anchor="seqyclean-discarded",
            description="This plot shows the breakdown of reasons for why reads were discarded.",
            plot=self.seqyclean_discarded(),
        )

        # Write the results to a file
        self.write_data_file(self.seqyclean_data, "multiqc_seqyclean")

        # Adding to the general statistics table
        self.seqyclean_general_stats_table()

    def seqyclean_summary(self):
        config = {"id": "seqyclean-summary-plot", "title": "SeqyClean: Summary", "ylab": "Number of Reads"}
        keys = [
            "PairsKept",  # for paired end
            "PairsDiscarded",
            "PE1DiscardedTotal",
            "PE2DiscardedTotal",
            "SEReadsKept",  # single end
            "SEDiscardedTotal",
            "ReadsKept",  # 454
            "DiscardedTotal",
        ]
        return bargraph.plot(self.seqyclean_data, self._clean_keys(keys), config)

    def seqyclean_analysis(self):
        config = {
            "id": "seqyclean-read-annotation-plot",
            "title": "SeqyClean: Read Annotations",
            "ylab": "Number of Reads",
        }
        keys = [
            "PE1TruSeqAdap_found",  # for paired end
            "PE1ReadsWVector_found",
            "PE1ReadsWContam_found",
            "PE2TruSeqAdap_found",
            "PE2ReadsWVector_found",
            "PE2ReadsWContam_found",
            "SETruSeqAdap_found",  # for single end
            "SEReadsWVector_found",
            "SEReadsWContam_found",
            "left_mid_tags_found",  # 454 data
            "right_mid_tags_found",
            "ReadsWithVector_found",
            "ReadsWithContam_found",
        ]
        return bargraph.plot(self.seqyclean_data, self._clean_keys(keys), config)

    def seqyclean_discarded(self):
        config = {
            "id": "seqyclean-discarded-reads-plot",
            "title": "SeqyClean: Discarded Reads",
            "ylab": "Number of Reads",
        }
        keys = [
            "SEDiscByContam",  # single end
            "SEDiscByLength",
            "PE1DiscByContam",  # paired end
            "PE1DiscByLength",
            "PE2DiscByContam",
            "PE2DiscByLength",
            "DiscByContam",  # 454 data
            "DiscByLength",
        ]
        return bargraph.plot(self.seqyclean_data, self._clean_keys(keys), config)

    def seqyclean_general_stats_table(self):
        headers = {
            "Perc_Kept": {
                "title": "% Kept",
                "description": "The percentage of reads remaining after cleaning",
                "scale": "YlGn",
                "suffix": "%",
                "max": 100,
                "min": 0,
            },
            "PercentageKept": {
                "title": "% Kept",
                "description": "The percentage of reads remaining after cleaning",
                "scale": "YlGn",
                "suffix": "%",
                "max": 100,
                "min": 0,
            },
        }
        self.general_stats_addcols(self.seqyclean_data, headers)

    def _clean_keys(self, keys):
        """Given a list of keys, make them easier to read for plot labels"""
        cats = {}
        for k in keys:
            nice_name = re.sub(r"([a-z])([A-Z])", r"\g<1> \g<2>", k)  # CamelCase > Camel Case
            nice_name = re.sub(r"([PS]E\d?)", r"\g<1> ", nice_name)  # PE1Label > PE1 Label
            nice_name = re.sub(r"W([A-Z])", r"W \g<1>", nice_name)  # WContam > W Contam
            nice_name = nice_name.replace("_", " ")  # tags_found > tags found
            nice_name = nice_name.title()  # Title Case
            nice_name = nice_name.replace("Pe", "PE").replace("Se", "SE")
            nice_name = nice_name.replace("Tru SEq", "TruSeq")
            cats[k] = {"name": nice_name}
        return cats
