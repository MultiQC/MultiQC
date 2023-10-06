""" MultiQC module to parse output from QUAST """


import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="QUAST",
            anchor="quast",
            href="http://quast.bioinf.spbau.ru/",
            info="is a quality assessment tool for genome assemblies, written by the Center for Algorithmic Biotechnology.",
            doi="10.1093/bioinformatics/btt086",
        )

        # Get modifiers from config file
        qconfig = getattr(config, "quast_config", {})

        self.contig_length_multiplier = qconfig.get("contig_length_multiplier", 0.001)
        self.contig_length_suffix = qconfig.get("contig_length_suffix", "Kbp")

        self.total_length_multiplier = qconfig.get("total_length_multiplier", 0.000001)
        self.total_length_suffix = qconfig.get("total_length_suffix", "Mbp")

        self.total_number_contigs_multiplier = qconfig.get("total_number_contigs_multiplier", 0.001)
        self.total_number_contigs_suffix = qconfig.get("total_number_contigs_suffix", "K")

        # Find and load any QUAST reports
        self.quast_data = dict()
        for f in self.find_log_files("quast"):
            self.parse_quast_log(f)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        # Filter to strip out ignored sample names
        self.quast_data = self.ignore_samples(self.quast_data)

        if len(self.quast_data) == 0:
            raise ModuleNoSamplesFound

        log.info("Found {} reports".format(len(self.quast_data)))

        # Write parsed report data to a file
        self.write_data_file(self.quast_data, "multiqc_quast")

        # Basic Stats Table
        self.quast_general_stats_table()

        # Quast Stats Table
        self.add_section(name="Assembly Statistics", anchor="quast-stats", plot=self.quast_table())
        # Number of contigs plot
        self.add_section(
            name="Number of Contigs",
            anchor="quast-contigs",
            description="""This plot shows the number of contigs found for each assembly, broken
                    down by length.""",
            plot=self.quast_contigs_barplot(),
        )
        # Number of genes plot
        ng_pdata = self.quast_predicted_genes_barplot()
        if ng_pdata:
            self.add_section(
                name="Number of Predicted Genes",
                anchor="quast-genes",
                description="""This plot shows the number of predicted genes found for each
                          assembly, broken down by length.""",
                plot=ng_pdata,
            )
        # Number of partial genes plot
        ng_pdata = self.quast_predicted_genes_barplot(partial=True)
        if ng_pdata:
            self.add_section(
                name="Number of Partially Predicted Genes",
                anchor="quast-partial-genes",
                description="""This plot shows the number of partially predicted genes found for each
                          assembly, broken down by length.""",
                plot=ng_pdata,
            )

    def parse_quast_log(self, f):
        lines = f["f"].splitlines()

        # Pull out the sample names from the first row
        s_names = lines[0].split("\t")
        # Prepend directory name(s) to sample names as configured
        s_names = [self.clean_s_name(s_name, f) for s_name in s_names]
        for s_name in s_names[1:]:
            if s_name in self.quast_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.quast_data[s_name] = dict()

        # Parse remaining stats for each sample
        for l in lines[1:]:
            s = l.split("\t")
            k = s[0]
            for i, v in enumerate(s[1:]):
                s_name = s_names[i + 1]
                partials = re.search(r"(\d+) \+ (\d+) part", v)
                if partials:
                    whole = partials.group(1)
                    partial = partials.group(2)
                    try:
                        self.quast_data[s_name][k] = float(whole)
                        self.quast_data[s_name]["{}_partial".format(k)] = float(partial)
                    except ValueError:
                        self.quast_data[s_name][k] = whole
                        self.quast_data[s_name]["{}_partial".format(k)] = partial
                else:
                    try:
                        self.quast_data[s_name][k] = float(v)
                    except ValueError:
                        self.quast_data[s_name][k] = v

    def quast_general_stats_table(self):
        """Take the parsed stats from the QUAST report and add some to the
        General Statistics table at the top of the report"""

        headers = OrderedDict()
        headers["N50"] = {
            "title": "N50 ({})".format(self.contig_length_suffix),
            "description": "N50 is the contig length such that using longer or equal length contigs produces half (50%) of the bases of the assembly (kilo base pairs)",
            "min": 0,
            "suffix": self.contig_length_suffix,
            "scale": "RdYlGn",
            "modify": lambda x: x * self.contig_length_multiplier,
        }
        headers["Total length"] = {
            "title": "Assembly Length ({})".format(self.total_length_suffix),
            "description": "The total number of bases in the assembly ({}).".format(self.total_length_suffix),
            "min": 0,
            "suffix": self.total_length_suffix,
            "scale": "YlGn",
            "modify": lambda x: x * self.total_length_multiplier,
        }
        self.general_stats_addcols(self.quast_data, headers)

    def quast_table(self):
        """Write some more statistics about the assemblies in a table."""
        headers = OrderedDict()

        headers["N50"] = {
            "title": "N50 ({})".format(self.contig_length_suffix),
            "description": "N50 is the contig length such that using longer or equal length contigs produces half (50%) of the bases of the assembly.",
            "min": 0,
            "suffix": self.contig_length_suffix,
            "scale": "RdYlGn",
            "modify": lambda x: x * self.contig_length_multiplier,
        }

        headers["N75"] = {
            "title": "N75 ({})".format(self.contig_length_suffix),
            "description": "N75 is the contig length such that using longer or equal length contigs produces 75% of the bases of the assembly",
            "min": 0,
            "suffix": self.contig_length_suffix,
            "scale": "RdYlGn",
            "modify": lambda x: x * self.contig_length_multiplier,
        }

        headers["L50"] = {
            "title": "L50 ({})".format(self.total_number_contigs_suffix) if self.total_number_contigs_suffix else "L50",
            "description": "L50 is the number of contigs larger than N50, i.e. the minimum number of contigs comprising 50% of the total assembly length.",
            "min": 0,
            "suffix": self.total_number_contigs_suffix,
            "scale": "RdYlGn-rev",
            "modify": lambda x: x * self.total_number_contigs_multiplier,
        }

        headers["L75"] = {
            "title": "L75 ({})".format(self.total_number_contigs_suffix) if self.total_number_contigs_suffix else "L75",
            "description": "L75 is the number of contigs larger than N75, i.e. the minimum number of contigs comprising 75% of the total assembly length.",
            "min": 0,
            "suffix": self.total_number_contigs_suffix,
            "scale": "RdYlGn-rev",
            "modify": lambda x: x * self.total_number_contigs_multiplier,
        }
        headers["Largest contig"] = {
            "title": "Largest contig ({})".format(self.contig_length_suffix),
            "description": "The size of the largest contig of the assembly",
            "min": 0,
            "suffix": self.contig_length_suffix,
            "scale": "YlGn",
            "modify": lambda x: x * self.contig_length_multiplier,
        }

        headers["Total length"] = {
            "title": "Length ({})".format(self.total_length_suffix),
            "description": "The total number of bases in the assembly.",
            "min": 0,
            "suffix": self.total_length_suffix,
            "scale": "YlGn",
            "modify": lambda x: x * self.total_length_multiplier,
        }

        headers["# misassemblies"] = {
            "title": "Misassemblies",
            "description": "The number of positions in the assembled contigs where the left flanking sequence aligns over 1 kbp away from the right flanking sequence on the reference (relocation) or they overlap on more than 1 kbp (relocation) or flanking sequences align on different strands (inversion) or different chromosomes (translocation).",
            "scale": "RdYlGn-rev",
            "format": "{,:.0f}",
        }
        headers["# mismatches per 100 kbp"] = {
            "title": "Mismatches/100kbp",
            "description": "The number of mismatches per 100 kbp",
            "scale": "YlOrRd",
            "format": "{:,.2f}",
        }
        headers["# indels per 100 kbp"] = {
            "title": "Indels/100kbp",
            "description": "The number of indels per 100 kbp",
            "scale": "YlOrRd",
            "format": "{:,.2f}",
        }
        headers["# genes"] = {
            "title": "Genes",
            "description": "# Genes",
            "scale": "YlGnBu",
            "format": "{:,.0f}",
            "shared_key": "gene_count",
        }
        headers["# genes_partial"] = {
            "title": "Genes (Partial)",
            "description": "# Genes (Partial)",
            "scale": "YlGnBu",
            "format": "{:,.0f}",
            "shared_key": "gene_count",
        }
        headers["# predicted genes (unique)"] = {
            "title": "Genes",
            "description": "# Predicted Genes (Unique)",
            "scale": "YlGnBu",
            "format": "{:,.0f}",
            "shared_key": "gene_count",
        }
        headers["Genome fraction (%)"] = {
            "title": "Genome Fraction",
            "description": "The total number of aligned bases in the reference, divided by the genome size.",
            "max": 100,
            "suffix": "%",
            "scale": "YlGn",
        }
        config = {
            "id": "quast_table",
            "min": 0,
        }
        return table.plot(self.quast_data, headers, config)

    def quast_contigs_barplot(self):
        """Make a bar plot showing the number and length of contigs for each assembly"""

        # Prep the data
        data = dict()
        categories = []
        for s_name, d in self.quast_data.items():
            nums_by_t = dict()
            for k, v in d.items():
                m = re.match(r"# contigs \(>= (\d+) bp\)", k)
                if m and v != "-":
                    nums_by_t[int(m.groups()[0])] = int(v)

            tresholds = sorted(nums_by_t.keys(), reverse=True)
            p = dict()
            cats = []
            for i, t in enumerate(tresholds):
                if i == 0:
                    c = ">= " + str(t) + " bp"
                    cats.append(c)
                    p[c] = nums_by_t[t]
                else:
                    c = str(t) + "-" + str(tresholds[i - 1]) + " bp"
                    cats.append(c)
                    p[c] = nums_by_t[t] - nums_by_t[tresholds[i - 1]]
            if not categories:
                categories = cats
            data[s_name] = p

        pconfig = {
            "id": "quast_num_contigs",
            "title": "QUAST: Number of Contigs",
            "ylab": "# Contigs",
            "yDecimals": False,
        }

        return bargraph.plot(data, categories, pconfig)

    def quast_predicted_genes_barplot(self, partial=False):
        """
        Make a bar plot showing the number and length of predicted genes
        for each assembly
        """

        # Prep the data
        # extract the ranges given to quast with "--gene-thresholds"
        # keys look like:
        #   `# predicted genes (>= 300 bp)`
        #   `# predicted genes (>= 300 bp)_partial`
        pattern = re.compile(r"# predicted genes \(>= (\d+) bp\)" + ("_partial" if partial else ""))

        data = {}
        all_categories = []
        data_key = "# predicted genes (>= {} bp)" + ("_partial" if partial else "")
        for s_name, d in self.quast_data.items():
            thresholds = []
            for k in d.keys():
                m = re.match(pattern, k)
                if m:
                    thresholds.append(int(m.groups()[0]))
            thresholds = sorted(list(set(thresholds)))
            if len(thresholds) < 2:
                continue

            highest_threshold = thresholds[-1]
            highest_cat = (highest_threshold, ">= {} bp".format(highest_threshold))  # tuple (key-for-sorting, label)
            all_categories.append(highest_cat)
            plot_data = {highest_cat[1]: d[data_key.format(highest_threshold)]}

            # converting >=T1, >=T2,.. into 0-T1, T1-T2,..
            for low, high in zip(thresholds[:-1], thresholds[1:]):
                cat = (low, "{}-{} bp".format(low, high))
                all_categories.append(cat)
                plot_data[cat[1]] = d[data_key.format(low)] - d[data_key.format(high)]

            data[s_name] = plot_data

        all_categories = [label for k, label in sorted(list(set(all_categories)))]

        if len(all_categories) > 0:
            return bargraph.plot(
                data,
                all_categories,
                {
                    "id": "quast_" + ("partially_" if partial else "") + "predicted_genes",
                    "title": "QUAST: Number of " + ("partially " if partial else "") + "predicted genes",
                    "ylab": "Number of " + ("partially " if partial else "") + "predicted genes",
                },
            )
        else:
            return None
