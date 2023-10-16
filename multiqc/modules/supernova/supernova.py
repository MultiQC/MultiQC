""" MultiQC module to parse output from Supernova"""


import json
import logging
import re
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Supernova",
            anchor="supernova",
            href="https://www.10xgenomics.com/",
            info="is a de novo genome assembler 10X Genomics linked-reads.",
            doi="10.1101/gr.214874.116",
        )

        # Headers for the supernova Table
        self.headers = OrderedDict()
        self.headers["Asm size"] = {
            "description": "assembly size (in megabases) ;only scaffolds >= 10 kb",
            "modify": lambda x: x / 1000000.0,
            "suffix": "Mb",
            "scale": "YlGn",
        }
        self.headers["% missing 10Kb"] = {
            "rid": "pct_missing_10Kb",
            "description": "% of base assembly missing from scaffolds >= 10 kb",
            "suffix": "%",
            "scale": "YlGn",
        }
        self.headers["# Long scaffs"] = {
            "rid": "num_long_scaffs",
            "description": "number of scaffolds >= 10 kb",
            "scale": "YlGn",
            "format": "{:,.0f}",
            "hidden": True,
        }
        self.headers["Scaff N50"] = {
            "description": "N50 scaffold size (in kilobases)",
            "modify": lambda x: x / 1000.0,
            "suffix": "Kb",
            "scale": "RdYlGn",
        }
        self.headers["Phase N50"] = {
            "description": "N50 phase block size (in kilobases)",
            "modify": lambda x: x / 1000.0,
            "suffix": "Kb",
            "scale": "RdYlGn",
            "hidden": True,
        }
        self.headers["Contig N50"] = {
            "description": "N50 contig size (in kilobases)",
            "modify": lambda x: x / 1000.0,
            "suffix": "Kb",
            "scale": "RdYlGn",
            "hidden": True,
        }
        self.headers["Edge N50"] = {
            "description": "N50 edge size (in kilobases)",
            "modify": lambda x: x / 1000.0,
            "suffix": "Kb",
            "scale": "RdYlGn",
            "hidden": True,
        }
        self.headers["Mol size"] = {
            "description": "weighted mean molecule size (in kilobases); ideal 50-100",
            "modify": lambda x: x / 1000.0,
            "suffix": "Kb",
            "scale": "BuGn",
        }
        self.headers["Read len"] = {
            "description": "mean read length (in bases) after trimming; ideal 140",
            "suffix": "b",
            "scale": "PuBu",
            "format": "{:,.0f}",
            "hidden": True,
        }
        self.headers["# Reads"] = {
            "rid": "num_reads",
            "description": "number of reads (in millions); ideal 800M-1200M for human",
            "modify": lambda x: x / 1000000.0,
            "suffix": "M",
            "scale": "PuBu",
        }
        self.headers["Raw coverage"] = {
            "description": "raw coverage; ideal ~56",
            "suffix": "x",
            "scale": "PuBu",
            "hidden": True,
        }
        self.headers["Coverage"] = {
            "description": "effective read coverage; ideal ~42 for nominal 56x cov",
            "suffix": "x",
            "scale": "PuBu",
        }
        self.headers["% Dup"] = {
            "rid": "pct_Dup",
            "description": "fraction of reads that are duplicates",
            "suffix": "%",
            "scale": "OrRd",
        }
        self.headers["% R2 Q30"] = {
            "rid": "pct_R2_Q30",
            "description": "fraction of Q30 bases in read 2; ideal 75-85%",
            "suffix": "%",
            "scale": "OrRd",
        }
        self.headers["Insert size"] = {
            "description": "median insert size (in bases); ideal 0.35-0.40 Kb",
            "suffix": "b",
            "scale": "OrRd",
            "format": "{:,.0f}",
            "hidden": True,
        }
        self.headers["% proper"] = {
            "rid": "pct_proper",
            "description": "fraction of proper read pairs; ideal >= 75%",
            "suffix": "%",
            "scale": "OrRd",
            "hidden": True,
        }
        self.headers["BC usage"] = {
            "description": "fraction of barcodes used; between 0 and 1",
            "scale": "OrRd",
            "hidden": True,
        }
        self.headers["Est size"] = {
            "description": "estimated genome size",
            "modify": lambda x: x / 1000000.0,
            "suffix": "Mb",
            "scale": "YlGn",
            "hidden": True,
        }
        self.headers["% repeats"] = {
            "rid": "pct_repeats",
            "description": "Estimated repetitive fraction (of genome)",
            "scale": "YlGn",
            "suffix": "%",
            "hidden": True,
        }
        self.headers["% AT"] = {
            "rid": "pct_AT",
            "description": "high AT index (of genome)",
            "scale": "YlGn",
            "suffix": "%",
            "hidden": True,
        }
        self.headers["Het dist"] = {
            "description": "mean distance between heterozygous SNPs (in kilobases)",
            "modify": lambda x: x / 1000.0,
            "suffix": "Kb",
            "scale": "YlGn",
            "format": "{:,.0f}",
            "hidden": True,
        }
        self.headers["p10"] = {
            "description": "molecule count extending 10 kb on both sides",
            "scale": "BuGn",
            "hidden": True,
        }
        self.headers["% missing BC"] = {
            "rid": "pct_missing_BC",
            "description": "fraction of reads that are not barcoded",
            "suffix": "%",
            "scale": "BuGn",
        }
        self.headers["Barcode N50"] = {
            "description": "N50 reads per barcode (in bases)",
            "suffix": "b",
            "scale": "BuGn",
            "format": "{:,.0f}",
        }
        self.headers["% Phased"] = {
            "rid": "pct_Phased",
            "description": "nonduplicate and phased reads; ideal 45-50%",
            "suffix": "%",
            "scale": "BuGn",
            "hidden": True,
        }

        reports = OrderedDict()
        summaries = OrderedDict()
        molecules = OrderedDict()
        kmers = OrderedDict()
        root_summary = {}

        ### Parse the input log files
        # report.txt files
        for f in self.find_log_files("supernova/report"):
            log.debug("Found report in: {}".format(f["root"]))
            sid, data = self.parse_report(f["f"])
            s_name = self.clean_s_name(sid, f)
            if s_name in reports.keys():
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            reports[s_name] = data
            self.add_data_source(f, s_name=s_name, section="supernova-table")

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

        # summary.json files
        for f in self.find_log_files("supernova/summary"):
            log.debug("Found summary.json in: {}".format(f["root"]))
            try:
                sid, data = self.parse_summary(f["f"])
            except ValueError:
                log.debug("Error parsing JSON file in {}".format(f["root"]))
                continue
            except RuntimeError:
                log.debug("Could not find sample_id in JSON file in {}".format(f["root"]))
                continue

            s_name = self.clean_s_name(sid, f)
            if s_name in summaries.keys():
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            summaries[s_name] = data
            self.add_data_source(f, s_name=s_name, section="supernova-table")
            # The plot json files do not contain sample IDs, sadly. So we need to store it somewhere.
            root_summary[f["root"]] = sid

        # histogram_molecules.json files
        for f in self.find_log_files("supernova/molecules"):
            log.debug("Found histogram_molecules.json in: {}".format(f["root"]))
            try:
                if f["root"] in root_summary.keys():
                    data = self.parse_histogram(f["f"])
                    sid = root_summary[f["root"]]
                    s_name = self.clean_s_name(sid, f)
                    molecules[s_name] = data
                    self.add_data_source(f, s_name=s_name, section="supernova-molecules")
            except RuntimeError:
                log.debug("Could not parse JSON file in {}".format(f["root"]))
                continue

        # histogram_kmer_count.json files
        for f in self.find_log_files("supernova/kmers"):
            log.debug("Found histogram_kmer_count.json in: {}".format(f["root"]))
            try:
                if f["root"] in root_summary.keys():
                    data = self.parse_histogram(f["f"], 400)
                    sid = root_summary[f["root"]]
                    s_name = self.clean_s_name(sid, f)
                    kmers[s_name] = data
                    self.add_data_source(f, s_name=s_name, section="supernova-kmers")
            except RuntimeError:
                log.debug("Could not parse JSON file in {}".format(f["root"]))
                continue

        # Data from summary.json supersedes data from report.txt
        for sample_id, sum_data in summaries.items():
            if sample_id in reports.keys():
                log.debug("Found summary data for sample {} which supersedes report data".format(sample_id))
                reports[sample_id] = sum_data
        # Ignore cmd-line specified samples
        reports = self.ignore_samples(reports)
        molecules = self.ignore_samples(molecules)
        kmers = self.ignore_samples(kmers)

        if len(reports) == 0:
            raise ModuleNoSamplesFound
        else:
            log.info("Found {} reports".format(len(reports.keys())))

        ### Write the report
        self.write_data_file(reports, "multiqc_supernova")
        config_table = {"id": "supernova_table", "namespace": "supernova"}
        self.add_section(
            name="Assembly statistics",
            anchor="supernova-table",
            description="Statistics gathered from the summary report(s) of Supernova. Note! "
            "There are more columns available but they are hidden by default.",
            helptext="As a bare minimum these numbers are generated from the file report.txt, "
            "found in the folder `sampleID/outs/`. If available the stats in the report "
            "file will be superseded by the higher precision numbers found in the file "
            "`sampleID/outs/assembly/stats/summary.json`",
            plot=table.plot(reports, self.headers, config_table),
        )

        # N50 barcharts
        n50_cats = [
            {"Scaff N50": {"name": "Scaffold N50", "color": "#66c2a5"}},
            {"Contig N50": {"name": "Contig N50", "color": "#fc8d62"}},
            {"Edge N50": {"name": "Edge N50", "color": "#8da0cb"}},
            {"Phase N50": {"name": "Phase block N50", "color": "#e78ac3"}},
        ]
        config_n50 = {
            "id": "supernova_n50",
            "title": "Supernova: N50 statistics",
            "ylab": "Scaffold N50",
            "cpswitch": False,
            "data_labels": ["Scaffold N50", "Contig N50", "Edge N50", "Phase block N50"],
        }
        self.add_section(
            name="N50 statistics",
            anchor="supernova-n50",
            description="Assembly N50 values - the shortest sequence length at 50% of the genome when sorted by size (see [wikipedia](https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics#N50)).",
            helptext="Note that assembly size and N50 values are computed after removing scaffolds &le; 10 kb and do not count `N`s: \n\n"
            "* **Scaffold N50** - N50 size of scaffolds in bases, \n"
            "* **Contig N50** - N50 size of contigs in bases, \n"
            "* **Edge N50** - N50 size of raw graph assembly edges in bases, \n"
            "* **Phase block N50** - N50 size of phase blocks in bases. \n\n"
            "[(source)](https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/output/asm-stats)",
            plot=bargraph.plot([reports, reports, reports, reports], n50_cats, config_n50),
        )

        # Conditional sections
        if len(molecules) > 0:
            # Remove the long tail, or fail if this is a legacy empty json file
            try:
                max_x = self.trim_tail(molecules, 100000)
            except IndexError:
                log.debug("The histogram file is empty. Skipping molecule length section")
                return
            # Add molecules plot
            config_molecules = {
                "id": "supernova_molecules",
                "title": "Supernova: Molecule Lengths",
                "xlab": "Inferred molecule length (bp)",
                "ylab": "# molecules",
                "smooth_points": 300,
                "smooth_points_sumcounts": True,
                "xmax": max_x,
            }
            self.add_section(
                name="Molecule Lengths",
                anchor="supernova-molecules",
                description="Shows the inferred molecule lengths of the input 10X library.",
                helptext="Inferred in the `patch` step of the Supernova pipeline. It is worth "
                "keeping in mind that the mean molecule length from the report is a length-weighted mean. "
                "See the [source code](https://github.com/10XGenomics/supernova/search?q=lw_mean_mol_len&type=) "
                "for how this value is calculated.",
                plot=linegraph.plot(molecules, config_molecules),
            )
        if len(kmers) > 0:
            # Remove the long tail, or fail if this is a legacy empty json file
            try:
                max_x = self.trim_tail(kmers, 50)
            except IndexError:
                log.debug("The histogram file is empty. Skipping kmers section")
                return
            # Add kmers plot
            config_kmers = {
                "id": "supernova_kmers",
                "title": "Supernova: Kmer Counts",
                "xlab": "Filtered kmer multiplicity",
                "ylab": "Counts",
                "smooth_points_sumcounts": False,
                "xmax": max_x,
            }
            self.add_section(
                name="K-mer counts",
                anchor="supernova-kmers",
                description="Shows the k-mer frequencies of the input data to Supernova (after filtering).",
                helptext="This data is generated from k-merizing the input read data, where the sequences are "
                "transformed in to the set of all possible sub-sequences of a fixed length of `K` (Supernova uses `K=48`). "
                "The plot shows on the x-axis the multiplicity (i.e. how many times are they repeated) of these k-mers "
                "and the y-axis the number of k-mers at this level of multiplicity. "
                "A careful reading of this plot can give some insights into the levels of heterozygosity and repeats "
                "in the genome that was sequenced and indications if the sequencing experiment was successful.",
                plot=linegraph.plot(kmers, config_kmers),
            )

    def parse_summary(self, content):
        stats = {
            "assembly_size": "Asm size",
            "bases_per_read": "Read len",
            "contig_N50": "Contig N50",
            "dup_perc": "% Dup",
            "edge_N50": "Edge N50",
            "effective_coverage": "Coverage",
            "hetdist": "Het dist",
            "lw_mean_mol_len": "Mol size",
            "median_ins_sz": "Insert size",
            "nreads": "# Reads",
            "phase_block_N50": "Phase N50",
            "placed_perc": "% Phased",
            "placed_frac": "% Phased",
            "proper_pairs_perc": "% proper",
            "q30_r2_perc": "% R2 Q30",
            "rpb_N50": "Barcode N50",
            "scaffold_N50": "Scaff N50",
            "scaffolds_10kb_plus": "# Long scaffs",
            "valid_bc_perc": "% missing BC",
            "m10": "% missing 10Kb",
            "high_AT_index": "% AT",
            "raw_coverage": "Raw coverage",
            "barcode_fraction": "BC usage",
            "repfrac": "Repeats",
            "est_genome_size": "Est size",
            "p10": "p10",
        }

        try:
            cdict = json.loads(content)
        except ValueError as e:
            raise e

        data = {}
        # Try to find sample_id
        sid = ""
        if "CS_SAMPLE_ID" in cdict.keys():
            sid = cdict["CS_SAMPLE_ID"]  # supernova 1.2
        elif "sample_id" in cdict.keys():
            sid = cdict["sample_id"]
        else:
            raise RuntimeError

        for key, value in cdict.items():
            if key in stats.keys():
                # Some trickery for supernova 1.1.4 compatability
                if key == "placed_frac":
                    value = value * 100
                if key == "valid_bc_perc":
                    value = 100 - value
                data[stats[key]] = value

        return (sid, data)

    def parse_report(self, content):
        # Some short-hands for converting the report numbers (pi is exactly three!)
        exp = {
            "K": 1000.0,
            "Kb": 1000.0,
            "kb": 1000.0,
            "Mb": 1000000.0,
            "M": 1000000.0,
            "Gb": 1000000000.0,
            "G": 1000000000.0,
        }
        stats = {
            "READS": "# Reads",
            "MEAN READ LEN": "Read len",
            "EFFECTIVE COV": "Coverage",
            "READ TWO Q30": "% R2 Q30",
            "MEDIAN INSERT": "Insert size",
            "PROPER PAIRS": "% proper",
            "MOLECULE LEN": "Mol size",
            "HETDIST": "Het dist",
            "UNBAR": "% missing BC",
            "BARCODE N50": "Barcode N50",
            "DUPS": "% Dup",
            "PHASED": "% Phased",
            "LONG SCAFFOLDS": "# Long scaffs",
            "EDGE N50": "Edge N50",
            "CONTIG N50": "Contig N50",
            "PHASEBLOCK N50": "Phase N50",
            "SCAFFOLD N50": "Scaff N50",
            "ASSEMBLY SIZE": "Asm size",
            "MISSING 10KB": "% missing 10Kb",
            "HIGH AT FRACTION": "% AT",
            "RAW COV": "Raw coverage",
            "BARCODE FRACTION": "BC usage",
            "REPETITIVE FRAC": "Repeats",
            "EST GENOME SIZE": "Est size",
            "P10": "p10",
        }

        data = {}
        # Find the sample ID
        sid = ""
        sid_pat = re.compile("- \[(.+)\]")
        # [number, unit, category]
        stat_pat = re.compile("-\s+(\d+\.\d+)\s+(\S+|.)\s+= (.+) =")

        for l in content.splitlines():
            sid_m = re.match(sid_pat, l)
            stat_m = re.match(stat_pat, l)

            if sid_m is not None:
                sid = sid_m.groups()[0]
            if stat_m is None:
                continue
            stat_val = stat_m.groups()
            stat_type = stat_val[2].strip()
            # Parse the lines containing statistics
            if stat_type in stats.keys():
                try:
                    if stat_val[1] in exp.keys():
                        data[stats[stat_type]] = float(stat_val[0]) * exp[stat_val[1]]
                    else:
                        data[stats[stat_type]] = float(stat_val[0])
                except ValueError:
                    log.debug('Error in parsing sample {}, on line "{}"'.format(sid, stat_val))

        return (sid, data)

    def parse_histogram(self, content, cutoff=None):
        try:
            cdict = json.loads(content)
        except ValueError as e:
            raise e

        numbins = cdict["numbins"] + 1
        xdata = [i * cdict["binsize"] for i in range(0, numbins)]
        return {i: j for (i, j) in zip(xdata, cdict["vals"][:cutoff])}

    def trim_tail(self, plot, min_x=50, pct=0.99):
        join_plot = {}
        cuml_plot = {}
        for sample, plot_data in plot.items():
            for key, value in plot_data.items():
                join_plot[key] = join_plot.get(key, 0) + value
        max_i = 0
        for key in join_plot.keys():
            max_i += join_plot[key]
            cuml_plot[key] = max_i
        max_x = [i for i, j in cuml_plot.items() if j <= max_i * pct][-1]
        # xlim = {, 50} at minimum
        if max_x < min_x:
            max_x = min_x
        return max_x
