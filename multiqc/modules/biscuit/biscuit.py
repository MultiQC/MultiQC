""" MultiQC module to parse output from BISCUITqc """


import logging
import re

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, violin

# Initialize the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Inherits from base module class
    Initialize data structures and prepare BISCUIT report
    Inputs:
        No inputs
    Returns:
        BISCUIT report for MultiQC
    """

    def __init__(self):
        # Initialize the parent object
        super(MultiqcModule, self).__init__(
            name="BISCUIT",
            anchor="biscuit",
            href="https://github.com/huishenlab/biscuit",
            info="is a tool to map bisulfite converted DNA sequence reads and determine cytosine methylation states.",
            # Can't find a DOI // doi=
        )

        # Set up data structures
        self.mdata = {
            # General statistics
            "align_mapq": {},
            "align_strand": {},
            "align_isize": {},
            # Duplicate reporting
            "dup_report": {},
            # Uniformity
            "qc_cv": {},
            # Base coverage
            "covdist_all_base_botgc": {},
            "covdist_all_base": {},
            "covdist_all_base_topgc": {},
            "covdist_q40_base_botgc": {},
            "covdist_q40_base": {},
            "covdist_q40_base_topgc": {},
            # CpG coverage
            "covdist_all_cpg_botgc": {},
            "covdist_all_cpg": {},
            "covdist_all_cpg_topgc": {},
            "covdist_q40_cpg_botgc": {},
            "covdist_q40_cpg": {},
            "covdist_q40_cpg_topgc": {},
            # Cytosine retention
            "cpg_retention_readpos": {},
            "cph_retention_readpos": {},
            "base_avg_retention_rate": {},
            "read_avg_retention_rate": {},
        }

        # NB: Cleaning filenames like this means that some MultiQC functionality, like -s / --fullnames doesn't work.
        # However, because some of the parsing relies on cleaned filenames, the above options break the
        # module if we use the centralised MultiQC functions.
        file_suffixes = [
            # General statistics
            ".txt",
            "_mapq_table",
            "_strand_table",
            "_isize_table",
            # Duplicate reporting
            "_dup_report",
            # Uniformity
            "_cv_table",
            # Base coverage
            "_covdist_all_base_botgc_table",
            "_covdist_all_base_table",
            "_covdist_all_base_topgc_table",
            "_covdist_q40_base_botgc_table",
            "_covdist_q40_base_table",
            "_covdist_q40_base_topgc_table",
            # CpG coverage
            "_covdist_all_cpg_botgc_table",
            "_covdist_all_cpg_table",
            "_covdist_all_cpg_topgc_table",
            "_covdist_q40_cpg_botgc_table",
            "_covdist_q40_cpg_table",
            "_covdist_q40_cpg_topgc_table",
            # Cytosine retention
            "_CpGRetentionByReadPos",
            "_CpHRetentionByReadPos",
            "_totalBaseConversionRate",
            "_totalReadConversionRate",
        ]

        # Find and parse alignment reports
        for k in self.mdata:
            for f in self.find_log_files(f"biscuit/{k}"):
                s_name = f["fn"]
                for suffix in file_suffixes:
                    s_name = s_name.replace(suffix, "")
                s_name = self.clean_s_name(s_name, f)

                # Add source file to multiqc_sources.txt
                self.add_data_source(f, s_name=s_name, section=k)

                if s_name in self.mdata[k]:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")

                self.mdata[k][s_name] = getattr(self, f"parse_logs_{k}")(f["f"], f["fn"])

        for k in self.mdata:
            self.mdata[k] = self.ignore_samples(self.mdata[k])

        n_samples = max([len(self.mdata[k]) for k in self.mdata])

        if n_samples == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {n_samples} samples")

        # Basic stats table
        self.biscuit_stats_table()

        # Write data to file
        self.write_data_file(self.mdata, "biscuit")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Make report sections
        for k in self.mdata:
            if len(self.mdata[k]) > 0:
                log.debug(f"Found {len(self.mdata[k])} {k} reports")
                getattr(self, f"chart_{k}")()

    def biscuit_stats_table(self):
        """
        Create general statistics table for BISCUIT data
        Inputs:
            Uses mdata['align_mapq'] and mdata['dup_report']
        Returns:
            Add columns to MultiQC general statistics table
        """
        pd = {}

        # Calculate % aligned
        for s_name, dd in self.mdata["align_mapq"].items():
            if len(dd) > 0:
                pd[s_name] = {"aligned": dd["frc_align"]}

        # Calculate % duplicated
        for s_name, dd in self.mdata["dup_report"].items():
            if s_name not in pd:
                pd[s_name] = {}
            if "all" in dd and dd["all"] != -1:
                pd[s_name]["dup_all"] = dd["all"]
            if "q40" in dd and dd["q40"] != -1:
                pd[s_name]["dup_q40"] = dd["q40"]

        pheader = {
            "dup_q40": {
                "title": "Dup. % for Q40 Reads",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "YlOrBr",
                "hidden": True,
            },
            "dup_all": {"title": "Dup. % for All Reads", "min": 0, "max": 100, "suffix": "%", "scale": "Reds"},
            "aligned": {
                "title": "% Aligned",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
            },
        }

        self.general_stats_addcols(pd, pheader)

    ########################################
    #####  General Mapping Information #####
    ########################################
    @staticmethod
    def parse_logs_align_mapq(f, fn):
        """
        Parse _mapq_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of aligned mapq data
        """
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug(f"No data available in {fn}. Will not fill corresponding entries.")
            return {}

        mapq = {}
        for line in file_data:
            s = line.split()
            mapq[s[0]] = s[1]  # mapq[MAPQ] = number of reads

        data = {
            "frc_align": 0,
            "opt_align": 0,
            "sub_align": 0,
            "not_align": 0,
            "mapqs": dict(zip(range(61), [0 for _ in range(61)])),
        }
        if len(mapq) > 0:
            total = sum([int(cnt) for _, cnt in mapq.items() if _ != "unmapped"])
            for mq, cnt in mapq.items():
                if mq == "unmapped":
                    data["not_align"] += int(cnt)
                else:
                    data["mapqs"][int(mq)] = 100.0 * float(cnt) / total
                    if int(mq) >= 40:
                        data["opt_align"] += int(cnt)
                    else:
                        data["sub_align"] += int(cnt)

            data["frc_align"] = (
                100
                * (data["opt_align"] + data["sub_align"])
                / (data["opt_align"] + data["sub_align"] + data["not_align"])
            )

        return data

    ########################################
    ####    Alignment Quality Report    ####
    ########################################
    def chart_align_mapq(self):
        """
        Chart _mapq_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Mapping Overview and Mapping Quality
            Distribution charts
        """

        #
        # Mapping Overview bar chart
        #
        pd = {}

        # Calculate alignment counts
        for s_name, dd in self.mdata["align_mapq"].items():
            if len(dd) > 0:
                pd[s_name] = {"opt_align": dd["opt_align"], "sub_align": dd["sub_align"], "not_align": dd["not_align"]}

        pheader = {
            "opt_align": {"color": "#1f78b4", "name": "Optimally Aligned Reads"},
            "sub_align": {"color": "#a6cee3", "name": "Suboptimally Aligned Reads"},
            "not_align": {"color": "#b2df8a", "name": "Unaligned Reads"},
        }

        pconfig = {
            "id": "biscuit-mapping-overview-plot",
            "title": "BISCUIT: Mapping Overview",
            "ylab": "Number of Reads",
            "cpswitch_counts_label": "# Reads",
        }

        self.add_section(
            name="Mapping Overview",
            anchor="biscuit-mapping-overview",
            description="""
                Number of optimally aligned reads (`MAPQ>=40`), suboptimally
                aligned reads (`MAPQ<40`), and unmapped reads. Primary alignments only.
            """,
            helptext="""
                A good library should have a high fraction of reads
                that are optimally aligned. Note, suboptimally aligned reads
                include both non-unique alignments and imperfect alignments.
            """,
            plot=bargraph.plot(pd, pheader, pconfig),
        )

        #
        # Mapping Quality Distribution line graph
        #

        # Calculate the % aligned for each mapping q score
        pd_mapq = {}
        for s_name, dd in self.mdata["align_mapq"].items():
            if len(dd) > 0:
                pd_mapq[s_name] = dd["mapqs"]

        pconfig = {
            "id": "biscuit_mapq",
            "title": "BISCUIT: Distribution of Mapping Qualities",
            "ymin": 0,
            "xmin": 0,
            "tt_label": "<strong>Q{point.x}:</strong> {point.y:.2f}% of mapped reads",
            "ysuffix": "%",
            "ylab": "% of primary mapped reads",
            "xlab": "Mapping quality score",
        }

        self.add_section(
            name="Mapping Quality Distribution",
            anchor="biscuit-mapq",
            description="""
                The percentage of the total number of mapped reads
                for each mapping quality score. Primary alignments only.
            """,
            helptext="""
                A good quality sample should have a high quality mapping score
                for the majority of alignments.
            """,
            plot=linegraph.plot(pd_mapq, pconfig),
        )

    ########################################
    ####     Strand Alignment Report    ####
    ########################################
    def parse_logs_align_strand(self, f, fn):
        """
        Parse _strand_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of strand data for reads 1 and 2
        """
        patterns = [
            r"(R1)\s+\((f)\)\:\s+(\d+)\s+(\d+)",
            r"(R1)\s+\((r)\)\:\s+(\d+)\s+(\d+)",
            r"(R2)\s+\((f)\)\:\s+(\d+)\s+(\d+)",
            r"(R2)\s+\((r)\)\:\s+(\d+)\s+(\d+)",
        ]

        data = {"read1": {}, "read2": {}}
        for pat in patterns:
            m = re.search(pat, f, re.MULTILINE)
            if m is not None:
                if m.group(1) == "R1":
                    if m.group(2) == "f":
                        data["read1"]["ff"] = int(m.group(3))
                        data["read1"]["fr"] = int(m.group(4))
                    else:
                        data["read1"]["rf"] = int(m.group(3))
                        data["read1"]["rr"] = int(m.group(4))
                else:
                    if m.group(2) == "f":
                        data["read2"]["ff"] = int(m.group(3))
                        data["read2"]["fr"] = int(m.group(4))
                    else:
                        data["read2"]["rf"] = int(m.group(3))
                        data["read2"]["rr"] = int(m.group(4))

        return data

    def chart_align_strand(self):
        """
        Chart _strand_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Mapping Strand Distribution chart
        """

        pd1 = {}
        pd2 = {}
        for s_name, dd in self.mdata["align_strand"].items():
            if len(dd["read1"]) > 0:
                pd1[s_name] = dd["read1"]
            if len(dd["read2"]) > 0:
                pd2[s_name] = dd["read2"]

        pheader = {
            "ff": {"color": "#F53855", "name": "ff: Watson-Aligned, Watson-Bisulfite Conversion"},
            "fr": {"color": "#E37B40", "name": "fr: Watson-Aligned, Crick-Bisulfite Conversion"},
            "rf": {"color": "#46B29D", "name": "rf: Crick-Aligned, Watson-Bisulfite Conversion"},
            "rr": {"color": "#324D5C", "name": "rr: Crick-Aligned, Crick-Bisulfite Conversion"},
        }
        pconfig = {
            "id": "biscuit_strands",
            "title": "BISCUIT: Mapping Strand Distribution",
            "ylab": "Number of Reads",
            "cpswitch_counts_label": "# Reads",
            "data_labels": [{"name": "Read 1"}, {"name": "Read 2"}],
        }

        # TODO: When PBAT mode is implemented, add comment in help text about
        #       how to interpret PBAT mode results
        if max(len(pd1), len(pd2)) > 0:
            self.add_section(
                name="Mapping Strand Distribution",
                anchor="biscuit-strands",
                description="For primary alignments, shows the number of reads mapped to each strand.",
                helptext="""
                    Most bisulfite libraries typically map Read 1 to the parent
                    strand (`ff`, `rr`) and Read 2 to the daughter / synthesized
                    strand (`fr`, `rf`).

                    Note that PBAT and many single-cell / low input
                    libraries may not follow this assumption.
                """,
                plot=bargraph.plot([pd1, pd2], [pheader, pheader], pconfig),
            )

    ########################################
    ####       Insert Size Report       ####
    ########################################
    @staticmethod
    def parse_logs_align_isize(f, fn):
        """
        Parse _isize_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of insert size data
        """
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug(f"No data available in {fn}. Will not fill corresponding entries.")
            return {"no_data_available": 1}

        data = {"percent": {}, "readcnt": {}}
        for line in file_data:
            fields = line.split("\t")
            data["percent"][int(fields[0])] = 100.0 * float(fields[1])
            data["readcnt"][int(fields[0])] = float(fields[2])

        return data

    def chart_align_isize(self):
        """
        Chart _isize_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Insert Size Distribution chart
        """

        pd_p = {}
        pd_r = {}
        for s_name, dd in self.mdata["align_isize"].items():
            if "no_data_available" not in dd.keys():
                pd_p[s_name] = dd["percent"]
                pd_r[s_name] = dd["readcnt"]

        pconfig = {
            "id": "biscuit_isize",
            "title": "BISCUIT: Insert Size Distribution",
            "ymin": 0,
            "xmin": 0,
            "smooth_points": 1000,  # limit number of points / smooth data
            "xlab": "Insert Size",
            "data_labels": [
                {
                    "name": "% of Mapped Reads",
                    "ylab": "% of Mapped Reads",
                    "ysuffix": "%",
                    "tt_label": "<strong>IS%{x} bp:</strong> %{y:.2f}%",
                },
                {
                    "name": "Mapped Reads",
                    "ylab": "Mapped Reads",
                    "ysuffix": "",
                    "tt_label": "<strong>IS%{x} bp:</strong> %{y:,.0f}",
                },
            ],
        }

        self.add_section(
            name="Insert Size Distribution",
            anchor="biscuit-isize",
            description="Shows the distribution of insert sizes.",
            helptext="""
                Insert size is defined as:

                ```
                (right-most coordinate of reverse-mate read) - (left-most coordinate of forward-mate read)
                ```

                Insert sizes are calculated for reads with a _"mapped in
                proper pair"_ `samtools` flag, and `MAPQ >= 40`.
            """,
            plot=linegraph.plot([pd_p, pd_r], pconfig),
        )

    ########################################
    ####        Duplicate Report        ####
    ########################################
    @staticmethod
    def parse_logs_dup_report(f, fn):
        """
        Parses _dup_report.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of duplicate fractions
        """
        patterns = [
            (r"Number of duplicate reads:\s+(\d+)", r"Number of reads:\s+(\d+)", "all"),
            (r"Number of duplicate q40-reads:\s+(\d+)", r"Number of q40-reads:\s+(\d+)", "q40"),
        ]

        data = {}
        for pat_dup, pat_tot, key in patterns:
            m1 = re.search(pat_dup, f, re.MULTILINE)
            m2 = re.search(pat_tot, f, re.MULTILINE)
            if m1 is not None and m2 is not None:
                data[key] = 100.0 * float(m1.group(1)) / float(m2.group(1))
            else:
                data[key] = -1

        return data

    def chart_dup_report(self):
        """
        Charts _dup_report.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Duplicate Rates chart
        """

        pd1 = {}  # Overall duplicate rate
        pd2 = {}  # MAPQ>=40 duplicate rate
        for s_name, dd in self.mdata["dup_report"].items():
            if "all" in dd and dd["all"] != -1:
                pd1[s_name] = {"dup_rate": dd["all"]}
            if "q40" in dd and dd["q40"] != -1:
                pd2[s_name] = {"dup_rate": dd["q40"]}

        pheader = {
            "dup_rate": {"color": "#a50f15", "name": "Duplicate Rate"},
        }
        pconfig = {
            "id": "biscuit_dup_report",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "BISCUIT: Percentage of Duplicate Reads",
            "data_labels": [{"name": "Overall Duplicate Rate"}, {"name": "MAPQ>=40 Duplicate Rate"}],
            "ylab": "Duplicate rate",
            "ymin": 0,
            "ymax": 100,
            "y_clipmax": 110,
            "use_legend": False,
            "tt_decimals": 1,
            "tt_suffix": "%",
        }

        if len(pd1) > 0:
            self.add_section(
                name="Duplicate Rates",
                anchor="biscuit-dup-report",
                description="Shows the percentage of total reads that are duplicates.",
                helptext="""
                    `MAPQ >= 40` shows the duplicate rate for just the reads
                    with a mapping quality score of `MAPQ >= 40`.
                """,
                plot=bargraph.plot([pd1, pd2], [pheader, pheader], pconfig),
            )

    ########################################
    ####      Depths and Uniformity     ####
    ########################################
    @staticmethod
    def parse_logs_qc_cv(f, fn):
        """
        Parses _cv_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of depth uniformity measures
        """

        data = {}
        targets = [
            "all_base",
            "all_cpg",
            "q40_base",
            "q40_cpg",
            "all_base_botgc",
            "all_cpg_botgc",
            "q40_base_botgc",
            "q40_cpg_botgc",
            "all_base_topgc",
            "all_cpg_topgc",
            "q40_base_topgc",
            "q40_cpg_topgc",
        ]
        for t in targets:
            m = re.search(rf"{t}\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)", f, re.MULTILINE)
            if m is not None:
                data[t] = {"mu": float(m.group(1)), "sigma": float(m.group(2)), "cv": float(m.group(3))}
            else:
                data[t] = {"mu": -1, "sigma": -1, "cv": -1}

        return data

    def chart_qc_cv(self):
        """
        Charts _cv_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Sequencing Depth - Whole Genome chart
        """

        cats = [
            ("all_base", "a_b"),
            ("q40_base", "q_b"),
            ("all_base_botgc", "a_b_b"),
            ("q40_base_botgc", "q_b_b"),
            ("all_base_topgc", "a_b_t"),
            ("q40_base_topgc", "q_b_t"),
            ("all_cpg", "a_c"),
            ("q40_cpg", "q_c"),
            ("all_cpg_botgc", "a_c_b"),
            ("q40_cpg_botgc", "q_c_b"),
            ("all_cpg_topgc", "a_c_t"),
            ("q40_cpg_topgc", "q_c_t"),
        ]

        pd = dict()
        for s_name, dd in self.mdata["qc_cv"].items():
            data = dict()
            for cat, key in cats:
                if cat in dd:
                    if dd[cat]["mu"] != -1:
                        data["mu_" + key] = dd[cat]["mu"]
                        data["cv_" + key] = dd[cat]["cv"]
            if len(data) > 0:
                pd[s_name] = data

        shared_mean = {"min": 0, "format": "{:,3f}", "minrange": 10}
        shared_cofv = {"min": 0, "format": "{:,3f}", "minrange": 50}

        pheader = {
            "mu_a_b": dict(
                shared_mean, **{"title": "All Genome Mean", "description": "Mean Sequencing Depth for All Reads"}
            ),
            "mu_q_b": dict(
                shared_mean, **{"title": "Q40 Genome Mean", "description": "Mean Sequencing Depth for Q40 Reads"}
            ),
            "mu_a_b_b": dict(
                shared_mean,
                **{
                    "title": "Low GC All Gen. Mean",
                    "description": "Mean Sequencing Depth for All Reads in Low GC-Content Regions",
                },
            ),
            "mu_q_b_b": dict(
                shared_mean,
                **{
                    "title": "Low GC Q40 Gen. Mean",
                    "description": "Mean Sequencing Depth for Q40 Reads in Low GC-Content Regions",
                },
            ),
            "mu_a_b_t": dict(
                shared_mean,
                **{
                    "title": "High GC All Gen. Mean",
                    "description": "Mean Sequencing Depth for All Reads in High GC-Content Regions",
                },
            ),
            "mu_q_b_t": dict(
                shared_mean,
                **{
                    "title": "High GC Q40 Gen. Mean",
                    "description": "Mean Sequencing Depth for Q40 Reads in High GC-Content Regions",
                },
            ),
            "cv_a_b": dict(
                shared_cofv, **{"title": "All Genome CoV", "description": "Sequencing Depth CoV for All Reads"}
            ),
            "cv_q_b": dict(
                shared_cofv, **{"title": "Q40 Genome CoV", "description": "Sequencing Depth CoV for Q40 Reads"}
            ),
            "cv_a_b_b": dict(
                shared_cofv,
                **{
                    "title": "Low GC All Gen. CoV",
                    "description": "Sequencing Depth CoV for All Reads in Low GC-Content Regions",
                },
            ),
            "cv_q_b_b": dict(
                shared_cofv,
                **{
                    "title": "Low GC Q40 Gen. CoV",
                    "description": "Sequencing Depth CoV for Q40 Reads in Low GC-Content Regions",
                },
            ),
            "cv_a_b_t": dict(
                shared_cofv,
                **{
                    "title": "High GC All Gen. CoV",
                    "description": "Sequencing Depth CoV for All Reads in High GC-Content Regions",
                },
            ),
            "cv_q_b_t": dict(
                shared_cofv,
                **{
                    "title": "High GC Q40 Gen. CoV",
                    "description": "Sequencing Depth CoV for Q40 Reads in High GC-Content Regions",
                },
            ),
            "mu_a_c": dict(
                shared_mean, **{"title": "All CpGs Mean", "description": "Mean Sequencing Depth for All CpGs"}
            ),
            "mu_q_c": dict(
                shared_mean, **{"title": "Q40 CpGs Mean", "description": "Mean Sequencing Depth for Q40 CpGs"}
            ),
            "mu_a_c_b": dict(
                shared_mean,
                **{
                    "title": "Low GC All CpGs Mean",
                    "description": "Mean Sequencing Depth for All CpGs in Low GC-Content Regions",
                },
            ),
            "mu_q_c_b": dict(
                shared_mean,
                **{
                    "title": "Low GC Q40 CpGs Mean",
                    "description": "Mean Sequencing Depth for Q40 CpGs in Low GC-Content Regions",
                },
            ),
            "mu_a_c_t": dict(
                shared_mean,
                **{
                    "title": "High GC All CpGs Mean",
                    "description": "Mean Sequencing Depth for All CpGs in High GC-Content Regions",
                },
            ),
            "mu_q_c_t": dict(
                shared_mean,
                **{
                    "title": "High GC Q40 CpGs Mean",
                    "description": "Mean Sequencing Depth for Q40 CpGs in High GC-Content Regions",
                },
            ),
            "cv_a_c": dict(
                shared_cofv, **{"title": "All CpGs CoV", "description": "Sequencing Depth CoV for All CpGs"}
            ),
            "cv_q_c": dict(
                shared_cofv, **{"title": "Q40 CpGs CoV", "description": "Sequencing Depth CoV for Q40 CpGs"}
            ),
            "cv_a_c_b": dict(
                shared_cofv,
                **{
                    "title": "Low GC All CpGs CoV",
                    "description": "Sequencing Depth CoV for All CpGs in Low GC-Content Regions",
                },
            ),
            "cv_q_c_b": dict(
                shared_cofv,
                **{
                    "title": "Low GC Q40 CpGs CoV",
                    "description": "Sequencing Depth CoV for Q40 CpGs in Low GC-Content Regions",
                },
            ),
            "cv_a_c_t": dict(
                shared_cofv,
                **{
                    "title": "High GC All CpGs CoV",
                    "description": "Sequencing Depth CoV for All CpGs in High GC-Content Regions",
                },
            ),
            "cv_q_c_t": dict(
                shared_cofv,
                **{
                    "title": "High GC Q40 CpGs CoV",
                    "description": "Sequencing Depth CoV for Q40 CpGs in High GC-Content Regions",
                },
            ),
        }
        pconfig = {
            "id": "biscuit_seq_depth",
            "table_title": "BISCUIT: Sequencing Depth",
            "sort_rows": False,
        }

        if len(pd) > 0:
            self.add_section(
                name="Sequencing Depth Statistics",
                anchor="biscuit-seq-depth",
                description="""
                    Shows the sequence depth mean and uniformity measured by the Coefficient of Variation
                    (`CoV`, defined as `stddev/mean`).
                """,
                helptext="""
                    The plot shows coverage across different selections:

                    * _Genome_ (Gen.) - Statistics for all bases across the entire genome
                    * _CpGs_ - Statistics for CpGs
                    * _All_ - Statistics for any mapped bases/CpGs
                    * _Q40_ - Statistics only those bases/CpGs with mapping quality `MAPQ >= 40`
                    * _High GC_ - Bases / CpGs that overlap with the top 10% of 100bp windows for GC-content
                    * _Low GC_ - Bases / CpGs that overlap with the bottom 10% of 100bp windows for GC-content

                """,
                plot=violin.plot(pd, pheader, pconfig),
            )

    ########################################
    #### Base Coverage and CpG Coverage ####
    ########################################
    @staticmethod
    def parse_logs_covdist_all_base(f, fn):
        """
        Parses _covdist_all_base_botgc_table.txt
               _covdist_all_base_table.txt
               _covdist_all_base_topgc_table.txt
               _covdist_all_cpg_botgc_table.txt
               _covdist_all_cpg_table.txt
               _covdist_all_cpg_topgc_table.txt
               _covdist_q40_base_botgc_table.txt
               _covdist_q40_base_table.txt
               _covdist_q40_base_topgc_table.txt
               _covdist_q40_cpg_botgc_table.txt
               _covdist_q40_cpg_table.txt
               _covdist_q40_cpg_topgc_table.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of coverage distributions up to 30X data
        """
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug(f"No data available in {fn}. Will not fill corresponding entries.")
            return dict(zip([i for i in range(31)], [-1 for _ in range(31)]))

        dd = {}
        for line in file_data:
            fields = line.split()
            dd[int(float(fields[0]))] = int(float(fields[1]))

        covs = sorted([k for k in dd])[:31]
        _ccov_cnt = sum(dd.values())

        ccov_cnts = []
        for cov in covs:
            ccov_cnts.append(_ccov_cnt / 1000000.0)
            _ccov_cnt -= dd[cov]

        return dict(zip(covs, ccov_cnts))

    def parse_logs_covdist_all_base_botgc(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_all_base_topgc(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_base(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_base_botgc(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_base_topgc(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_all_cpg(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_all_cpg_botgc(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_all_cpg_topgc(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_cpg(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_cpg_botgc(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def parse_logs_covdist_q40_cpg_topgc(self, f, fn):
        """Handled by parse_logs_covdist_all_base()"""
        return self.parse_logs_covdist_all_base(f, fn)

    def chart_covdist_all_base(self):
        """
        Charts _covdist_all_base_botgc_table.txt
               _covdist_all_base_table.txt
               _covdist_all_base_topgc_table.txt
               _covdist_all_cpg_botgc_table.txt
               _covdist_all_cpg_table.txt
               _covdist_all_cpg_topgc_table.txt
               _covdist_q40_base_botgc_table.txt
               _covdist_q40_base_table.txt
               _covdist_q40_base_topgc_table.txt
               _covdist_q40_cpg_botgc_table.txt
               _covdist_q40_cpg_table.txt
               _covdist_q40_cpg_topgc_table.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Cumulative Coverage chart
        """

        pd = [
            self.mdata["covdist_all_base"],
            self.mdata["covdist_q40_base"],
            self.mdata["covdist_all_cpg"],
            self.mdata["covdist_q40_cpg"],
            self.mdata["covdist_all_base_botgc"],
            self.mdata["covdist_q40_base_botgc"],
            self.mdata["covdist_all_cpg_botgc"],
            self.mdata["covdist_q40_cpg_botgc"],
            self.mdata["covdist_all_base_topgc"],
            self.mdata["covdist_q40_base_topgc"],
            self.mdata["covdist_all_cpg_topgc"],
            self.mdata["covdist_q40_cpg_topgc"],
        ]

        pconfig = {
            "id": "biscuit_cumulative",
            "title": "BISCUIT: Cumulative Coverage",
            "ymin": 0,
            "tt_label": "<strong>{point.x}X:</strong> {point.y:.2f}M",
            "xlab": "Coverage",
            "ylab": "Millions of Bases",
            "data_labels": [
                {"name": "All Bases", "ylab": "Millions of Bases"},
                {"name": "Q40 Bases", "ylab": "Millions of Bases"},
                {"name": "All CpGs", "ylab": "Millions of CpGs"},
                {"name": "Q40 CpGs", "ylab": "Millions of CpGs"},
                {"name": "Low GC All Bases", "ylab": "Millions of Bases"},
                {"name": "Low GC Q40 Bases", "ylab": "Millions of Bases"},
                {"name": "Low GC All CpGs", "ylab": "Millions of CpGs"},
                {"name": "Low GC Q40 CpGs", "ylab": "Millions of CpGs"},
                {"name": "High GC All Bases", "ylab": "Millions of Bases"},
                {"name": "High GC Q40 Bases", "ylab": "Millions of Bases"},
                {"name": "High GC All CpGs", "ylab": "Millions of CpGs"},
                {"name": "High GC Q40 CpGs", "ylab": "Millions of CpGs"},
            ],
        }

        self.add_section(
            name="Cumulative Coverage",
            anchor="biscuit-cumulative-coverage",
            description="Shows the number of bases or CpGs covered by a given number of reads.",
            helptext="""
                * _All_ - Coverage for any mapped reads
                * _Q40_ - Coverage for reads with mapping quality `MAPQ >= 40`
                * _High GC_ - Coverage for reads that overlap with the top 10% of 100bp windows for GC-content
                * _Low GC_ - Coverage for reads that overlap with the bottom 10% of 100bp windows for GC-content
            """,
            plot=linegraph.plot(pd, pconfig),
        )

    def chart_covdist_all_base_botgc(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_all_base_topgc(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_q40_base(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_q40_base_botgc(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_q40_base_topgc(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_all_cpg(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_all_cpg_botgc(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_all_cpg_topgc(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_q40_cpg(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_q40_cpg_botgc(self):
        """Handled by chart_covdist_all_base()"""
        pass

    def chart_covdist_q40_cpg_topgc(self):
        """Handled by chart_covdist_all_base()"""
        pass

    ########################################
    ####          CpG Retention         ####
    ########################################
    def parse_logs_cpg_retention_readpos(self, f, fn):
        """
        Parses _CpGRetentionByReadPos.txt
               _CpHRetentionByReadPos.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of fraction of retained cytosines for reads 1 and 2
            in either a CpH or CpG context
        """
        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug(f"No data available in {fn}. Will not fill corresponding entries.")
            return {"no_data_available": 1}

        r1 = {"C": {}, "R": {}}
        r2 = {"C": {}, "R": {}}
        for line in file_data:
            fields = line.strip().split("\t")

            if fields[0] not in ["1", "2"] or fields[2] not in ["C", "R"]:
                return {}
            if fields[0] == "1":
                r1[fields[2]][int(fields[1])] = int(fields[3])
            elif fields[0] == "2":
                r2[fields[2]][int(fields[1])] = int(fields[3])

        r1rate = dict()
        for k in sorted(r1["C"].keys()):
            if k in r1["R"]:
                r1rate[k] = 100.0 * float(r1["R"][k]) / (r1["R"][k] + r1["C"][k])

        r2rate = dict()
        for k in sorted(r2["C"].keys()):
            if k in r2["R"]:
                r2rate[k] = 100.0 * float(r2["R"][k]) / (r2["R"][k] + r2["C"][k])

        return {"1": r1rate, "2": r2rate}

    def chart_cpg_retention_readpos(self):
        """
        Charts _CpGRetentionByReadPos.txt
               _CpHRetentionByReadPos.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Retenion vs. Base Position in Read chart
        """

        pd = [
            dict(
                [
                    (s_name, dd["1"])
                    for s_name, dd in self.mdata["cpg_retention_readpos"].items()
                    if "no_data_available" not in dd.keys()
                ]
            ),
            dict(
                [
                    (s_name, dd["2"])
                    for s_name, dd in self.mdata["cpg_retention_readpos"].items()
                    if "no_data_available" not in dd.keys()
                ]
            ),
            dict(
                [
                    (s_name, dd["1"])
                    for s_name, dd in self.mdata["cph_retention_readpos"].items()
                    if "no_data_available" not in dd.keys()
                ]
            ),
            dict(
                [
                    (s_name, dd["2"])
                    for s_name, dd in self.mdata["cph_retention_readpos"].items()
                    if "no_data_available" not in dd.keys()
                ]
            ),
        ]

        pconfig = {
            "id": "biscuit_retention_cytosine",
            "title": "BISCUIT: Retention vs. Base Position in Read",
            "xlab": "Position in Read",
            "xsuffix": "bp",
            "ylab": "CpG Retention Rate (%)",
            "ymin": 0,
            "ymax": 100,
            "y_minrange": 0,
            "y_clipmin": 0,
            "tt_label": "<strong>Position {point.x}:</strong> {point.y:.2f}%",
            "data_labels": [
                {"name": "CpG Read 1", "ylab": "CpG Retention Rate (%)"},
                {"name": "CpG Read 2", "ylab": "CpG Retention Rate (%)"},
                {"name": "CpH Read 1", "ylab": "CpH Retention Rate (%)"},
                {"name": "CpH Read 2", "ylab": "CpH Retention Rate (%)"},
            ],
        }

        self.add_section(
            name="Retention vs. Base Position in Read",
            anchor="biscuit-retention-cytosine",
            description="Distribution of cytosine retention rates across base positions in the read (a.k.a. _M-bias_ plot).",
            plot=linegraph.plot(pd, pconfig),
        )

    def parse_logs_cph_retention_readpos(self, f, fn):
        """Handled by parse_logs_cpg_retention_readpos()"""
        return self.parse_logs_cpg_retention_readpos(f, fn)

    def chart_cph_retention_readpos(self):
        """Handled by chart_cpg_retention_readpos()"""
        pass

    def parse_logs_read_avg_retention_rate(self, f, fn):
        """
        Parses _totalReadConversionRate.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of read averaged fraction of retainied cytosines by context
        """

        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug(f"No data available in {fn}. Will not fill corresponding entries.")
            return {"no_data_available": 1}

        data = {}
        for line in file_data:
            fields = line.split("\t")
            # Skip rows that have NaNs as something went wrong in processing
            if "nan" in fields:
                log.debug(f"Found NaN in {fn}. Skipping.")
                continue

            # BISCUIT returns -1 if insufficient data. Only fill fields with value >= 0.
            if float(fields[0]) >= 0:
                data["rca"] = 100.0 * float(fields[0])
            if float(fields[1]) >= 0:
                data["rcc"] = 100.0 * float(fields[1])
            if float(fields[2]) >= 0:
                data["rcg"] = 100.0 * float(fields[2])
            if float(fields[3]) >= 0:
                data["rct"] = 100.0 * float(fields[3])

        return data

    def chart_read_avg_retention_rate(self):
        """
        Charts _totalReadConversionRate.txt
               _totalBaseConversionRate.txt
        Inputs:
            No inputs
        Returns:
            No returns, generates Cytosine Retention chart
        """

        pdata_byread = {}
        for s_name, dd in self.mdata["read_avg_retention_rate"].items():
            if "no_data_available" not in dd.keys():
                pdata_byread[s_name] = dd

        pdata_bybase = {}
        for s_name, dd in self.mdata["base_avg_retention_rate"].items():
            if "no_data_available" not in dd.keys():
                pdata_bybase[s_name] = dd

        pheader_byread = {
            "rca": {"color": "#D81B60", "name": "CpA Retention"},
            "rcc": {"color": "#1E88E5", "name": "CpC Retention"},
            "rcg": {"color": "#A0522D", "name": "CpG Retention"},
            "rct": {"color": "#004D40", "name": "CpT Retention"},
        }
        pheader_bybase = {
            "bca": {"color": "#D81B60", "name": "CpA Retention"},
            "bcc": {"color": "#1E88E5", "name": "CpC Retention"},
            "bcg": {"color": "#A0522D", "name": "CpG Retention"},
            "bct": {"color": "#004D40", "name": "CpT Retention"},
        }
        pconfig = {
            "id": "biscuit_retention",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "BISCUIT: Cytosine Retention",
            "data_labels": [{"name": "Read-averaged Retention"}, {"name": "Base-averaged Retention"}],
            "ylab": "Percent Retained",
            "ymin": 0,
            "ymax": 100,
            "yCeiling": 110,
            "stacking": None,
            "tt_decimals": 1,
            "tt_suffix": "%",
        }

        self.add_section(
            name="Cytosine Retention",
            anchor="biscuit-retention",
            description="Shows the cytosine retention rate for different contexts.",
            helptext="""
                The cytosine retention rate is calculated as `1 - (cytosine conversion rate)`.

                Assuming complete, but not over, bisulfite conversion, the cytosine retention rate
                is the average cytosine modification (including 5mC, 5hmC, etc) rate.

                Note, if a sample is missing from the Base-averaged Retention table,
                there wasn't sufficient data to plot that sample.
            """,
            plot=bargraph.plot([pdata_byread, pdata_bybase], [pheader_byread, pheader_bybase], pconfig),
        )

    def parse_logs_base_avg_retention_rate(self, f, fn):
        """
        Parses _totalBaseConversionRate.txt
        Inputs:
            f - current matched file
            fn - filename
        Returns:
            data - dictionary of base averaged fraction of retainied cytosines by context
        """

        file_data = f.splitlines()[2:]

        # Handle missing data
        if len(file_data) == 0:
            log.debug(f"No data available in {fn}. Will not fill corresponding entries.")
            return {"no_data_available": 1}

        data = {}
        for line in file_data:
            fields = line.split("\t")
            # Skip rows that have NaNs as something went wrong in processing
            if "nan" in fields:
                log.debug(f"Found NaN in {fn}. Skipping.")
                continue

            # BISCUIT returns -1 if insufficient data. Only fill fields with value >= 0.
            if float(fields[0]) >= 0:
                data["bca"] = 100.0 * float(fields[0])
            if float(fields[1]) >= 0:
                data["bcc"] = 100.0 * float(fields[1])
            if float(fields[2]) >= 0:
                data["bcg"] = 100.0 * float(fields[2])
            if float(fields[3]) >= 0:
                data["bct"] = 100.0 * float(fields[3])

        return data

    def chart_base_avg_retention_rate(self):
        """Handled by chart_read_avg_retention_rate()"""
        pass
