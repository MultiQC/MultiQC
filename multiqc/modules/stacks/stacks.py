""" MultiQC module to parse Stacks 2 denovo output"""


import logging
import os
import re

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Stacks",
            anchor="stacks",
            href="http://catchenlab.life.illinois.edu/stacks/",
            info="A software for analyzing restriction enzyme-based data (e.g. RAD-seq).",
            doi="10.1111/mec.12354",
        )

        gsheaders = {
            "n_loci": {
                "title": "# loci",
                "description": "Number of loci built",
                "format": "{:,.i}",
                "scale": "RdYlGn",
            },
            "n_used_fw_reads": {
                "title": "K reads used",
                "modify": lambda x: float(x) / 1000.0,
                "description": "Number of thousand reads used",
                "scale": "BuGn",
            },
            "mean_cov": {
                "title": "cov",
                "suffix": "X",
                "description": "Mean sequence coverage at locus",
                "scale": "BuPu",
            },
            "mean_cov_ns": {
                "title": "weighted cov",
                "suffix": "X",
                "description": "The coverage at each locus is weighted by the number of samples present at that locus (i.e. coverage at shared loci counts more)",
                "scale": "YlGn",
            },
        }

        sheaders = {
            "# Pop ID": {
                "title": "PopID",
                "description": "Population ID as defined in the Population Map file.",
                "scale": False,
                "format": "{:,.s}",
            },
            "Private": {
                "title": "Private",
                "description": "Number of private alleles in this population.",
                "scale": "PuBu",
                "hidden": True,
            },
            "Num_Indv": {
                "title": "# Indv",
                "description": "Mean number of individuals per locus in this population.",
                "scale": "YlGn",
            },
            "P": {
                "title": "P",
                "description": "Mean frequency of the most frequent allele at each locus in this population.",
                "scale": "PuBu",
                "min": 0,
                "max": 1,
            },
            "Obs_Het": {
                "title": "Obs Het",
                "description": "Mean observed heterozygosity in this population.",
                "scale": "YlGn",
                "min": 0,
                "max": 1,
            },
            "Obs_Hom": {
                "title": "Obs Hom",
                "description": "Mean observed homozygosity in this population.",
                "scale": "PuBu",
                "min": 0,
                "max": 1,
                "hidden": True,
            },
            "Exp_Hom": {
                "title": "Exp_Hom",
                "description": "Mean expected homozygosity in this population.",
                "scale": "YlGn",
                "min": 0,
                "max": 1,
                "hidden": True,
            },
            "Exp_Het": {
                "title": "Exp Het",
                "description": "Mean expected heterozygosity in this population.",
                "scale": "PuBu",
                "min": 0,
                "max": 1,
            },
            "Pi": {
                "title": "Pi",
                "description": "Mean value of &#960; in this population.",
                "scale": "YlGn",
                "min": 0,
                "max": 1,
            },
            "Fis": {
                "title": "Fis",
                "description": "Mean measure of Fis in this population.",
                "scale": "PuOr",
                "min": -1,
                "max": 1,
            },
        }

        num_files = 0
        # Parse gstacks data
        cov_data = {}
        for f in self.find_log_files("stacks/gstacks"):
            run_name = os.path.dirname(f["root"])
            s_name = self.clean_s_name(os.path.basename(f["root"]), f, root=run_name)
            try:
                cov_data.update(self.parse_gstacks(f["f"], s_name))
                self.add_data_source(f, section="gstacks")
                num_files += 1
            except Exception as e:
                log.error(f"Could not parse gstacks.distribs file in {f['s_name']}: {e}")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Parse populations data
        distribs_loci = {}
        distribs_snps = {}
        for f in self.find_log_files("stacks/populations"):
            run_name = os.path.dirname(f["root"])
            s_name = self.clean_s_name(os.path.basename(f["root"]), f, root=run_name)
            i, j = self.parse_populations(f["f"], s_name)
            try:
                distribs_loci.update(i)
                distribs_snps.update(j)
                self.add_data_source(f, section="populations")
                num_files += 1
            except Exception as e:
                log.error(f"Could not parse population.log.distribs file in {f['s_name']}: {e}")

        # Parse sumstats file
        sumstats_data = {}
        for f in self.find_log_files("stacks/sumstats"):
            run_name = os.path.dirname(f["root"])
            s_name = self.clean_s_name(os.path.basename(f["root"]), f, root=run_name)
            try:
                sumstats_data.update(self.parse_sumstats(f["f"], s_name))
                self.add_data_source(f, section="sumstats")
                num_files += 1
            except Exception as e:
                log.error(f"Could not parse populations.sumstats_summary file in {f['s_name']}: {e}")

        # Ignore samples
        cov_data = self.ignore_samples(cov_data)
        distribs_loci = self.ignore_samples(distribs_loci)
        distribs_snps = self.ignore_samples(distribs_snps)
        sumstats_data = self.ignore_samples(sumstats_data)

        if len(cov_data) == 0 and len(sumstats_data) == 0 and len(distribs_loci) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {num_files} reports")

        # Write parsed report data to a file
        self.write_data_file(cov_data, "multiqc_stacks_cov")
        self.write_data_file(sumstats_data, "multiqc_stacks_sumstats")

        # Write the sample table
        config_table = {
            "id": "gstacks_table",
            "namespace": "stacks",
            "title": "Stacks: Sample statistics",
        }
        self.add_section(
            name="Sample statistics",
            anchor="stacks-gstacks",
            description="The sample specific statistics for Stacks",
            helptext="""**Note!** The sample names have the following scheme `<run folder name> | <input fastq file prefix>`.
                        This data is obtained from the gstacks program run after builing sample and catalog loci merge
                        paired-ends and call variants.
                        These numbers are obtained from the `gstacks.log.distribs` file""",
            plot=table.plot(cov_data, gsheaders, config_table),
        )
        # Write population sumstats table
        config_table = {
            "id": "sumstats_table",
            "namespace": "stacks",
            "title": "Stacks: Population summary statistics",
        }
        self.add_section(
            name="Population summary statistics",
            anchor="stacks-sumstats",
            description="Population statistics as calculated from variant sites found in this run",
            helptext="""**Note!** The sample names have the following scheme `<run folder name> | <population ID>`,
                        where the population ID is defined in the input population map file.
                        This information is obtained from the Stacks program `population` and the file populations.sumstats_summary.tsv
                        """,
            plot=table.plot(sumstats_data, sheaders, config_table),
        )
        config_distribs = {
            "id": "distribs_plot",
            "title": "Stacks: Population plots",
            "namespace": "stacks",
            "tt_label": "{point.y} loci, {point.x} samples/SNPs",
            "ylab": "# loci",
            "data_labels": [
                {"name": "Samples per loci", "ylab": "# loci", "xlab": "# samples"},
                {"name": "SNPs per loci", "ylab": "# loci", "xlab": "# SNPs"},
            ],
        }
        self.add_section(
            name="Population plots",
            anchor="stacks-distribs",
            description="Plots showing, 1) the number of loci shared by number of samples and 2) the number of SNPs per sample",
            helptext="""The distributions are obtained from the Stacks program `populations` and it's output file `populations.log.distribs`.
            These numbers are Stacks' post-filtering.""",
            plot=linegraph.plot([distribs_loci, distribs_snps], config_distribs),
        )

    @staticmethod
    def parse_gstacks(file_contents, s_name):
        headers = None
        out_dict = {}

        for line in file_contents.splitlines():
            if line.startswith("sample	n_loci	n_used_fw_reads"):
                headers = ["n_loci", "n_used_fw_reads", "mean_cov", "mean_cov_ns"]
            elif line.startswith("END effective_coverages_per_sample"):
                break
            elif headers is not None:
                cdict = dict.fromkeys(headers)
                content = line.split("\t")
                for i in range(1, len(content)):
                    cdict[list(cdict.keys())[i - 1]] = content[i]
                if s_name is not None:
                    out_dict[s_name + " | " + content[0]] = cdict
                else:
                    out_dict[content[0]] = cdict

        return out_dict

    @staticmethod
    def parse_sumstats(file_contents, s_name):
        out_dict = dict()
        # ["# Pop ID","Private","Num_Indv","P","Obs_Het","Obs_Hom","Exp_Het","Exp_Hom","Pi","Fis"]
        fields = [0, 1, 2, 5, 8, 11, 14, 17, 20, 23]
        fl = file_contents.splitlines()
        var_lines = fl[fl.index("# Variant positions") + 1 : fl.index("# All positions (variant and fixed)")]
        headers = None
        for line in var_lines:
            if line.startswith("# Pop ID	Private	Num_Indv"):
                headers = line.split("\t")
            elif headers is not None:
                cdict = {}
                content = line.split("\t")
                for i in fields:
                    cdict[headers[i]] = content[i]
                out_dict[s_name + " | " + content[0]] = cdict

        return out_dict

    @staticmethod
    def parse_populations(file_contents, s_name):
        loci_dict = {}
        snps_dict = {}
        pat = re.compile("BEGIN (.*)\n\n{0,1}#.*\n.+\n((?:.+\n)+)END", re.MULTILINE)
        for match in pat.finditer(file_contents):
            title, table = match.groups()
            title = title.strip()
            if title == "missing_samples_per_loc_postfilters":
                if s_name is not None:
                    title = s_name
                loci_dict[title] = {}
                for line in table.splitlines():
                    row = line.split()
                    loci_dict[title][int(row[0])] = int(row[1])
            elif title == "snps_per_loc_postfilters":
                if s_name is not None:
                    title = s_name
                snps_dict[title] = {}
                for line in table.splitlines():
                    row = line.split()
                    snps_dict[title][int(row[0])] = int(row[1])

        return loci_dict, snps_dict
