import csv
import logging
import random
from collections import defaultdict
from math import isinf, isnan

import spectra

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, heatmap, scatter, table
from multiqc.utils import mqc_colour

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Somalier",
            anchor="somalier",
            href="https://github.com/brentp/somalier",
            info="Genotype to pedigree correspondence checks from sketches derived from BAM/CRAM or VCF",
            extra="""
            Somalier can be used to find sample swaps or duplicates in cancer
            projects, where there is often no jointly-called VCF across samples.
        
            It is also extremely efficient and so can be used to find all-vs-all
            relatedness estimates for thousands of samples.
        
            It also outputs information on sex, depth, heterozgyosity, and ancestry
            to be used for general QC.
            """,
            doi="10.1186/s13073-020-00761-2",
        )

        # Find and load any somalier reports
        self.somalier_data = dict()
        self.somalier_background_pcs = dict()
        self.somalier_ancestry_cats = list()
        self.somalier_length_counts = dict()
        self.somalier_length_exp = dict()
        self.somalier_length_obsexp = dict()

        # parse somalier sample file
        for f in self.find_log_files("somalier/samples"):
            parsed_data = self.parse_somalier_samples(f)
            if parsed_data is not None:
                for s_name_raw in parsed_data:
                    s_name = "*".join([self.clean_s_name(s, f) for s in s_name_raw.split("*")])
                    if s_name in self.somalier_data.keys():
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name)
                    self.somalier_data[s_name] = parsed_data[s_name_raw]

        # parse somalier CSV files
        for f in self.find_log_files("somalier/pairs"):
            parsed_data = self.parse_somalier_pairs_tsv(f)
            if parsed_data is not None:
                for s_name_raw in parsed_data:
                    s_name = "*".join([self.clean_s_name(s, f) for s in s_name_raw.split("*")])
                    if s_name in self.somalier_data.keys():
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name)
                    self.somalier_data[s_name] = parsed_data[s_name_raw]

        # parse somalier ancestry files
        for f in self.find_log_files("somalier/somalier-ancestry", filehandles=True):
            self.parse_somalier_ancestry(f)

        # Filter to strip out ignored sample names
        self.somalier_data = self.ignore_samples(self.somalier_data)

        if len(self.somalier_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.somalier_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.somalier_data, "multiqc_somalier")

        # Somalier Stats Table
        self.somalier_stats_table()

        # Relatedness plot
        self.somalier_relatedness_plot()

        self.somalier_relatedness_heatmap_plot()

        # hetcheck plot
        self.somalier_het_check_plot()

        self.somalier_sex_check_plot()

        # ancestry plots
        self.somalier_ancestry_barplot()

        self.somalier_ancestry_pca_plot()

    @staticmethod
    def parse_somalier_samples(f):
        """Go through log file looking for somalier output"""
        parsed_data = dict()
        headers = None
        sample_i = -100
        for line in f["f"].splitlines():
            s = line.split("\t")
            if headers is None:
                s[0] = s[0].lstrip("#")
                headers = s
                sample_i = headers.index("sample_id")
            else:
                parsed_data[s[sample_i]] = dict()
                for i, v in enumerate(s):
                    if i != sample_i:
                        try:
                            parsed_data[s[sample_i]][headers[i]] = float(v)
                        except ValueError:
                            parsed_data[s[sample_i]][headers[i]] = v
        if len(parsed_data) == 0:
            return None
        return parsed_data

    @staticmethod
    def parse_somalier_pairs_tsv(f):
        """Parse csv output from somalier"""
        parsed_data = dict()
        headers = None
        s_name_idx = None
        for line in f["f"].splitlines():
            s = line.lstrip("#").split("\t")
            if headers is None:
                headers = s
                try:
                    s_name_idx = [headers.index("sample_a"), headers.index("sample_b")]
                except ValueError:
                    log.warning(f"Could not find sample name in somalier output: {f['fn']}")
                    return None
            else:
                s_name = "*".join([s[idx] for idx in s_name_idx])  # not safe to hard code, but works
                parsed_data[s_name] = dict()
                for i, v in enumerate(s):
                    if i not in s_name_idx:  # Skip if (i == 0 or 1); i.e. sample_a, sample_b
                        if isnan(float(v)) or isinf(float(v)):
                            # Inf or NaN indicate the absence of data
                            v = None
                        else:
                            # add the pattern as a suffix to key
                            v = float(v)
                        parsed_data[s_name][headers[i]] = v

        if len(parsed_data) == 0:
            return None
        return parsed_data

    def parse_somalier_ancestry(self, f):
        # dict for parsed data, ancestry prediction probabilities and PCs
        parsed_data = dict()

        # list for background principal components and associated ancestry
        bg_pc1 = []
        bg_pc2 = []
        bg_ancestry = []

        reader: csv.DictReader = csv.DictReader(f["f"], dialect="excel-tab")
        idx = "#sample_id"

        # check file not empty, else parse file
        if reader.fieldnames is not None:
            # get ancestry categories from header fieldnames
            for c in reader.fieldnames:
                # use _prob substring to identify ancestry cats and add to object
                if "_prob" in c:
                    self.somalier_ancestry_cats.append(c.replace("_prob", ""))

            # parse rows of tsv file
            for row in reader:
                # only background has given_ancestry
                # i.e. row is sample, not background
                if len(row["given_ancestry"]) == 0:
                    d = {}
                    for k, v in row.items():
                        # add ancestry prediction to d
                        if "_prob" in k and v is not None:
                            d[k.replace("_prob", "")] = float(v)
                        # add principal component to d
                        elif "PC" in k and v is not None:
                            d[k] = float(v)
                        # else: do nothing

                    # extract predicted ancestry and probability
                    # for general stats table
                    d["ancestry"] = row["predicted_ancestry"]
                    d["p_ancestry"] = d[d["ancestry"]]

                    parsed_data[row[idx]] = d
                else:  # row is background, parse PC's
                    bg_pc1.append(float(row["PC1"]))
                    bg_pc2.append(float(row["PC2"]))
                    bg_ancestry.append(row["given_ancestry"])
                    # background pc's are added to parsed_data later

            # check that something was parsed:
            if len(parsed_data) > 0:
                # add background principal components
                self.somalier_background_pcs["background_pcs"] = {"PC1": bg_pc1, "PC2": bg_pc2, "ancestry": bg_ancestry}

                # cycle over keys, i.e. sample names
                # safely add new data to object data
                # warn when overwriting
                for s_name_raw in parsed_data:
                    s_name = "*".join([self.clean_s_name(s, f) for s in s_name_raw.split("*")])
                    if s_name in self.somalier_data.keys():
                        intersect_keys = parsed_data[s_name_raw].keys() & self.somalier_data.keys()
                        if len(intersect_keys) > 0:
                            log.debug(f"Duplicate sample name found! Overwriting: {s_name} : {intersect_keys}")
                    self.add_data_source(f, s_name)
                    try:
                        self.somalier_data[s_name].update(parsed_data[s_name_raw])
                    except KeyError:
                        self.somalier_data[s_name] = parsed_data[s_name_raw]
        else:
            log.warning(f"Detected empty file: {f['fn']}")

    def somalier_stats_table(self):
        """Add data to somalier stats table

        Bigger table within the somalier module, showing more stats"""

        headers = {
            "phenotype": {
                "title": "Phenotype",
                "description": "Sample's phenotype from pedigree info",
                "hidden": True,
            },
            "original_pedigree_sex": {
                "title": "Sex",
                "description": "Sample's sex from pedigree info",
                "scale": False,
            },
            "paternal_id": {
                "title": "Father ID",
                "description": "ID of sample's father ",
                "scale": False,
                "hidden": True,
            },
            "maternal_id": {
                "title": "Mother ID",
                "description": "ID of sample's mother",
                "scale": False,
                "hidden": True,
            },
            "family_id": {
                "title": "Family ID",
                "description": "ID of sample's family",
                "scale": False,
                "hidden": True,
            },
            "sex": {
                "title": "Inferred sex",
                "description": "Sample's inferred sex",
                "scale": False,
                "hidden": True,
            },
            "ancestry": {"title": "Ancestry", "description": "Most probable ancestry background", "scale": False},
            "p_ancestry": {
                "title": "P(Ancestry)",
                "description": "Ancestry probability",
                "max": 1,
                "min": 0,
                "scale": "RdYlGn",
                "format": "{:,.2f}",
            },
            "n_het": {
                "title": "HetVar",
                "description": "Heterozygous variants",
                "shared_key": "variant_count",
                "format": "{:,.0f}",
            },
            "n_hom_ref": {
                "title": "HomRefVar",
                "description": "Homozygous reference variants",
                "shared_key": "variant_count",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "n_hom_alt": {
                "title": "HomAltVar",
                "description": "Homozygous alternate variants",
                "shared_key": "variant_count",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "n_unknown": {"title": "NA sites", "description": "Unknown sites", "format": "{:,.0f}"},
            "depth_mean": {
                "title": "Mean depth",
                "description": "Mean depth of all sites",
                "scale": "RdYlGn",
                "suffix": " X",
                "hidden": True,
            },
            "depth_sd": {
                "title": "Depth std",
                "description": "Depth's standard deviation of all sites",
                "scale": "RdYlGn",
                "hidden": True,
            },
            "gt_depth_mean": {
                "title": "Sites depth",
                "description": "Mean depth of genotyped sites",
                "scale": "RdYlGn",
                "suffix": " X",
            },
            "gt_depth_sd": {
                "title": "Genot depth std",
                "description": "Depth's standard deviation of genotype sites",
                "scale": "RdYlGn",
                "suffix": " X",
                "hidden": True,
            },
            "ab_mean": {
                "title": "Allele balance",
                "description": "Mean allele balance",
                "scale": "RdYlGn",
            },
            "ab_std": {
                "title": "Allele balance std",
                "description": "Standard deviation of allele balance",
                "scale": "RdYlGn",
                "hidden": True,
            },
            "p_middling_ab": {
                "title": "Allele balance < 0.2, > 0.8",
                "description": "Proportion of sites with allele balance < 0.2 or > 0.8",
                "max": 1,
                "min": 0,
                "scale": "RdYlGn",
                "format": "{:,.2f}",
            },
            "X_het": {
                "title": "HetVar X",
                "description": "Heterozygous variants on X chromosome",
                "shared_key": "variant_count_xy",
                "format": "{:,.0f}",
            },
            "X_hom_ref": {
                "title": "HomRefVar X",
                "description": "Homozygous reference variants on X chromosome",
                "shared_key": "variant_count_xy",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "X_hom_alt": {
                "title": "HomAltVar X",
                "description": "Homozygous alternate variants on X chromosome",
                "shared_key": "variant_count_xy",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "X_n": {
                "title": "Sites X",
                "description": "Total sites on X chromosome",
                "shared_key": "variant_count_xy",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "X_depth_mean": {
                "title": "Mean depth X",
                "description": "Mean depth of sites on X chromosome",
                "scale": "RdYlGn",
                "suffix": " X",
            },
            "Y_n": {
                "title": "Sites Y",
                "description": "Total sites on Y chromosome",
                "shared_key": "variant_count_xy",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "Y_depth_mean": {
                "title": "Mean depth Y",
                "description": "Mean depth of sites on Y chromosome",
                "scale": "RdYlGn",
                "suffix": " X",
            },
        }

        t_config = {
            "id": "somalier_stats",
            "namespace": "Somalier",
            "title": "Somalier: Statistics",
            "no_violin": True,
            "raw_data_fn": "multiqc_somalier_stats",
        }

        self.add_section(
            name="Statistics",
            anchor="somalier-stats",
            description="Various statistics from the somalier report.",
            plot=table.plot(self.somalier_data, headers, t_config),
        )

    def somalier_relatedness_plot(self):
        alpha = 0.6
        relatedness_groups = {
            0: {
                "name": "Unrelated",
                "color": f"rgba(74, 124, 182, {alpha})",
            },
            0.49: {
                "name": "Sib-sib",
                "color": f"rgba(243, 123, 40, {alpha})",
            },
            0.5: {
                "name": "Parent-child",
                "color": f"rgba(159, 84, 47, {alpha})",
            },
        }

        # Get index colour scale
        cscale = mqc_colour.mqc_colour_scale()
        extra_colours = cscale.get_colours("Dark2")
        extra_colours = _make_col_alpha(extra_colours, alpha)
        extra_colour_idx = 0
        data = dict()
        for pair, d in self.somalier_data.items():
            if "expected_relatedness" not in d:
                continue

            relatedness = d["expected_relatedness"]
            # -1 is not the same family, 0 is same family but unrelated
            # @brentp says he usually bundles them together
            if relatedness == -1:
                relatedness = 0

            # New unique value that we've not seen before
            if relatedness not in relatedness_groups:
                relatedness_groups[relatedness] = {
                    "name": str(relatedness),
                    "color": extra_colours[extra_colour_idx],
                }
                extra_colour_idx += 0
                if extra_colour_idx > len(extra_colours):
                    extra_colour_idx = 0

            data[pair] = {
                "x": d["ibs0"],
                "y": d["ibs2"],
                "color": relatedness_groups[relatedness]["color"],
                "group": relatedness_groups[relatedness]["name"],
            }

        if len(data) == 0:
            return

        pconfig = {
            "id": "somalier_relatedness_plot",
            "title": "Somalier: Sample Shared Allele Rates (IBS)",
            "xlab": "IBS0 (no alleles shared)",
            "ylab": "IBS2 (both alleles shared)",
            "marker_line_width": 0,
        }

        colours_legend = ""
        for rel, group in sorted(relatedness_groups.items()):
            col = group["color"].replace(str(alpha), "1.0")
            colours_legend += f'<span style="color:{col}">{group["name"]}</span>, '

        self.add_section(
            name="Relatedness",
            anchor="somalier-relatedness",
            description=f"""
            Shared allele rates between sample pairs.
            Points are coloured by degree of expected relatedness: {colours_legend}""",
            plot=scatter.plot(data, pconfig),
        )

    def somalier_relatedness_heatmap_plot(self):
        # inspiration: MultiQC/modules/vcftools/relatedness2.py

        data = []
        labels = set()
        rels = defaultdict(dict)
        for s_name, d in self.somalier_data.items():
            if "relatedness" in d:
                a, b = s_name.split("*")
                labels.add(a)
                labels.add(b)
                rels[a][b] = rels[b][a] = d["relatedness"]
                rels[a][a] = rels[b][b] = 1.0

        # impose alphabetical order and avoid json serialisation errors in utils.report
        labels = sorted(labels)

        for x in labels:
            line = []
            for y in labels:
                try:
                    line.append(rels[x][y])
                except KeyError:
                    line.append(None)
            data.append(line)

        if len(data) > 0:
            pconfig = {
                "id": "somalier_relatedness_heatmap_plot",
                "title": "Somalier: Sample Relatedness",
                "xlab": "Sample A",
                "ylab": "Sample B",
            }

            self.add_section(
                name="Relatedness Heatmap",
                anchor="somalier-relatedness-heatmap",
                description="Heatmap displaying relatedness of sample pairs.",
                plot=heatmap.plot(
                    data=data,
                    xcats=labels,
                    ycats=labels,
                    pconfig=pconfig,
                ),
            )

    def somalier_het_check_plot(self):
        """plot the het_check scatter plot"""
        # empty dictionary to add sample names, and dictionary of values
        data = {}

        # for each sample, and list in self.somalier_data
        for s_name, d in self.somalier_data.items():
            # check the sample contains the required columns
            if "gt_depth_mean" in d and "ab_std" in d:
                # add sample to dictionary with value as a dictionary of points to plot
                data[s_name] = {"x": d["gt_depth_mean"], "y": d["ab_std"]}

        if len(data) > 0:
            pconfig = {
                "id": "somalier_het_check_plot",
                "title": "Somalier: Sample Observed Heterozygosity",
                "xlab": "Mean depth",
                "xsuffix": "x",
                "ylab": "Standard deviation of allele-balance",
            }

            self.add_section(
                name="Heterozygosity",
                description="Standard deviation of heterozygous allele balance against mean depth.",
                helptext="A high standard deviation in allele balance suggests contamination.",
                anchor="somalier-hetcheck",
                plot=scatter.plot(data, pconfig),
            )

    def somalier_sex_check_plot(self):
        data = {}
        sex_index = {"female": 0, "male": 1, "unknown": 2}

        for s_name, d in self.somalier_data.items():
            if "X_depth_mean" in d and "original_pedigree_sex" in d:
                if d["gt_depth_mean"] == 0:
                    y = 0
                else:
                    y = 2 * d["X_depth_mean"] / d["gt_depth_mean"]
                data[s_name] = {
                    "x": sex_index.get(d["original_pedigree_sex"], 2) + (random.random() - 0.5) * 0.1,
                    "y": y,
                }

        if len(data) > 0:
            pconfig = {
                "id": "somalier_sex_check_plot",
                "title": "Somalier: Sample Predicted Sex",
                "xlab": "Sex from pedigree",
                "ylab": "Scaled mean depth on X",
                "categories": ["Female", "Male", "Unknown"],
            }

            self.add_section(
                name="Sex",
                description="Predicted sex against scaled depth on X",
                helptext="Higher values of depth, low values suggest male.",
                anchor="somalier-sexcheck",
                plot=scatter.plot(data, pconfig),
            )

    def somalier_ancestry_barplot(self):
        data = dict()
        c_scale = mqc_colour.mqc_colour_scale(name="Paired").colours
        cats = dict()
        anc_cats = self.somalier_ancestry_cats

        # use Paired color scale, unless number of categories exceed colors
        if len(anc_cats) <= len(c_scale):
            for i in range(len(anc_cats)):
                c = anc_cats[i]
                if i < (len(c_scale) - 1):
                    col = c_scale[i]
                else:
                    # default col if more cats than cols
                    col = "rgb(211,211,211,0.5)"
                cats[c] = {"name": c, "color": col}
        else:
            cats = None

        for s_name, d in self.somalier_data.items():
            # ensure that only relevant items are added,
            # i.e. only ancestry category values
            ls = {k: v * 100.0 for k, v in d.items() if (k in self.somalier_ancestry_cats)}
            if len(ls) > 0:  # only add dict, if it contains values
                data[s_name] = ls

        if len(data) > 0:
            pconfig = {
                "id": "somalier_ancestry_barplot",
                "title": "Somalier: Sample Predicted Ancestry Proportions",
                "cpswitch": False,
                "hide_empty": False,
                "ylab": "Predicted Ancestry",
                "tt_suffix": "%",
            }

            self.add_section(
                name="Ancestry Barplot",
                description="Predicted ancestries of samples.",
                helptext="""
                Shows the percentwise predicted probability of each
                ancestry. A sample might contain traces of several ancestries.
                If the number of samples is too high, the plot is rendered as a
                non-interactive flat image.
                """,
                anchor="somalier-ancestry",
                plot=bargraph.plot(data=data, cats=cats, pconfig=pconfig),
            )

    def somalier_ancestry_pca_plot(self):
        data = dict()

        # add background
        # N.B. this must be done after samples to have samples on top
        d = self.somalier_background_pcs.pop("background_pcs", {})
        if d:
            # generate color scale to match the number of categories
            c_scale = mqc_colour.mqc_colour_scale(name="Paired").colours
            cats = self.somalier_ancestry_cats
            ancestry_colors = dict(zip(cats, c_scale[: len(cats)]))
            default_background_color = "rgba(255,192,203,0.3)"

            # Make colours semi-transparent
            ancestry_colors = dict(zip(ancestry_colors.keys(), _make_col_alpha(ancestry_colors.values(), 0.3)))

            background = [
                {
                    "x": pc1,
                    "y": pc2,
                    "color": ancestry_colors.get(ancestry, default_background_color),
                    "name": ancestry,
                    "marker_size": 3,
                    "marker_line_width": 0,
                    "annotate": False,
                }
                for pc1, pc2, ancestry in zip(d["PC1"], d["PC2"], d["ancestry"])
            ]
            data["background"] = background

        # cycle over samples and add PC coordinates to data dict
        for s_name, d in self.somalier_data.items():
            if "PC1" in d and "PC2" in d:
                data[s_name] = {
                    "x": d["PC1"],
                    "y": d["PC2"],
                    "color": "rgba(0, 0, 0, 0.6)",
                }

        # generate section and plot
        if len(data) > 0:
            pconfig = {
                "id": "somalier_ancestry_pca_plot",
                "title": "Somalier: Sample Predicted Ancestry",
                "xlab": "PC1",
                "ylab": "PC2",
                "marker_size": 5,
            }

            self.add_section(
                name="Ancestry PCA",
                description="Principal components of samples against background PCs.",
                helptext="""
                Sample PCs are plotted against background PCs from the
                background data supplied to somalier.
                Color indicates predicted ancestry of sample. Data points in close
                proximity are predicted to be of similar ancestry. Consider whether
                the samples cluster as expected.
                """,
                anchor="somalier-ancestry-pca",
                plot=scatter.plot(data, pconfig),
            )


def _make_col_alpha(cols, alpha):
    """Take an HTML colour value and return a rgba string with alpha"""
    cols_return = []
    for col in cols:
        col_srgb = spectra.html(col)
        cols_rgb = [c * 255.0 for c in col_srgb.clamped_rgb]
        cols_return.append("rgba({},{},{},{})".format(*cols_rgb, alpha))
    return cols_return
