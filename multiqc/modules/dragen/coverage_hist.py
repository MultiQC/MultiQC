import logging
import re
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule
from multiqc.modules.qualimap.QM_BamQC import coverage_histogram_helptext, genome_fraction_helptext
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class DragenCoverageHist(BaseMultiqcModule):
    def add_coverage_hist(self):
        data_by_phenotype_by_sample = defaultdict(dict)
        for f in self.find_log_files("dragen/wgs_fine_hist"):
            data_by_phenotype = parse_wgs_fine_hist(f)
            s_name = f["s_name"]
            if s_name in data_by_phenotype_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="wgs_fine_hist")
            data_by_phenotype_by_sample[s_name].update(data_by_phenotype)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

        # Filter to strip out ignored sample names:
        data_by_phenotype_by_sample = self.ignore_samples(data_by_phenotype_by_sample)

        # Merge tumor and normal data:
        data_by_sample = defaultdict(dict)
        for sn in data_by_phenotype_by_sample:
            for phenotype in data_by_phenotype_by_sample[sn]:
                new_sn = sn
                if phenotype == "normal":
                    new_sn = sn + "_normal"
                data_by_sample[new_sn] = data_by_phenotype_by_sample[sn][phenotype]
        if not data_by_sample:
            return set()

        # Only plot data, don't want to write this to a file
        # (can do so with --export-plots already)
        # self.write_data_file(data_by_sample, "dragen_cov_hist")

        dist_data = {sn: dist for sn, (dist, cum, depth_1pc) in data_by_sample.items()}
        cum_data = {sn: cum for sn, (dist, cum, depth_1pc) in data_by_sample.items()}
        depth_1pc = max(depth_1pc for (dist, cum, depth_1pc) in data_by_sample.values())

        self.add_section(
            name="Coverage distribution",
            anchor="dragen-coverage-distribution",
            description="Number of locations in the reference genome with a given depth of coverage.",
            helptext=coverage_histogram_helptext,
            plot=linegraph.plot(
                dist_data,
                {
                    "id": "dragen_coverage_dist",
                    "title": "Dragen: Coverage distribution",
                    "xlab": "Depth (x)",
                    "ylab": "Number of bases in genome covered by X reads",
                    "ymin": 0,
                    "xmin": 0,
                    "xmax": depth_1pc,  # trim long flat tail
                    "tt_label": "<b>{point.x}</b>: {point.y} loci",
                    "cpswitch": True,
                },
            ),
        )

        self.add_section(
            name="Cumulative coverage hist",
            anchor="dragen-cum-coverage-histogram",
            description="Number of locations in the reference genome with at least given depth of coverage.",
            helptext=genome_fraction_helptext,
            plot=linegraph.plot(
                cum_data,
                {
                    "id": "dragen_cumulative_coverage_hist",
                    "title": "Dragen: Cumulative coverage hist",
                    "xlab": "Depth (x)",
                    "ylab": "% of bases in genome covered by at least X reads",
                    "ymin": 0,
                    "ymax": 100,
                    "xmin": 0,
                    "xmax": depth_1pc,  # trim long flat tail
                    "tt_label": "<b>{point.x}</b>: {point.y:.2f}%",
                },
            ),
        )
        return data_by_sample.keys()


def parse_wgs_fine_hist(f):
    """
    T_SRR7890936_50pc.wgs_fine_hist_normal.csv
    T_SRR7890936_50pc.wgs_fine_hist_tumor.csv

    Depth,Overall
    0,104231614
    1,9430586
    2,5546235
    ...
    998,208
    999,177
    1000+,201801

    Contains two columns: Depth and Overall. The value in the Depth column ranges from 0 to 1000+
    and the Overall column indicates the number of loci covered at the corresponding depth.

    Parsing all values except for 1000+ and plotting a distribution histogram and cumulative histogram
    """

    # first pass to calculate total number of bases to calculate percentages
    parsed_data = dict()
    for line in f["f"].splitlines():
        if line.startswith("Depth,Overall"):
            continue
        key, cnt = line.split(",")
        try:
            cnt = int(cnt)
        except ValueError:
            continue
        parsed_data[key] = cnt

    total_cnt = sum(parsed_data.values())

    data = dict()
    cum_data = dict()
    cum_cnt = 0
    depth_1pc = None

    for key, cnt in reversed(list(parsed_data.items())):
        try:
            depth = int(key)
        except ValueError:
            continue
        cum_cnt += cnt
        if total_cnt > 0:
            cum_pct = cum_cnt / total_cnt * 100.0
        else:
            cum_pct = 0
        if cum_pct < 1:  # to trim long flat tail
            depth_1pc = depth
        data[depth] = cnt
        cum_data[depth] = cum_pct

    m = re.search(r"(tumor|normal).csv", f["fn"])
    if m:
        phenotype = m.group(1)
    else:
        phenotype = "unknown"
    return {phenotype: (data, cum_data, depth_1pc)}
