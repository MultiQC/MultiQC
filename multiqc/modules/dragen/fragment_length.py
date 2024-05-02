# Initialise the logger
import logging
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


MIN_CNT_TO_SHOW_ON_PLOT = 5


class DragenFragmentLength(BaseMultiqcModule):
    def add_fragment_length_hist(self):
        data_by_rg_by_sample = defaultdict(dict)

        for f in self.find_log_files("dragen/fragment_length_hist"):
            data_by_rg = parse_fragment_length_hist_file(f)
            s_name = f["s_name"]
            if s_name in data_by_rg_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="fragment_length_hist")

            for rg, data in data_by_rg.items():
                if any(rg in d_rg for sn, d_rg in data_by_rg_by_sample.items()):
                    log.debug(f"Duplicate read group name {rg} found for {s_name}! Overwriting")
            data_by_rg_by_sample[s_name].update(data_by_rg)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

        # Filter to strip out ignored sample names:
        data_by_rg_by_sample = self.ignore_samples(data_by_rg_by_sample)
        if not data_by_rg_by_sample:
            return set()

        # Write data to file
        self.write_data_file(data_by_rg_by_sample, "dragen_frag_len")

        # Merging all data
        data_by_rg = {}
        for sn, d_rg in data_by_rg_by_sample.items():
            for rg, d in d_rg.items():
                if rg in d_rg:
                    rg = rg + " (" + sn + ")"
                data_by_rg[rg] = d

        # Exit early if we have no valid data, such as from FastQcOnly runs
        if not data_by_rg:
            return set()

        smooth_points = 300
        self.add_section(
            name="Fragment length hist",
            anchor="dragen-fragment-length-histogram",
            description="""
            Distribution of estimated fragment lengths of mapped reads per read group.
            Only points supported by at least {} reads are shown to prevent long flat tail.
            The plot is also smoothed down to showing 300 points on the X axis to reduce noise.
            """.format(MIN_CNT_TO_SHOW_ON_PLOT),
            plot=linegraph.plot(
                data_by_rg,
                {
                    "id": "dragen_fragment_length",
                    "title": "Dragen: Fragment length hist",
                    "ylab": "Number of reads",
                    "xlab": "Fragment length (bp)",
                    "ymin": 0,
                    "xmin": 0,
                    "tt_label": "<b>{point.x} bp</b>: {point.y} reads",
                    "smooth_points": smooth_points,
                },
            ),
        )
        return data_by_rg_by_sample.keys()


def parse_fragment_length_hist_file(f):
    """
    T_SRR7890936_50pc.fragment_length_hist.csv

    #Sample: N_SRR7890889
    FragmentLength,Count
    36,1
    37,0
    38,0
    39,0
    40,0
    41,1
    ...
    39203,0
    39204,0
    39205,1
    #Sample: T_SRR7890936_50pc
    FragmentLength,Count
    53,2
    54,0
    ...
    39316,0
    39317,1
    """

    data_by_rg = defaultdict(dict)

    read_group = None
    for line in f["f"].splitlines():
        if line.startswith("#Sample"):
            read_group = line.split("#Sample: ")[1]
        else:
            assert read_group is not None
            frag_len, cnt = line.split(",")
            try:
                frag_len = int(frag_len)
                cnt = int(cnt)
            except ValueError:
                assert line == "FragmentLength,Count", line
            else:
                if cnt >= MIN_CNT_TO_SHOW_ON_PLOT:  # to prevent long flat tail
                    data_by_rg[read_group][frag_len] = cnt

    return data_by_rg
