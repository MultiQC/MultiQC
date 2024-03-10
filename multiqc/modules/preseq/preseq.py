""" MultiQC module to parse output from Preseq """


import logging

import numpy as np

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph
from multiqc.utils import mqc_colour

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Preseq",
            anchor="preseq",
            href="http://smithlabresearch.org/software/preseq/",
            info="""estimates the complexity of a library, showing how many additional
                    unique reads are sequenced for increasing total read count.
                    A shallow curve indicates complexity saturation. The dashed line
                    shows a perfectly complex library where total reads = unique reads.""",
            doi="10.1038/nmeth.2375",
        )

        # Find and load any Preseq reports
        reads_data = dict()
        bases_data = dict()
        for f in self.find_log_files("preseq"):
            sample_data_raw, sample_data_is_bases = _parse_preseq_logs(f)
            if sample_data_raw is None:
                continue

            if f["s_name"] in sample_data_raw:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")

            if sample_data_is_bases:
                bases_data[f["s_name"]] = sample_data_raw
            else:
                reads_data[f["s_name"]] = sample_data_raw

            self.add_data_source(f)

        # Filter to strip out ignored sample names
        bases_data = self.ignore_samples(bases_data)
        reads_data = self.ignore_samples(reads_data)
        all_data = {**bases_data, **reads_data}
        if not all_data:
            raise ModuleNoSamplesFound

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        if bases_data and reads_data:
            log.warning("Mixed 'TOTAL_READS' and 'TOTAL_BASES' reports. Will build two separate plots")
        log.info(f"Found {len(all_data)} reports")
        # Write data to file
        self.write_data_file(all_data, "preseq")

        # Preseq plot. If mixed logs are found, build two separate plots.
        if bases_data:
            self._make_preseq_length_trimmed_plot(bases_data, True)
        if reads_data:
            self._make_preseq_length_trimmed_plot(reads_data, False)

    def _make_preseq_length_trimmed_plot(self, data_raw, is_basepairs):
        """Generate the preseq plot.

        For Y axis, plot coverages if `config.preseq.read_length` and `config.preseq.genome_size`
        specified (unless requested otherwise in `config.preseq.y_axis`).

        For X axis, plot counts (or base pairs) in X axis (unless coverages requested
        explicitly in `config.preseq.x_axis`).
        """
        counts_in_1x = _get_counts_in_1x(is_basepairs)
        x_axis = getattr(config, "preseq", {}).get("x_axis", "counts")
        y_axis = getattr(config, "preseq", {}).get("y_axis", "coverage" if counts_in_1x else "counts")

        # Modify counts
        d_cnts = {sn: _modify_raw_data(sample_data, is_basepairs) for sn, sample_data in data_raw.items()}

        # Convert counts (base pairs) -> depths
        d_covs = {
            sn: _counts_to_coverages(sample_data, counts_in_1x) if counts_in_1x else None
            for sn, sample_data in data_raw.items()
        }

        # Prepare final dataset for plotting
        data = dict()
        for sn, s_d_cnts, s_d_covs in zip(d_cnts, d_cnts.values(), d_covs.values()):
            keys = s_d_covs.keys() if x_axis == "coverage" else s_d_cnts.keys()
            values = s_d_covs.values() if y_axis == "coverage" else s_d_cnts.values()
            data[sn] = dict(zip(keys, values))

        # Count maximum values to draw the "ideal" line
        max_y_raw, max_sn = max((max(sd.values()), sn) for sn, sd in data_raw.items())
        max_y_cnt = list(_modify_raw_data({max_y_raw: max_y_raw}, is_basepairs).items())[0][0]
        max_y_cov = list(_counts_to_coverages({max_y_raw: max_y_raw}, counts_in_1x).items())[0][0]
        max_y = max_y_cov if y_axis == "coverage" else max_y_cnt
        max_yx = max_y_cov if x_axis == "coverage" else max_y_cnt

        # Preparing axis and tooltip labels
        x_suffix, y_tt_lbl, x_axis_name, y_suffix, x_tt_lbl, y_axis_name = _prepare_labels(
            is_basepairs, max_y_cov, x_axis, y_axis
        )

        name = "Complexity curve"
        description = ""
        section_id = "preseq_plot"
        pconfig = {
            "id": "preseq_complexity_plot",
            "title": "Preseq: Complexity curve",
            "xlab": x_axis_name,
            "ylab": y_axis_name,
            "xmin": 0,
            "ymin": 0,
            "tt_label": "<b>" + y_tt_lbl + "</b>: " + x_tt_lbl,
            "xsuffix": x_suffix,
            "ysuffix": y_suffix,
            "extra_series": [],
        }
        if not is_basepairs:
            pconfig["title"] += " (molecule count)"
            pconfig["id"] += "_molecules"
            name += " (molecule count)"
            section_id += "_molecules"

        # Parse the real counts if we have them
        real_cnts_all, real_cnts_unq = self._parse_real_counts(data.keys())
        real_vals_all, real_vals_unq = _prep_real_counts(
            real_cnts_all, real_cnts_unq, is_basepairs, counts_in_1x, x_axis, y_axis
        )
        pconfig["extra_series"].extend(
            _real_counts_to_plot_series(data, real_vals_unq, real_vals_all, x_suffix, y_suffix, y_tt_lbl)
        )
        if real_vals_unq:
            description += "<p>Points show read count versus deduplicated read counts (externally calculated).</p>"
        elif real_vals_all:
            description += "<p>Points show externally calculated read counts on the curves.</p>"

        # Trim the data to not have a ridiculous x-axis (10Gbp anyone?)
        if getattr(config, "preseq", {}).get("notrim", False) is not True:
            max_y *= 0.8
            max_yx *= 0.8
            max_x = 0
            for x in sorted(list(data[max_sn].keys())):
                max_x = max(max_x, x)
                if data[max_sn][x] > max_y and x > real_vals_all.get(max_sn, 0) and x > real_vals_unq.get(max_sn, 0):
                    break
            pconfig["xmax"] = max_x
            description += "<p>Note that the x-axis is trimmed at the point where all the datasets \
                show 80% of their maximum y-value, to avoid ridiculous scales.</p>"

        # Plot perfect library as dashed line
        pconfig["extra_series"].append(
            {
                "name": "A perfect library where each read is unique",
                "data": [[0, 0], [max_yx, max_y]],
                "dash": "dash",
                "line": {"width": 1},
                "color": "#000000",
                "showlegend": False,
            }
        )

        self.add_section(name=name, description=description, anchor=section_id, plot=linegraph.plot(data, pconfig))

    def _parse_real_counts(self, sample_names):
        real_counts_file_raw = None
        real_counts_file_name = None
        for f in self.find_log_files("preseq/real_counts"):
            if real_counts_file_raw is not None:
                log.warning(f"Multiple Preseq real counts files found, now using {f['fn']}")
            real_counts_file_raw = f["f"]
            real_counts_file_name = f["fn"]

        real_counts_total = {}
        real_counts_unique = {}

        if real_counts_file_raw is not None:
            try:
                for line in real_counts_file_raw.splitlines():
                    if not line.startswith("#"):
                        cols = line.strip().split()  # Split on any whitespace
                        sn = self.clean_s_name(cols[0], None)
                        if sn in sample_names:
                            if len(cols) >= 2:
                                if cols[1].isdigit():
                                    real_counts_total[sn] = int(cols[1])
                            if len(cols) >= 3:
                                if cols[2].isdigit():
                                    real_counts_unique[sn] = int(cols[2])
            except IOError as e:
                log.error(f"Error loading real counts file {real_counts_file_name}: {str(e)}")
            else:
                log.debug(f"Found {len(real_counts_total)} matching sets of counts from {real_counts_file_name}")

        return real_counts_total, real_counts_unique


def _parse_preseq_logs(f):
    """Go through log file looking for preseq output"""

    lines = f["f"].splitlines()
    header = lines.pop(0)

    data_is_bases = False
    if header.startswith("TOTAL_READS	EXPECTED_DISTINCT"):
        pass
    elif header.startswith("TOTAL_BASES	EXPECTED_DISTINCT"):
        data_is_bases = True
    elif header.startswith("total_reads	distinct_reads"):
        pass
    else:
        log.debug(f"First line of preseq file {f['fn']} did not look right")
        return None, None

    data = dict()
    for line in lines:
        s = line.split()
        # Sometimes the Expected_distinct count drops to 0, not helpful
        if float(s[1]) == 0 and float(s[0]) > 0:
            continue
        data[float(s[0])] = float(s[1])

    return data, data_is_bases


def _modify_raw_data(sample_data, is_basepairs):
    """Modify counts or base pairs according to `read_count_multiplier`
    or `base_count_multiplier`.
    """
    return {_modify_raw_val(x, is_basepairs): _modify_raw_val(y, is_basepairs) for x, y in sample_data.items()}


def _modify_raw_val(val, is_basepairs):
    """Modify counts or base pairs according to `read_count_multiplier`
    or `base_count_multiplier`.
    """
    return float(val) * (config.base_count_multiplier if is_basepairs else config.read_count_multiplier)


def _counts_to_coverages(sample_data, counts_in_1x):
    """If the user specified read length and genome size in the config,
    convert the raw counts/bases into the depth of coverage.
    """
    if not counts_in_1x:
        return {None: None}

    return {_count_to_coverage(x, counts_in_1x): _count_to_coverage(y, counts_in_1x) for x, y in sample_data.items()}


def _count_to_coverage(val, counts_in_1x):
    return val / counts_in_1x


def _get_counts_in_1x(data_is_basepairs):
    """Read length and genome size from the config and calculate
    the approximate number of counts (or base pairs) in 1x of depth
    """
    read_length = float(getattr(config, "preseq", {}).get("read_length", 0))
    genome_size = getattr(config, "preseq", {}).get("genome_size")
    if genome_size:
        try:
            genome_size = float(genome_size)
        except ValueError:
            presets = {"hg19_genome": 2897310462, "hg38_genome": 3049315783, "mm10_genome": 2652783500}
            if genome_size in presets:
                genome_size = presets[genome_size]
            else:
                log.warning(
                    "The size for genome "
                    + genome_size
                    + " is unknown to MultiQC, "
                    + "please specify it explicitly or choose one of the following: "
                    + ", ".join(presets.keys())
                    + ". Falling back to molecule counts."
                )
                genome_size = None

    if genome_size:
        if data_is_basepairs:
            return genome_size
        elif read_length:
            return genome_size / read_length
    else:
        return None


def _prepare_labels(is_basepairs, max_y_cov, x_axis, y_axis):
    cov_suffix = "x"

    cov_lbl = None
    if x_axis == "coverage" or y_axis == "coverage":
        cov_precision = "2"
        if max_y_cov > 30:  # no need to be so precise when the depth numbers are high
            cov_precision = "1"
        if max_y_cov > 300:  # when the depth are very high, decimal digits are excessive
            cov_precision = "0"
        cov_lbl = "{value:,." + cov_precision + "f}x"

    cnt_lbl = "{value:,.2f} " + config.read_count_prefix
    cnt_suffix = config.read_count_prefix
    if config.read_count_multiplier == 1:
        cnt_lbl = "{value:,.0f}"
    if is_basepairs:
        cnt_lbl = "{value:,.2f} " + config.base_count_prefix
        cnt_suffix = " " + config.base_count_prefix
        if config.base_count_multiplier == 1:
            cnt_lbl = "{value:,.0f}"

    if x_axis == "coverage":
        x_suffix = cov_suffix
        x_tt_lbl = cov_lbl.replace("value", "point.x")
        x_axis_name = "Total coverage (including duplicates)"
    else:
        x_suffix = cnt_suffix
        x_tt_lbl = cnt_lbl.replace("value", "point.x")
        if is_basepairs:
            x_tt_lbl += " pairs (total)"
            x_axis_name = "Base pairs (total)"
        else:
            x_tt_lbl += " total molecules"
            x_axis_name = "Total molecules (including duplicates)"

    if y_axis == "coverage":
        y_suffix = cov_suffix
        y_tt_lbl = cov_lbl.replace("value", "point.y") + " depth"
        y_axis_name = "Unique coverage"
    else:
        y_suffix = cnt_suffix
        y_tt_lbl = cnt_lbl.replace("value", "point.y")
        if is_basepairs:
            y_tt_lbl += " pairs (unique)"
            y_axis_name = "Base pairs (unique reads)"
        else:
            y_tt_lbl += " unique molecules"
            y_axis_name = "Unique molecules"
    return x_suffix, y_tt_lbl, x_axis_name, y_suffix, x_tt_lbl, y_axis_name


def _prep_real_counts(real_cnts_all, real_cnts_unq, is_basepairs, counts_in_1x, x_axis, y_axis):
    if x_axis == "coverage":
        real_vals_all = {sn: _count_to_coverage(x, counts_in_1x) for sn, x in real_cnts_all.items()}
    else:
        real_vals_all = {sn: _modify_raw_val(x, is_basepairs) for sn, x in real_cnts_all.items()}

    if y_axis == "coverage":
        real_vals_unq = {sn: _count_to_coverage(y, counts_in_1x) for sn, y in real_cnts_unq.items()}
    else:
        real_vals_unq = {sn: _modify_raw_val(y, is_basepairs) for sn, y in real_cnts_unq.items()}

    return real_vals_all, real_vals_unq


def _real_counts_to_plot_series(data, yx_by_sample, xs_by_sample, x_suffix, y_suffix, y_tt_lbl):
    scale = mqc_colour.mqc_colour_scale("plot_defaults")

    series = []
    for si, sn in enumerate(sorted(data.keys())):
        if sn in xs_by_sample:
            x = float(xs_by_sample[sn])
            point = {
                "color": scale.get_colour(si),
                "showlegend": False,
                "marker": {
                    "symbol": "diamond",
                    "line": {"color": "black", "width": 1},
                },
            }
            if sn in yx_by_sample:
                y = float(yx_by_sample[sn])
                point["data"] = [[x, y]]
                point["name"] = f"{sn}: actual read count vs. deduplicated read count (externally calculated)"
                series.append(point)
                y = y_tt_lbl.replace("point.y", "y").format(y=y)
                log.debug(f"Real counts for {sn}: {x}{x_suffix}, ({y}{y_suffix})")
            else:
                xs = sorted(data[sn].keys())
                ys = sorted(data[sn].values())
                if x > max(xs):
                    log.warning(
                        f"Total reads for {sn} ({x}{x_suffix}) > max preseq value ({max(xs)}{x_suffix}): "
                        "skipping this point"
                    )
                else:
                    interp_y = np.interp(x, xs, ys)
                    point["data"] = [[x, interp_y]]
                    point["name"] = sn + ": actual read count (externally calculated)"
                    series.append(point)
                    y = y_tt_lbl.replace("point.y", "y").format(y=interp_y)
                    log.debug(f"Real count for {sn}: {x}{{x_suffix}} ({y}{{y_suffix}})")
    return series
