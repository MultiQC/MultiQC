import logging
import re
import traceback

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

log = logging.getLogger(__name__)

VERSION_REGEX = r"\[FLASH\] FLASH v([\d\.]+) complete!"


class MultiqcModule(BaseMultiqcModule):
    """
    To create a log file suitable for the module, you can use `tee`. From the FLASh help:

    ```bash
    flash reads_1.fq reads_2.fq 2>&1 | tee logfilename.log
    ```

    The sample name is set by the first input filename listed in the log. However, this can be changed to using the first output filename (i.e. if you used FLASh's `--output-prefix=PREFIX` option) by using the following config:

    ```yaml
    flash:
      use_output_name: true
    ```

    The module can also parse the `.hist` numeric histograms output by FLASh.

    Note that the histogram's file format and extension are too generic by themselves which could result in the accidental parsing a file output by another tool. To get around this, the MultiQC module only parses files with the filename pattern `*flash*.hist`.

    To customise this (for example, enabling for any file ending in `*.hist`), use the following config change:

    ```yaml
    sp:
      flash/hist:
        fn: "*.hist"
    ```
    """

    """FLASh MultiQC module

    Options:
      use_output_name: true - use first output filename as sample name
        default uses first input filename in log
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="FLASh",
            anchor="flash",
            href="https://ccb.jhu.edu/software/FLASH/",
            info="Merges paired-end reads from next-generation sequencing experiments.",
            doi="10.1093/bioinformatics/btr507",
        )

        flash_results = self.parse_flash()

        hist_results = self.hist_results()

        if not flash_results and not hist_results:
            raise ModuleNoSamplesFound

    def parse_flash(self):
        # Find all log files with flash msgs
        flash_data = {}
        for logfile in self.find_log_files("flash/log"):
            flash_data.update(self.parse_flash_log(logfile))
            self.add_data_source(logfile)

        # Ignore sample names
        flash_data = self.ignore_samples(flash_data)
        if not flash_data:
            return 0

        log.info(f"Found {len(flash_data)} log reports")

        self.general_stats_table(flash_data)
        self.add_section(
            name="Read combination statistics", anchor="flash-bargraph", plot=self.summary_plot(flash_data)
        )

        self.write_data_file(flash_data, "multiqc_flash_combo_stats")
        return len(flash_data)

    @staticmethod
    def split_log(logf):
        """split concat log into individual samples"""
        flashpatt = re.compile(
            r"\[FLASH\] Fast Length Adjustment of SHort reads\n(.+?)\[FLASH\] FLASH", flags=re.DOTALL
        )
        return flashpatt.findall(logf)

    @staticmethod
    def get_field(field, slog, fl=False):
        """parse sample log for field
        set fl=True to return a float
        otherwise, returns int
        """
        field += r"\:\s+([\d\.]+)"
        match = re.search(field, slog)
        if match:
            if fl:
                return float(match.group(1))
            return int(match.group(1))
        return 0

    def clean_pe_name(self, nlog, logf):
        """additional name cleaning for paired end data"""
        use_output_name = getattr(config, "flash", {}).get("use_output_name", False)
        if use_output_name:
            name = re.search(r"Output files\:\n\[FLASH\]\s+(.+?)\n", nlog)
        else:
            name = re.search(r"Input files\:\n\[FLASH\]\s+(.+?)\n", nlog)
        if not name:
            return None
        name = name.group(1)
        name = self.clean_s_name(name, logf)
        return name

    def parse_flash_log(self, logf):
        """parse flash logs"""
        data = {}
        samplelogs = self.split_log(logf["f"])
        version_match = re.search(VERSION_REGEX, logf["f"])

        for slog in samplelogs:
            try:
                sample = dict()
                ## Sample name ##
                s_name = self.clean_pe_name(slog, logf)
                if s_name is None:
                    continue
                sample["s_name"] = s_name

                if version_match:
                    self.add_software_version(version_match.group(1), s_name)

                ## Log attributes ##
                sample["totalpairs"] = self.get_field("Total pairs", slog)
                sample["discardpairs"] = self.get_field("Discarded pairs", slog)
                sample["percdiscard"] = self.get_field("Percent Discarded", slog, fl=True)
                sample["combopairs"] = self.get_field("Combined pairs", slog)
                sample["inniepairs"] = self.get_field("Innie pairs", slog)
                sample["outiepairs"] = self.get_field("Outie pairs", slog)
                sample["uncombopairs"] = self.get_field("Uncombined pairs", slog)
                sample["perccombo"] = self.get_field("Percent combined", slog, fl=True)

                data[s_name] = sample
            except Exception as err:
                log.warning(f"Error parsing record in {logf['fn']}. {err}")
                log.debug(traceback.format_exc())
                continue
        return data

    def general_stats_table(self, data):
        """Add percent combined to general stats table"""
        headers = {
            "combopairs": {
                "title": "Combined pairs",
                "description": "Num read pairs combined",
                "shared_key": "read_count",
                "hidden": True,
                "scale": False,
            },
            "perccombo": {
                "title": "% Combined",
                "description": "% read pairs combined",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "PiYG",
            },
        }
        self.general_stats_addcols(data, headers)

    @staticmethod
    def summary_plot(data):
        """Barplot of combined pairs"""
        cats = {
            "combopairs": {"name": "Combined pairs", "color": "#191970"},
            "inniepairs": {"name": "Combined innie pairs", "color": "#191970"},
            "outiepairs": {"name": "Combined outie pairs", "color": "#00A08A"},
            "uncombopairs": {"name": "Uncombined pairs", "color": "#cd1076"},
            "discardpairs": {"name": "Discarded pairs", "color": "#ffd700"},
        }
        splotconfig = {
            "id": "flash_combo_stats_plot",
            "title": "FLASh: Read combination statistics",
            "ylab": "Number of read pairs",
        }
        # Only plot the combopairs category if we don't have inniepairs and outiepairs
        for s_name, d in data.items():
            if (d["inniepairs"] > 0 or d["outiepairs"] > 0) and d["combopairs"] > 0:
                del data[s_name]["combopairs"]

        return bargraph.plot(data, cats, splotconfig)

    @staticmethod
    def parse_hist_files(histf):
        """parse histogram files"""
        nameddata = dict()
        data = dict()
        try:
            for line in histf["f"].splitlines():
                s = line.split()
                if s:
                    if len(s) != 2:
                        raise RuntimeError(
                            "invalid format: " + str(len(s)) + " column(s) found in row. must be exactly 2."
                        )
                    data[int(s[0])] = int(s[1])
        except Exception as err:
            log.warning("Error parsing %s. %s", histf["fn"], err)
            log.debug(traceback.format_exc())
        else:
            if data:
                nameddata[histf["s_name"]] = data
            else:
                log.debug("%s is empty.", histf["fn"])
        finally:
            return nameddata

    @staticmethod
    def get_colors(n):
        """get colors for freqpoly graph"""
        cb_palette = [
            "#E69F00",
            "#56B4E9",
            "#009E73",
            "#F0E442",
            "#0072B2",
            "#D55E00",
            "#CC79A7",
            "#001F3F",
            "#0074D9",
            "#7FDBFF",
            "#39CCCC",
            "#3D9970",
            "#2ECC40",
            "#01FF70",
            "#FFDC00",
            "#FF851B",
            "#FF4136",
            "#F012BE",
            "#B10DC9",
            "#85144B",
            "#AAAAAA",
            "#000000",
        ]

        whole = int(n / 22)
        extra = n % 22
        cols = cb_palette * whole
        if extra >= 0:
            cols.extend(cb_palette[0:extra])
        return cols

    @staticmethod
    def freqpoly_plot(data):
        """make freqpoly plot of merged read lengths"""
        rel_data = {}
        for key, val in data.items():
            tot = sum(val.values(), 0)
            rel_data[key] = {k: v / tot for k, v in val.items()}
        fplotconfig = {
            "id": "flash_freqpoly_plot",
            "title": "FLASh: Frequency of merged read lengths",
            "xlab": "Merged Read Length",
            "ylab": "Frequency",
            "data_labels": [
                {"name": "Absolute", "ylab": "Frequency", "xlab": "Merged Read Length"},
                {"name": "Relative", "ylab": "Relative Frequency", "xlab": "Merged Read Length"},
            ],
            "colors": dict(zip(data.keys(), MultiqcModule.get_colors(len(data)))),
        }
        return linegraph.plot([data, rel_data], fplotconfig)

    def hist_results(self):
        """process flash numeric histograms"""
        hist_data = {}
        for f in self.find_log_files("flash/hist"):
            hist_data.update(self.parse_hist_files(f))

        # ignore sample names
        hist_data = self.ignore_samples(hist_data)
        if not hist_data:
            return 0

        log.info("Found %d histogram reports", len(hist_data))

        self.add_section(
            name="Frequency polygons of merged read lengths",
            anchor="flash-histogram",
            description="This plot is made from the numerical histograms output by FLASh.",
            plot=self.freqpoly_plot(hist_data),
        )
        return len(hist_data)
