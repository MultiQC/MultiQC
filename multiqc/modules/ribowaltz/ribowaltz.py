"""MultiQC module for riboWaltz ribosome profiling QC"""

import logging
from typing import Dict, List, Optional

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph
from multiqc.plots.bargraph import BarPlotConfig
from multiqc.plots.linegraph import LinePlotConfig

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    riboWaltz computes P-site offsets and performs quality control for ribosome profiling data.

    The module parses QC output files from riboWaltz and generates visualizations for:

    - **P-site region distribution**: Shows where P-sites map across transcript regions
      (5' UTR, CDS, 3' UTR). Good Ribo-seq data shows >70% CDS enrichment.
    - **Reading frame distribution**: Shows P-site distribution across reading frames
      for each transcript region. Frame 0 should be >50% in CDS but not in UTRs.
    - **Metaprofiles**: Shows P-site frequency around start and stop codons.
      Good data shows trinucleotide periodicity with Frame 0 peaks.

    Supported input files:

    - `*ribowaltz*psite_region.tsv` - P-site region distribution
    - `*ribowaltz*frames.tsv` - Reading frame distribution
    - `*ribowaltz*metaprofile_psite.tsv` - Metaprofile around start/stop codons

    Files must contain "ribowaltz" in the filename since headers are generic.
    Both tab-delimited and comma-delimited files are supported.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="riboWaltz",
            anchor="ribowaltz",
            href="https://github.com/LabTranslationalArchitectomics/riboWaltz",
            info="Computes P-site offsets and performs quality control for ribosome profiling data.",
            doi="10.1371/journal.pcbi.1006169",
        )

        # Data storage
        self.psite_region_data: Dict[str, Dict] = {}
        self.frames_data: Dict[str, Dict[str, Dict[str, float]]] = {}
        self.metaprofile_start_data: Dict[str, Dict[int, float]] = {}
        self.metaprofile_stop_data: Dict[str, Dict[int, float]] = {}
        self.rna_reference_data: Optional[Dict[str, float]] = None

        # Parse all file types
        self.parse_psite_region_files()
        self.parse_frames_files()
        self.parse_metaprofile_files()

        # Filter ignored samples
        self.psite_region_data = self.ignore_samples(self.psite_region_data)
        self.frames_data = self.ignore_samples(self.frames_data)
        self.metaprofile_start_data = self.ignore_samples(self.metaprofile_start_data)
        self.metaprofile_stop_data = self.ignore_samples(self.metaprofile_stop_data)

        # Check if we found anything
        n_samples = max(
            len(self.psite_region_data),
            len(self.frames_data),
            len(self.metaprofile_start_data),
            len(self.metaprofile_stop_data),
        )
        if n_samples == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {n_samples} samples")

        # Add software version (riboWaltz doesn't include version in output files)
        self.add_software_version(None)

        # Add general stats table
        self.ribowaltz_general_stats_table()

        # Add sections (in order)
        if self.psite_region_data:
            self.psite_region_bargraph()
        if self.frames_data:
            self.frames_bargraph()
        if self.metaprofile_start_data:
            self.metaprofile_start_linegraph()
        if self.metaprofile_stop_data:
            self.metaprofile_stop_linegraph()

        # Write data files at the end
        if self.psite_region_data:
            self.write_data_file(self.psite_region_data, "multiqc_ribowaltz_psite_region")
        if self.frames_data:
            self.write_data_file(self.frames_data, "multiqc_ribowaltz_frames")
        if self.metaprofile_start_data:
            self.write_data_file(self.metaprofile_start_data, "multiqc_ribowaltz_metaprofile_start")
        if self.metaprofile_stop_data:
            self.write_data_file(self.metaprofile_stop_data, "multiqc_ribowaltz_metaprofile_stop")

    def detect_delimiter(self, line: str) -> Optional[str]:
        """Auto-detect delimiter from header line"""
        if "\t" in line:
            return "\t"
        elif "," in line:
            return ","
        return None

    def parse_psite_region_files(self):
        """Parse psite_region.tsv files"""
        for f in self.find_log_files("ribowaltz/psite_region", filehandles=True):
            self.parse_psite_region(f)

    def parse_psite_region(self, f):
        """Parse a single psite_region.tsv file"""
        lines = f["f"].readlines()
        if len(lines) < 2:
            return

        # Detect delimiter from header
        header_line = lines[0].strip()
        delim = self.detect_delimiter(header_line)
        if delim is None:
            log.warning(f"Could not detect delimiter in {f['fn']}")
            return

        headers = header_line.split(delim)

        # Verify expected headers
        required_headers = ["sample", "region", "scaled_count"]
        if not all(h in headers for h in required_headers):
            log.debug(f"Unexpected headers in {f['fn']}: {headers}")
            return

        # Parse data lines
        sample_data: Dict[str, Dict[str, float]] = {}
        rna_ref_data: Dict[str, float] = {}

        for line in lines[1:]:
            parts = line.strip().split(delim)
            if len(parts) != len(headers):
                continue

            row = dict(zip(headers, parts))
            sample = row["sample"]
            region = row["region"]
            try:
                scaled_count = float(row["scaled_count"])
            except ValueError:
                continue

            # Separate sample data from RNA reference
            if sample == "RNAs":
                rna_ref_data[region] = scaled_count
            else:
                if sample not in sample_data:
                    sample_data[sample] = {}
                sample_data[sample][region] = scaled_count

        # Store RNA reference data (only once)
        if rna_ref_data and self.rna_reference_data is None:
            self.rna_reference_data = rna_ref_data

        # Store parsed data (handle multiple samples per file)
        for sample, data in sample_data.items():
            clean_name = self.clean_s_name(sample, f)
            if clean_name in self.psite_region_data:
                log.debug(f"Duplicate sample name found! Overwriting: {clean_name}")

            self.add_data_source(f, clean_name)
            self.psite_region_data[clean_name] = data

    def parse_frames_files(self):
        """Parse frames.tsv files"""
        for f in self.find_log_files("ribowaltz/frames", filehandles=True):
            self.parse_frames(f)

    def parse_frames(self, f):
        """Parse a single frames.tsv file"""
        lines = f["f"].readlines()
        if len(lines) < 2:
            return

        # Detect delimiter from header
        header_line = lines[0].strip()
        delim = self.detect_delimiter(header_line)
        if delim is None:
            log.warning(f"Could not detect delimiter in {f['fn']}")
            return

        headers = header_line.split(delim)

        # Verify expected headers
        required_headers = ["sample", "region", "frame", "scaled_count"]
        if not all(h in headers for h in required_headers):
            log.debug(f"Unexpected headers in {f['fn']}: {headers}")
            return

        # Parse data lines: sample -> region -> frame -> scaled_count
        sample_data: Dict[str, Dict[str, Dict[str, float]]] = {}

        for line in lines[1:]:
            parts = line.strip().split(delim)
            if len(parts) != len(headers):
                continue

            row = dict(zip(headers, parts))
            sample = row["sample"]
            region = row["region"]
            try:
                frame = int(row["frame"])
                scaled_count = float(row["scaled_count"])
            except ValueError:
                continue

            frame_key = f"Frame {frame}"

            if sample not in sample_data:
                sample_data[sample] = {}
            if region not in sample_data[sample]:
                sample_data[sample][region] = {}
            sample_data[sample][region][frame_key] = scaled_count

        # Store parsed data
        for sample, data in sample_data.items():
            clean_name = self.clean_s_name(sample, f)
            if clean_name in self.frames_data:
                log.debug(f"Duplicate sample name found! Overwriting: {clean_name}")

            self.add_data_source(f, clean_name)
            self.frames_data[clean_name] = data

    def parse_metaprofile_files(self):
        """Parse metaprofile_psite.tsv files"""
        for f in self.find_log_files("ribowaltz/metaprofile", filehandles=True):
            self.parse_metaprofile(f)

    def parse_metaprofile(self, f):
        """Parse a single metaprofile_psite.tsv file"""
        lines = f["f"].readlines()
        if len(lines) < 2:
            return

        header_line = lines[0].strip()
        delim = self.detect_delimiter(header_line)
        if delim is None:
            log.warning(f"Could not detect delimiter in {f['fn']}")
            return

        headers = header_line.split(delim)

        # Verify expected headers
        required_headers = ["sample", "region", "x", "y"]
        if not all(h in headers for h in required_headers):
            log.debug(f"Unexpected headers in {f['fn']}: {headers}")
            return

        # Separate start and stop codon data
        start_data: Dict[str, Dict[int, float]] = {}
        stop_data: Dict[str, Dict[int, float]] = {}

        for line in lines[1:]:
            parts = line.strip().split(delim)
            if len(parts) != len(headers):
                continue

            row = dict(zip(headers, parts))
            sample = row["sample"]
            region = row["region"]
            try:
                x = int(float(row["x"]))  # Convert to int for linegraph x-axis
                y = float(row["y"])
            except ValueError:
                continue

            if "start" in region.lower():
                if sample not in start_data:
                    start_data[sample] = {}
                start_data[sample][x] = y
            elif "stop" in region.lower():
                if sample not in stop_data:
                    stop_data[sample] = {}
                stop_data[sample][x] = y

        # Store parsed data
        for sample, data in start_data.items():
            clean_name = self.clean_s_name(sample, f)
            self.add_data_source(f, clean_name)
            if clean_name in self.metaprofile_start_data:
                log.debug(f"Duplicate sample name found! Overwriting: {clean_name}")
            self.metaprofile_start_data[clean_name] = data

        for sample, data in stop_data.items():
            clean_name = self.clean_s_name(sample, f)
            if clean_name in self.metaprofile_stop_data:
                log.debug(f"Duplicate sample name found! Overwriting: {clean_name}")
            self.metaprofile_stop_data[clean_name] = data

    def ribowaltz_general_stats_table(self):
        """Add key metrics to the general stats table"""
        headers = {
            "cds_pct": {
                "title": "CDS %",
                "description": "Percentage of P-sites mapping to CDS",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.1f}",
            },
            "frame0_cds_pct": {
                "title": "Frame 0 %",
                "description": "Percentage of CDS P-sites in Frame 0",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.1f}",
            },
        }

        data: dict = {}

        # Add CDS percentage from psite_region_data
        for s_name, regions in self.psite_region_data.items():
            if s_name not in data:
                data[s_name] = {}
            if "CDS" in regions:
                data[s_name]["cds_pct"] = regions["CDS"]

        # Add Frame 0 percentage from frames_data
        for s_name, region_frames in self.frames_data.items():
            if s_name not in data:
                data[s_name] = {}
            cds_frames = region_frames.get("CDS", {})
            if "Frame 0" in cds_frames:
                data[s_name]["frame0_cds_pct"] = cds_frames["Frame 0"]

        if data:
            self.general_stats_addcols(data, headers)

    def psite_region_bargraph(self):
        """Create stacked bar graph for P-site region distribution"""

        # Prepare data for plotting
        plot_data: Dict[str, Dict[str, float]] = {}

        for s_name, regions in self.psite_region_data.items():
            plot_data[s_name] = regions

        # Add RNA reference as special sample if available
        if self.rna_reference_data:
            plot_data["RNA-seq reference"] = self.rna_reference_data

        # Define categories (order matters for stacked bar)
        cats = ["5' UTR", "CDS", "3' UTR"]

        pconfig = BarPlotConfig(
            id="ribowaltz_psite_regions",
            title="riboWaltz: P-site Region Distribution",
            ylab="% of P-sites",
            cpswitch=False,
            ymax=100,
            ymin=0,
        )

        self.add_section(
            name="P-site Region Distribution",
            anchor="ribowaltz_psite_regions",
            description=(
                "Distribution of P-sites across transcript regions. "
                "Good Ribo-seq data shows strong CDS enrichment (>70%). "
                "The RNA-seq reference shows expected distribution from uniform transcript coverage."
            ),
            helptext="""
This plot shows where ribosome-protected fragments (P-sites) map across different
transcript regions:

- **5' UTR**: Upstream of the coding sequence
- **CDS**: Coding sequence where translation occurs
- **3' UTR**: Downstream of the coding sequence

**Quality metrics:**

- Good Ribo-seq: >70% CDS enrichment
- RNA-seq reference: More evenly distributed (~46% CDS)

Strong CDS enrichment indicates that the ribosome profiling captured actively
translating ribosomes rather than random RNA fragments.
            """,
            plot=bargraph.plot(plot_data, cats, pconfig),
        )

    def frames_bargraph(self):
        """Create grouped stacked bar graph for reading frame distribution"""

        # Prepare data with special sample keys for grouping
        plot_data: Dict[str, Dict[str, float]] = {}
        samples = sorted(self.frames_data.keys())
        regions = ["5' UTR", "CDS", "3' UTR"]

        for s_name in samples:
            for region in regions:
                if region in self.frames_data[s_name]:
                    key = f"{s_name}_{region}"
                    plot_data[key] = self.frames_data[s_name][region]

        # Define frame categories
        cats = ["Frame 0", "Frame 1", "Frame 2"]

        # Build sample groups for visual grouping by region
        # Format: {"Group Name": [[sample_key, display_name], ...], ...}
        sample_groups: Dict[str, List[List[str]]] = {}
        for region in regions:
            group_samples: List[List[str]] = []
            for s_name in samples:
                key = f"{s_name}_{region}"
                if key in plot_data:
                    group_samples.append([key, s_name])
            if group_samples:
                sample_groups[region] = group_samples

        pconfig = BarPlotConfig(
            id="ribowaltz_frames",
            title="riboWaltz: Reading Frame Distribution",
            ylab="% of P-sites",
            cpswitch=False,
            ymax=100,
            ymin=0,
            sample_groups=sample_groups,
            use_legend=True,
        )

        self.add_section(
            name="Reading Frame Distribution",
            anchor="ribowaltz_frames",
            description=(
                "Distribution of P-sites across reading frames for each transcript region. "
                "Good Ribo-seq data shows Frame 0 enrichment (>50%) in the CDS but not in UTRs."
            ),
            helptext="""
This plot shows the distribution of P-sites across the three possible reading frames (0, 1, 2)
for each transcript region.

**Quality metrics:**

- **CDS Frame 0**: Should be >50% (ideally >60%)
- **UTR Frame 0**: Should NOT be enriched (should be ~33%)

Strong Frame 0 enrichment in the CDS indicates proper ribosome positioning and
successful P-site identification. UTRs should not show frame preference because
they are not translated.

Samples are grouped by region (5' UTR, CDS, 3' UTR) for easy comparison.
            """,
            plot=bargraph.plot(plot_data, cats, pconfig),
        )

    def metaprofile_start_linegraph(self):
        """Create line graph for metaprofile around start codon"""

        pconfig = LinePlotConfig(
            id="ribowaltz_metaprofile_start",
            title="riboWaltz: Metaprofile (Start Codon)",
            xlab="Distance from start codon (nt)",
            ylab="P-site frequency",
            x_decimals=False,
            ymin=0,
        )

        self.add_section(
            name="Metaprofile (Start Codon)",
            anchor="ribowaltz_metaprofile_start",
            description=(
                "P-site frequency around the start codon. "
                "Good Ribo-seq data shows trinucleotide periodicity with peaks at frame 0 positions."
            ),
            helptext="""
This plot shows the distribution of P-sites relative to the start codon position.

**Quality metrics:**

- Should show clear **3-nucleotide periodicity**
- Peaks should align with **Frame 0** positions (0, 3, 6, 9, etc.)
- Strong peak at or near the start codon

Clear periodicity indicates:

1. Successful P-site offset calculation
2. Proper ribosome footprint protection
3. Active translation initiation
            """,
            plot=linegraph.plot(self.metaprofile_start_data, pconfig),
        )

    def metaprofile_stop_linegraph(self):
        """Create line graph for metaprofile around stop codon"""

        pconfig = LinePlotConfig(
            id="ribowaltz_metaprofile_stop",
            title="riboWaltz: Metaprofile (Stop Codon)",
            xlab="Distance from stop codon (nt)",
            ylab="P-site frequency",
            x_decimals=False,
            ymin=0,
        )

        self.add_section(
            name="Metaprofile (Stop Codon)",
            anchor="ribowaltz_metaprofile_stop",
            description=(
                "P-site frequency around the stop codon. "
                "Good Ribo-seq data shows trinucleotide periodicity and proper termination patterns."
            ),
            helptext="""
This plot shows the distribution of P-sites relative to the stop codon position.

**Quality metrics:**

- Should show **3-nucleotide periodicity** upstream of stop
- Periodicity should decrease downstream (after termination)
- May show accumulation at the stop codon

This pattern confirms proper translation termination and validates the
overall quality of the ribosome profiling experiment.
            """,
            plot=linegraph.plot(self.metaprofile_stop_data, pconfig),
        )
