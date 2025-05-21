from multiqc.base_module import BaseMultiqcModule
import re
import logging

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="BamUtil clipOverlap",
            anchor="bamutil_clipoverlap",
            href="http://genome.sph.umich.edu/wiki/BamUtil",
            info="bamUtil is a repository that contains several programs that perform operations on SAM/BAM files. The clipOverlap submodule clips overlapping read pairs to prevent double-counting of coverage.",
            doi="10.1186/s13059-014-0559-z",
        )
        self.bamutil_data = self.parse_logs()
        if len(self.bamutil_data) == 0:
            log.debug("Could not find any bamUtil clipOverlap reports")
            raise UserWarning

        # Write parsed data to a file
        self.write_data_file(self.bamutil_data, "multiqc_bamutil_clipoverlap")

        # Add software version placeholder
        self.add_software_version(None)

        self.bamutil_general_stats()

    def parse_logs(self):
        data = dict()
        for f in self.find_log_files("bamutil_clipoverlap", filehandles=True):
            s_name = self.clean_s_name(f["s_name"])

            # Initialize all potential metrics to None
            number_overlapping_pairs = None
            avg_ref_bases_overlapped = None
            variance_ref_bases = None
            orientation_additional_clipping = None
            forward_strand_clipped = None
            reverse_strand_clipped = None
            missing_overlapping_mates = None
            has_warnings = False
            for line in f["f"]:
                line = line.strip()

                # Check for warning messages
                if line.startswith("WARNING"):
                    has_warnings = True
                    if "did not find expected overlapping mates" in line:
                        # Extract the number of missing overlapping mates
                        match = re.search(r"did not find expected overlapping mates for (\d+) records", line)
                        if match:
                            missing_overlapping_mates = int(match.group(1))
                    continue

                # Skip problematic lines or info messages
                if line.startswith("Problems encountered") or "Command line parameter" in line or not line:
                    continue

                # Parse statistics
                if "Number of overlapping pairs:" in line:
                    number_overlapping_pairs = int(line.split(":")[1].strip())
                elif "Average # Reference Bases Overlapped:" in line:
                    avg_ref_bases_overlapped = float(line.split(":")[1].strip())
                elif "Variance of Reference Bases overlapped:" in line:
                    variance_ref_bases = float(line.split(":")[1].strip())
                elif "Number of times orientation causes additional clipping:" in line:
                    orientation_additional_clipping = int(line.split(":")[1].strip())
                elif "Number of times the forward strand was clipped:" in line:
                    forward_strand_clipped = int(line.split(":")[1].strip())
                elif "Number of times the reverse strand was clipped:" in line:
                    reverse_strand_clipped = int(line.split(":")[1].strip())

            # Only add the sample if we have at least the mandatory metrics
            if number_overlapping_pairs is not None and avg_ref_bases_overlapped is not None:
                data[s_name] = {
                    "bamutil_overlapping_pairs": number_overlapping_pairs,
                    "bamutil_avg_ref_bases_overlapped": avg_ref_bases_overlapped,
                }

                # Add optional metrics if they exist
                if variance_ref_bases is not None:
                    data[s_name]["bamutil_variance_ref_bases"] = variance_ref_bases
                if orientation_additional_clipping is not None:
                    data[s_name]["bamutil_orientation_additional_clipping"] = orientation_additional_clipping
                if forward_strand_clipped is not None:
                    data[s_name]["bamutil_forward_strand_clipped"] = forward_strand_clipped
                if reverse_strand_clipped is not None:
                    data[s_name]["bamutil_reverse_strand_clipped"] = reverse_strand_clipped
                if missing_overlapping_mates is not None:
                    data[s_name]["bamutil_missing_overlapping_mates"] = missing_overlapping_mates
                data[s_name]["bamutil_has_warnings"] = has_warnings

                # Add data source
                self.add_data_source(f, s_name)
            else:
                log.warning(f"Could not parse bamUtil clipOverlap stats for {s_name}")

        return data

    def bamutil_general_stats(self):
        """Add clipOverlap statistics to the general stats table"""
        headers = {
            "bamutil_overlapping_pairs": {
                "title": "Overlapping Pairs",
                "description": "Number of read pairs with overlapping segments (bamUtil clipOverlap)",
                "min": 0,
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "bamutil_avg_ref_bases_overlapped": {
                "title": "Avg Overlap",
                "description": "Average number of reference bases overlapped per read pair (bamUtil clipOverlap)",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Purples",
            },
            "bamutil_forward_strand_clipped": {
                "title": "Fwd Clipped",
                "description": "Number of times the forward strand was clipped (bamUtil clipOverlap)",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Greens",
                "hidden": True,
            },
            "bamutil_reverse_strand_clipped": {
                "title": "Rev Clipped",
                "description": "Number of times the reverse strand was clipped (bamUtil clipOverlap)",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Greens",
                "hidden": True,
            },
            "bamutil_orientation_additional_clipping": {
                "title": "Orientation Clip",
                "description": "Number of times orientation causes additional clipping (bamUtil clipOverlap)",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "OrRd",
                "hidden": True,
            },
            "bamutil_variance_ref_bases": {
                "title": "Overlap Variance",
                "description": "Variance of reference bases overlapped (bamUtil clipOverlap)",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "RdPu",
                "hidden": True,
            },
            "bamutil_missing_overlapping_mates": {
                "title": "Missing Mates",
                "description": "Number of reads where expected overlapping mates were not found (bamUtil clipOverlap)",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Reds",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.bamutil_data, headers)

        # If we have multiple samples, add a section with a plot
        if len(self.bamutil_data) > 1:
            self.add_clipping_barplot()

    def add_clipping_barplot(self):
        """Create a bar plot showing clipping statistics"""
        # Define categories and their display names
        categories = {
            "bamutil_forward_strand_clipped": {"name": "Forward Strand Clipped"},
            "bamutil_reverse_strand_clipped": {"name": "Reverse Strand Clipped"},
            "bamutil_orientation_additional_clipping": {"name": "Additional Clipping from Orientation"},
        }

        # Only include samples with all required data
        plot_data = {}
        for s_name, sample_data in self.bamutil_data.items():
            has_all_metrics = all(metric in sample_data for metric in categories.keys())
            if has_all_metrics:
                plot_data[s_name] = sample_data

        if len(plot_data) > 0:
            # Plot configuration
            config = {
                "id": "bamutil_clipping_plot",
                "title": "BamUtil: Clipping Distribution",
                "ylab": "Number of Reads",
                "cpswitch_counts_label": "Number of Reads Clipped",
                "hide_zero_cats": False,
            }

            # Create the section with the plot
            from multiqc.plots import bargraph

            self.add_section(
                name="Clipping Distribution",
                anchor="bamutil_clipping_dist",
                description="Distribution of forward vs. reverse strand clipping and orientation-based additional clipping",
                helptext="""
                This plot shows the distribution of where clipping was performed:
                
                * **Forward Strand Clipped**: Number of times the forward strand was clipped
                * **Reverse Strand Clipped**: Number of times the reverse strand was clipped  
                * **Additional Clipping from Orientation**: Number of times read orientation caused additional clipping
                
                Differences between forward and reverse strand clipping may indicate sequencing biases or reference issues.
                """,
                plot=bargraph.plot(plot_data, categories, config),
            )
