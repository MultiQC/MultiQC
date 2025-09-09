import logging
import pandas as pd

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, box, linegraph, scatter, table
from multiqc.plots.table_object import ColumnDict, TableConfig
from multiqc.utils import mqc_colour


log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
   
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Space Ranger",
            anchor="spaceranger",
            href=["XXX"],
            info="XXX",
            doi=["XX"],
        )

        data_by_sample = {}
        for f in self.find_log_files("spaceranger", filehandles=True):
            parsed_data = self.parse_spaceranger_metrics(f)
            if parsed_data:
                sample_name = parsed_data["Sample ID"]
                data_by_sample[sample_name] = parsed_data
                self.add_data_source(f, sample_name)
        
        data_by_sample = self.ignore_samples(data_by_sample)
        log.info(f"Found {len(data_by_sample)} Space Ranger reports")
        self.write_data_file(data_by_sample, "multiqc_spaceranger")
        self.spaceranger_general_stats_table(data_by_sample)

        headers: Dict[str, Dict[str, Any]] = {
            "Number of Reads": {
                "title": "Number of Reads",
                "description": "Total number of reads sequenced",
                "scale": "YlOrRd",
                "format": "{:,.0f}",
            },
             "Reads Mapped to Probe Set": {
                "title": "Reads Mapped to Probe Set",
                "description": "Reads Mapped to Probe Set",
                "scale": "YlOrRd",
            },
             "Number of Genes": {
                "title": "Number of Genes",
                "description": "Number of Genes",
                "scale": "YlOrRd",
            },
             "Estimated UMIs from Genomic DNA": {
                "title": "Estimated UMIs from Genomic DNA",
                "description": "Estimated UMIs from Genomic DNA",
                "scale": "YlOrRd",
                "format": "{:,.0f}",
            },
            }

        self.add_section(
                    name="Table",
                    anchor="Table",
                    description="XXX",
                    helptext="""
                    XXX
                    """,
                    plot=table.plot(
                                    data_by_sample,
                                    headers,
                                    pconfig=TableConfig(
                                        id="Space Ranger table",
                                        title="Space Ranger table: Data Quality",
                                    ))
        )
            
    def parse_spaceranger_metrics(self, f):
        """Parse Space Ranger metrics_summary.csv file"""
        lines = f["f"].read().splitlines()
        if len(lines) < 2:
            return {}

        # Get header and data row
        header = lines[0].split(",")
        data_row = lines[1].split(",")

        if len(header) != len(data_row):
            log.warning(f"Header and data row lengths don't match in {f['fn']}")
            return {}

        # Create mapping of headers to values
        metrics = dict(zip(header, data_row))

        # Convert numeric fields
        numeric_fields = [
            "Number of Reads",
            "Valid Barcodes",
            "Valid UMI Sequences",	
            "Sequencing Saturation", 
            "Q30 Bases in Barcode",
            "Q30 Bases in Probe Read",
            "Q30 Bases in UMI",
            "Reads Mapped to Probe Set",
            "Reads Mapped Confidently to Probe Set",
            "Fraction Reads in Squares Under Tissue",
            "Genes Detected",
            "Reads Mapped Confidently to the Filtered Probe Set",
            "Number of Genes", 
            "Reads Half-Mapped to Probe Set",
            "Reads Split-Mapped to Probe Set",
            "Estimated UMIs from Genomic DNA",
            "Estimated UMIs from Genomic DNA per Unspliced Probe",
            "Number of Squares Under Tissue 2 µm",
            "Mean Reads Under Tissue per Square 2 µm",	
            "Fraction of Squares Under Tissue 2 µm",
            "Mean Genes Under Tissue per Square 2 µm",	
            "Mean UMIs Under Tissue per Square 2 µm",
            "Total Genes Detected Under Tissue 2 µm",
            "Number of Bins Under Tissue 8 µm",
            "Mean Reads Under Tissue per Bin 8 µm",	
            "Fraction of Bins Under Tissue 8 µm",
            "Mean Genes Under Tissue per Bin 8 µm",	
            "Mean UMIs Under Tissue per Bin 8 µm",
            "Total Genes Detected Under Tissue 8 µm",
            "UMIs per sq mm of Tissue",
            "Reads per sq mm of Tissue",
            "Number of Bins Under Tissue 16 µm",
            "Mean Reads Under Tissue per Bin 16 µm",
            "Fraction of Bins Under Tissue 16 µm",
            "Mean Genes Under Tissue per Bin 16 µm",
            "Mean UMIs Under Tissue per Bin 16 µm",
            "Total Genes Detected Under Tissue 16 µm",
            "Number of Cells",
            "Reads in Cells",
            "UMIs in Cells",
            "Mean Reads per Cell",
            "Median Genes per Cell",
            "Median UMIs per Cell",
            "Median Cell Area (μm²)"
            "Median Nucleus Area (μm²)",
            "Maximum Nucleus Diameter (pixels)"
        ]

        parsed_metrics = {}
        for field in numeric_fields:
            if field in metrics and metrics[field] != "":
                try:
                    parsed_metrics[field] = float(metrics[field])
                except ValueError:
                    log.warning(f"Could not convert {field}='{metrics[field]}' to float")

        # Keep string fields
        string_fields = [
            "Sample ID"
        ]
        for field in string_fields:
            if field in metrics:
                parsed_metrics[field] = metrics[field]

        return parsed_metrics

    def spaceranger_general_stats_table(self, data_by_sample):
        """Add key Xenium metrics to the general statistics table"""
        headers: Dict[str, Dict[str, Any]] = {
            "Number of Reads": {
                "title": "Number of Reads",
                "description": "Total number of reads sequenced",
                "scale": "YlOrRd",
                "format": "{:,.0f}",
            },
             "Reads Mapped to Probe Set": {
                "title": "Reads Mapped to Probe Set",
                "description": "Reads Mapped to Probe Set",
                "scale": "YlOrRd",
            },
             "Number of Genes": {
                "title": "Number of Genes",
                "description": "Number of Genes",
                "scale": "YlOrRd",
            },
             "Estimated UMIs from Genomic DNA": {
                "title": "Estimated UMIs from Genomic DNA",
                "description": "Estimated UMIs from Genomic DNA",
                "scale": "YlOrRd",
                "format": "{:,.0f}",
            },
            }
        self.general_stats_addcols(data_by_sample, headers)
"""            ,
            "fraction_transcripts_assigned": {
                "title": "Transcripts Assigned",
                "description": "Fraction of transcripts assigned to cells",
                "suffix": "%",
                "scale": "RdYlGn",
                "modify": lambda x: x * 100.0,
                "max": 100.0,
            },
            "median_genes_per_cell": {
                "title": "Genes/Cell",
                "description": "Median number of genes per cell",
                "scale": "Purples",
                "format": "{:,.0f}",
            },
            "fraction_transcripts_decoded_q20": {
                "title": "Q20+ Transcripts",
                "description": "Fraction of transcripts decoded with Q20+",
                "suffix": "%",
                "scale": "Greens",
                "modify": lambda x: x * 100.0,
                "max": 100.0,
            },
            "cell_area_median": {
                "title": "Median Cell",
                "description": "Median cell area",
                "suffix": " μm²",
                "scale": "Blues",
                "format": "{:,.1f}",
                "shared_key": "xenium_cell_area",
            },
            "nucleus_area_median": {
                "title": "Median Nucleus",
                "description": "Median nucleus area",
                "suffix": " μm²",
                "scale": "Oranges",
                "format": "{:,.1f}",
                "shared_key": "xenium_cell_area",
            },
            "nucleus_to_cell_area_ratio_median": {
                "title": "Nucleus/Cell",
                "description": "Median nucleus to cell area ratio",
                "scale": "Greens",
                "format": "{:.3f}",
                "max": 1.0,
            }, """
