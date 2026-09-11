"""MultiQC module to parse output from OrthoFinder"""

import logging
import os

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table
from multiqc.plots.table_object import ColumnDict

log = logging.getLogger(__name__)

# OrthoFinder always writes both summary tables into a directory of this name.
STATS_DIRNAME = "Comparative_Genomics_Statistics"

# Rows of Statistics_PerSpecies.tsv, up to the first blank line. Everything below
# that blank line is the "genes per-species in orthogroup" histogram, which
# describes orthogroups rather than species and so has no per-sample meaning.
PER_SPECIES_KEYS = {
    "Number of genes": "genes",
    "Number of genes in orthogroups": "genes_in_orthogroups",
    "Number of unassigned genes": "unassigned_genes",
    "Percentage of genes in orthogroups": "percent_genes_in_orthogroups",
    "Percentage of unassigned genes": "percent_unassigned_genes",
    "Number of orthogroups containing species": "orthogroups_containing_species",
    "Percentage of orthogroups containing species": "percent_orthogroups_containing_species",
    "Number of species-specific orthogroups": "species_specific_orthogroups",
    "Number of genes in species-specific orthogroups": "genes_in_species_specific_orthogroups",
    "Percentage of genes in species-specific orthogroups": "percent_genes_in_species_specific_orthogroups",
}

# Run-level rows of Statistics_Overall.tsv. That file also lists sibling output
# filenames (for example "Orthogroups file"), which are skipped by accepting only
# the keys below.
OVERALL_KEYS = {
    "Number of species": "species",
    "Number of genes": "genes",
    "Number of genes in orthogroups": "genes_in_orthogroups",
    "Number of unassigned genes": "unassigned_genes",
    "Percentage of genes in orthogroups": "percent_genes_in_orthogroups",
    "Percentage of unassigned genes": "percent_unassigned_genes",
    "Number of orthogroups": "orthogroups",
    "Number of species-specific orthogroups": "species_specific_orthogroups",
    "Number of genes in species-specific orthogroups": "genes_in_species_specific_orthogroups",
    "Percentage of genes in species-specific orthogroups": "percent_genes_in_species_specific_orthogroups",
    "Mean orthogroup size": "mean_orthogroup_size",
    "Median orthogroup size": "median_orthogroup_size",
    "G50 (assigned genes)": "g50_assigned",
    "G50 (all genes)": "g50_all",
    "O50 (assigned genes)": "o50_assigned",
    "O50 (all genes)": "o50_all",
    "Number of orthogroups with all species present": "orthogroups_all_species",
    "Number of single-copy orthogroups": "single_copy_orthogroups",
}


def parse_number(value: str) -> float:
    """Return an int where the text is integral, else a float.

    Raises ValueError on anything else, so a format change surfaces instead of
    silently producing a zero.
    """
    value = value.strip()
    try:
        return int(value)
    except ValueError:
        return float(value)


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses the two summary tables that OrthoFinder writes into
    `Comparative_Genomics_Statistics/`:

    - `Statistics_PerSpecies.tsv` becomes one MultiQC sample per species, showing
      how completely each species' genes were assigned to orthogroups.
    - `Statistics_Overall.tsv` describes the run as a whole and is shown in its
      own table, one row per run.

    Both files end with orthogroup-size distribution tables. Those are not parsed,
    because they describe orthogroups rather than samples.

    The file format is the same in OrthoFinder 2 and 3, but OrthoFinder 3 only writes
    these files when the run continues past orthogroup inference. A run stopped with
    `-og` (`--only-groups`) produces no `Comparative_Genomics_Statistics/` directory
    and so nothing for this module to find.

    Because every run writes these files under the same fixed names, and into a
    directory whose name is also fixed, runs are named after the run directory above
    `Comparative_Genomics_Statistics`.
    """

    def __init__(self):
        super().__init__(
            name="OrthoFinder",
            anchor="orthofinder",
            href="https://github.com/OrthoFinder/OrthoFinder",
            info="Infers orthogroups, orthologues and rooted gene trees for comparative genomics.",
            extra="""
            OrthoFinder identifies the orthogroups containing all genes descended
            from a single gene in the last common ancestor of the species
            considered. The proportion of genes assigned to an orthogroup is the
            usual first check that a comparative genomics run behaved sensibly.
            """,
            doi="10.1186/s13059-019-1832-y",
        )

        self.per_species_data: dict[str, dict] = {}
        for f in self.find_log_files("orthofinder/per_species"):
            for s_name, data in self._parse_per_species(f).items():
                s_name = self.clean_s_name(s_name, f)
                if s_name in self.per_species_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name=s_name)
                self.per_species_data[s_name] = data
                self.add_software_version(None, s_name)

        self.overall_data: dict[str, dict] = {}
        for f in self.find_log_files("orthofinder/overall"):
            data = self._parse_overall(f)
            if not data:
                continue
            s_name = self._run_name(f)
            if s_name in self.overall_data:
                log.debug(f"Duplicate run name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name=s_name)
            self.overall_data[s_name] = data
            self.add_software_version(None, s_name)

        self.per_species_data = self.ignore_samples(self.per_species_data)
        self.overall_data = self.ignore_samples(self.overall_data)

        if len(self.per_species_data) == 0 and len(self.overall_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.per_species_data)} species and {len(self.overall_data)} runs")

        if self.per_species_data:
            self._add_per_species_section()
            self._add_general_stats()

        if self.overall_data:
            self._add_overall_section()

        # Must come after every section has been added.
        if self.per_species_data:
            self.write_data_file(self.per_species_data, "multiqc_orthofinder_per_species")
        if self.overall_data:
            self.write_data_file(self.overall_data, "multiqc_orthofinder_overall")

    def _run_name(self, f) -> str:
        """Name an OrthoFinder run after the directory a reader would recognise it by.

        Both summary files live in a directory OrthoFinder always calls
        `Comparative_Genomics_Statistics`, so that name distinguishes nothing. The
        run directory above it carries the name the user chose, or the
        `Results_<date>` that OrthoFinder generated.
        """
        root = f["root"]
        if not root:
            return f["s_name"]
        # `root` is relative when the analysis directory was given as a relative
        # path (`multiqc .`), and then the stats directory has no parent to fall
        # back on. Resolving it first makes both spellings name the same run.
        root = os.path.abspath(root)
        if os.path.basename(root) == STATS_DIRNAME:
            root = os.path.dirname(root)
        return self.clean_s_name(os.path.basename(root), f)

    def _parse_per_species(self, f) -> dict[str, dict[str, float]]:
        """Parse Statistics_PerSpecies.tsv, where species are columns and metrics rows."""
        lines = f["f"].splitlines()
        if not lines:
            return {}

        # Header is an empty leading field followed by the species names.
        species = [s.strip() for s in lines[0].split("\t")[1:]]
        if not any(species):
            log.debug(f"No species names in header of {f['fn']}, skipping")
            return {}

        parsed: dict[str, dict[str, float]] = {s: {} for s in species if s}
        for line in lines[1:]:
            # The per-species block ends at the first blank line. A distribution
            # table follows, which must not be read as further metrics.
            if not line.strip():
                break
            fields = line.split("\t")
            key = PER_SPECIES_KEYS.get(fields[0].strip())
            if key is None:
                continue
            for idx, s_name in enumerate(species, start=1):
                if s_name:
                    parsed[s_name][key] = parse_number(fields[idx])

        return {s_name: data for s_name, data in parsed.items() if data}

    def _parse_overall(self, f) -> dict[str, float]:
        """Parse Statistics_Overall.tsv, a flat key and value list, into one run row."""
        data: dict[str, float] = {}
        for line in f["f"].splitlines():
            # A distribution table follows the first blank line. Stopping here keeps
            # its repeated column headers out of the run statistics.
            if not line.strip():
                break
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            key = OVERALL_KEYS.get(fields[0].strip())
            if key is not None:
                data[key] = parse_number(fields[1])
        return data

    def _add_general_stats(self) -> None:
        headers: dict[str, ColumnDict] = {
            "percent_genes_in_orthogroups": {
                "title": "Genes in Orthogroups",
                "description": "Percentage of this species' genes assigned to an orthogroup",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "genes": {
                "title": "Genes",
                "description": "Number of genes in this species' proteome",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Blues",
                "hidden": True,
            },
            "species_specific_orthogroups": {
                "title": "Species-Specific Orthogroups",
                "description": "Orthogroups containing genes from this species only",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Purples",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.per_species_data, headers)

    def _add_per_species_section(self) -> None:
        cats = {
            "genes_in_orthogroups": {"name": "Genes in Orthogroups", "color": "#7cb5ec"},
            "unassigned_genes": {"name": "Unassigned Genes", "color": "#f7a35c"},
        }
        pconfig = {
            "id": "orthofinder_per_species_plot",
            "title": "OrthoFinder: Genes Assigned to Orthogroups",
            "ylab": "Number of genes",
            "cpswitch_counts_label": "Number of genes",
        }
        self.add_section(
            name="Genes per Species",
            anchor="orthofinder-per-species",
            description="How each species' genes were distributed across orthogroups.",
            helptext="""
            A species with a large unassigned fraction is either more distant from
            the rest of the set than intended, or has a proteome of a different
            quality. Fragmented assemblies produce short predicted proteins that
            fail to cluster, so a single outlying species here often points at the
            input rather than at the clustering.
            """,
            plot=bargraph.plot(self.per_species_data, cats, pconfig),
        )

    def _add_overall_section(self) -> None:
        headers: dict[str, ColumnDict] = {
            "species": {
                "title": "Species",
                "description": "Number of species in the run",
                "min": 0,
                "scale": "Blues",
            },
            "orthogroups": {
                "title": "Orthogroups",
                "description": "Total number of orthogroups inferred",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Purples",
            },
            "percent_genes_in_orthogroups": {
                "title": "Genes in Orthogroups",
                "description": "Percentage of all genes assigned to an orthogroup",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "single_copy_orthogroups": {
                "title": "Single-Copy Orthogroups",
                "description": "Orthogroups with exactly one gene from every species",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Greens",
            },
            "orthogroups_all_species": {
                "title": "Orthogroups with All Species",
                "description": "Orthogroups containing at least one gene from every species",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Oranges",
            },
            "mean_orthogroup_size": {
                "title": "Mean Orthogroup Size",
                "description": "Mean number of genes per orthogroup",
                "min": 0,
                "format": "{:,.1f}",
                "scale": "Blues",
            },
            "median_orthogroup_size": {
                "title": "Median Orthogroup Size",
                "description": "Median number of genes per orthogroup",
                "min": 0,
                "format": "{:,.1f}",
                "scale": "Purples",
                "hidden": True,
            },
            "g50_assigned": {
                "title": "G50",
                "description": "Half of assigned genes are in orthogroups of at least this size",
                "min": 0,
                "scale": "Greens",
                "hidden": True,
            },
            "o50_assigned": {
                "title": "O50",
                "description": "Smallest number of orthogroups containing half of the assigned genes",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Oranges",
                "hidden": True,
            },
        }
        pconfig = {
            "id": "orthofinder_overall_table",
            "title": "OrthoFinder: Run Summary",
        }
        self.add_section(
            name="Run Summary",
            anchor="orthofinder-overall",
            description="Whole-run statistics, one row per OrthoFinder run.",
            helptext="""
            Single-copy orthogroups (one gene from every species) are the usual
            input for species-tree inference, so their count is a quick read on how
            much phylogenetic signal a run yields. G50 and O50 describe orthogroup
            size the way N50 describes contig length: half of all assigned genes sit
            in orthogroups of at least G50 genes, and the largest O50 orthogroups
            together hold half of the assigned genes. Median size, G50 and O50 are
            hidden by default and can be enabled from the column toggles.
            """,
            plot=table.plot(data=self.overall_data, headers=headers, pconfig=pconfig),
        )
