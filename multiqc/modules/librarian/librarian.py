""" MultiQC module to parse output from Librarian """

import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import heatmap
from multiqc.utils import mqc_colour

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(
            name="Librarian",
            anchor="librarian",
            href="https://github.com/DesmondWillowbrook/Librarian",
            info=" - a tool to predict the sequencing library type from the base composition of a FastQ file.",
            doi="10.12688/f1000research.125325.1",
        )

        # To store the summary data
        self.librarian_data = dict()

        # Parse the output files
        self.parse_librarian_data()

        # Remove filtered samples
        self.librarian_data = self.ignore_samples(self.librarian_data)

        # Write the data files to disk
        if not self.librarian_data:
            raise ModuleNoSamplesFound
        log.info("Found {} samples".format(len(self.librarian_data)))

        if self.librarian_data:
            self.write_data_file(self.librarian_data, "multiqc_librarian_data")

        self.plot_librarian_heatmap()
        show_general_stats = getattr(config, "librarian", {}).get("show_general_stats", False)
        if show_general_stats:
            self.add_general_stats()

    def parse_librarian_data(self):
        """Loop through Librarian files and parse their data"""
        for f in self.find_log_files("librarian", filehandles=True):
            headers = f["f"].readline().strip().split("\t")
            for l in f["f"]:
                data = dict(zip(headers, l.strip().split("\t")))
                s_name = self.clean_s_name(data["sample_name"], f)
                del data["sample_name"]
                data_float = {}
                for k, v in data.items():
                    try:
                        data_float[k] = float(v)
                    except:
                        data_float[k] = v
                if s_name in self.librarian_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.librarian_data[s_name] = data_float
                self.add_data_source(f, s_name=s_name)

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, s_name)

    def plot_librarian_heatmap(self):
        """Make the heatmap plot"""
        # Get all library types
        hm_library_types = set()
        for d in self.librarian_data.values():
            hm_library_types.update(list(d.keys()))
        hm_library_types = sorted(list(hm_library_types))
        hm_data = []
        hm_sample_names = []
        for s_name in self.librarian_data:
            hm_sample_names.append(s_name)
            sample_data = []
            for library_type in hm_library_types:
                try:
                    sample_data.append(self.librarian_data[s_name][library_type])
                except KeyError:
                    sample_data.append(None)
            hm_data.append(sample_data)

        pconfig = {
            "title": "Librarian: Library Predictions",
            "xTitle": "Library type",
            "yTitle": "Sample name",
            "min": 0,
            "max": 100,
            "square": False,
            "xcats_samples": False,
            "ycats_samples": True,
            "colstops": [[0, "#FFFFFF"], [1, "#FF0000"]],
            "decimalPlaces": 0,
        }

        self.add_section(
            name="Library Type Prediction",
            anchor="librarian-library-type",
            description="""
                For each projected test library, the location on the Compositions/Probability Map is determined.
                This plot shows how published library types are represented at the same location.
            """,
            helptext="""
                Some regions on the map are very specific to a certain library type,
                others are more mixed. Therefore, for some test libraries the results
                will be much clearer than for others.

                The different plots are intended to provide a good overview of how similar
                the test library is to published data.

                The cause of any deviations should be inspected;
                the interpretation will be different depending on how characteristic
                the composition signature of the library type and how far off the
                projection of the test sample is.
            """,
            plot=heatmap.plot(hm_data, hm_library_types, hm_sample_names, pconfig),
        )

    def add_general_stats(self):
        """Add general stats column for most likely library type"""
        # Get most likely library type for each sample
        data = dict()
        lib_types = set()
        for s_name, d in self.librarian_data.items():
            lib_type = max(d, key=d.get)
            score = float(d[lib_type])
            data[s_name] = {
                "most_likely_library_type": lib_type,
                "most_likely_library_type_score": score,
            }
            lib_types.add(lib_type)

        # Generate colours for library types
        lib_types = list(lib_types)
        lib_type_colours = mqc_colour.mqc_colour_scale("Paired", 0, len(lib_types))
        bgcols = {}
        for idx, lib_type in enumerate(lib_types):
            bgcols[lib_type] = lib_type_colours.get_colour(idx)

        headers = {
            "most_likely_library_type": {
                "title": "Likely Type",
                "description": "Most likely library type.",
                "bgcols": bgcols,
            },
            "most_likely_library_type_score": {
                "title": "Type score",
                "description": "Library prediction type score",
                "format": "{:,.0f}",
                "scale": "RdYlGn",
                "min": 0,
                "max": 100,
            },
        }
        self.general_stats_addcols(data, headers)
