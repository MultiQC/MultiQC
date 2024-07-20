import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses _only_ `*_jf.hist` files. The general usage of jellyfish to be parsed by MultiQC module needs to be:

    - `gunzip -c file.fastq.gz | jellyfish count -o file.jf -m ...`
    - `jellyfish histo -o file_jf.hist -f file.jf`

    In case a user wants to customise the matching pattern for jellyfish, then multiqc can be run with the option `--cl-config "sp: { jellyfish: { fn: 'PATTERN' } }"` where `PATTERN` is the pattern to be matched. For example:

    ```bash
    multiqc . --cl-config "sp: { jellyfish: { fn: '*.hist' } }"
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Jellyfish",
            anchor="jellyfish",
            href="https://github.com/gmarcais/Jellyfish",
            info="Counting k-mers in DNA.",
            extra="""
            A k-mer is a substring of length k, and counting the occurrences of all such substrings 
            is a central step in many analyses of DNA sequence. JELLYFISH can count k-mers using an order of 
            magnitude less memory and an order of magnitude faster than other k-mer counting packages by using 
            an efficient encoding of a hash table and by exploiting the "compare-and-swap" CPU instruction to 
            increase parallelism.
            """,
            doi="10.1093/bioinformatics/btr011",
        )

        self.jellyfish_data = dict()
        self.jellyfish_max_x = 0
        for f in self.find_log_files("jellyfish", filehandles=True):
            self.parse_jellyfish_data(f)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        if self.jellyfish_max_x < 100:
            # the maximum is below 100, we display anyway up to 200
            self.jellyfish_max_x = 200
        else:
            # in this case the area plotted is a function of the maximum x
            self.jellyfish_max_x = 2 * self.max_key

        self.jellyfish_data = self.ignore_samples(self.jellyfish_data)

        if len(self.jellyfish_data) == 0:
            raise ModuleNoSamplesFound

        # Write data to file
        self.write_data_file(self.jellyfish_data, "jellyfish")

        log.info(f"Found {len(self.jellyfish_data)} reports")

        self.frequencies_plot(xmax=self.jellyfish_max_x)

    def parse_jellyfish_data(self, f):
        """Go through the hist file and memorise it"""
        histogram = {}
        occurence = 0
        for line in f["f"]:
            line = line.rstrip("\n")
            occurence = int(line.split(" ")[0])
            count = int(line.split(" ")[1])
            histogram[occurence] = occurence * count
        # delete last occurnece as it is the sum of all kmer occuring more often than it.
        del histogram[occurence]
        # sanity check
        self.max_key = max(histogram, key=histogram.get)
        self.jellyfish_max_x = max(self.jellyfish_max_x, self.max_key)
        if len(histogram) > 0:
            if f["s_name"] in self.jellyfish_data:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.add_data_source(f)
            self.jellyfish_data[f["s_name"]] = histogram

    def frequencies_plot(self, xmin=0, xmax=200):
        """Generate the qualities plot"""

        helptext = """
            A possible way to assess the complexity of a library even in
            absence of a reference sequence is to look at the kmer profile of the reads.
            The idea is to count all the kmers (_i.e._, sequence of length `k`) that occur
            in the reads. In this way it is possible to know how many kmers occur
            `1,2,.., N` times and represent this as a plot.
            This plot tell us for each x, how many k-mers (y-axis) are present in the
            dataset in exactly x-copies.

            In an ideal world (no errors in sequencing, no bias, no  repeated regions)
            this plot should be as close as  possible to a gaussian distribution.
            In reality we will always see a peak for `x=1` (_i.e._, the errors)
            and another peak close to the expected coverage. If the genome is highly
            heterozygous a second peak at half of the coverage can be expected."""

        pconfig = {
            "id": "Jellyfish_kmer_plot",
            "title": "Jellyfish: K-mer plot",
            "ylab": "Counts",
            "xlab": "k-mer frequency",
            "x_decimals": False,
            "xmin": xmin,
            "xmax": xmax,
        }

        self.add_section(
            anchor="jellyfish_kmer_plot",
            description="The K-mer plot lets you estimate library complexity and coverage from k-mer content.",
            helptext=helptext,
            plot=linegraph.plot(self.jellyfish_data, pconfig),
        )
