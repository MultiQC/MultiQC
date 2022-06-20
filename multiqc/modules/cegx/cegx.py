import logging
from multiqc.modules.base_module import BaseMultiqcModule
import re

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
          name='CEGX',
          anchor='cegx',
          href="https://cambridge-epigenetix.com",
          info="sample summary statistics from CEGX pipeline", #FIXME
          doi="" #FIXME
        )

        for f in self.find_log_files('cegx/mbias', filehandles=True):
            # f['f'] is now a filehandle instead of contents
            # for l in f['f']:
            #     print(l)
            self.parse_cegx_mbias(f['f'])

    def parse_cegx_mbias(self, file_content):
        """Parse cegx mbias log file"""
        data = dict()
        header = ['pos', 'modC', 'unmodC', 'perc_mod', 'coverage']


        for line in file_content:
            # Get rid of trailing newlines
            line = line.strip("\t")
            # file does not contain a header
            # pos  modC   unmodC  % modC  coverage
            #  1   2367	 44883	 0.0501	  47250
            #  2   1712	 46021	 0.0359	  47733

            linedata = self.parse_line(line, header)
            print(tuple(linedata))



        return data

    def parse_line(self, line, headers): #for tsv that has no defined headers
        """Parse a line from the cegx log file"""
        data = dict()

        # If we got an empty line to parse
        if not line.strip():
            return data

        # convert data into a dict applying header definition
        data = zip(headers, line)

        # # If the percentage is this weird nan, set it to 0
        # if percentage == "-nan%":
        #     data["perc_mod"] = 0.0
        # # Add the count as an integer
        # data["count"] = int(count)
        #
        # # If there is a single value
        # if len(values) == 1:
        #     value = values.pop()
        #     # Is the value a percentage
        #     if value.endswith("%"):
        #         data["percentage"] = float(value[:-1])
        #     # Otherwise, it is an integer
        #     else:
        #         data["count"] = int(value)
        # elif len(values) == 2:
        #     # If there are multiple values, they are in the format
        #     # count (percentage%)
        #     count = values[0]
        #     percentage = values[1]
        #
        #     # The percentage should be in the format:
        #     # (12.34%) or (100%) or (0%) or (-nan%)
        #     # So we can remove the bracets first
        #     percentage = percentage[1:-1]
        #     # Then we make sure it is one of the supported patterns
        #     assert re.fullmatch(r"(\d+\.\d+%|\d+%|-nan%)", percentage)
        #
        #

        return data
