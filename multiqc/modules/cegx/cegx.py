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

        return data
