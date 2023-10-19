def read_histogram(self, program_key, headers, formats, picard_tool, sentieon_algo=None):
    """
    Reads a Picard HISTOGRAM file.

    Args:
        self: the Picard QC module
        program_key: the key used to find the program (ex. picard/quality_by_cycle)
        headers: the list of expected headers for the histogram
        formats: the list of methods to apply to re-format each field (on a given row)
        picard_tool: the name of the Picard tool to be found in the header, e.g. MeanQualityByCycle
        sentieon_algo: the name of the Sentieon algorithm to be found in the header, e.g. MeanQualityByCycle
    """
    all_data = dict()

    assert len(formats) == len(headers)

    # Go through logs and find Metrics
    for f in self.find_log_files(program_key, filehandles=True):
        s_name = f["s_name"]
        lines = list(f["f"])
        for line in lines:
            maybe_s_name = self.extract_sample_name(line, f, picard_tool=picard_tool, sentieon_algo=sentieon_algo)
            if maybe_s_name:
                s_name = maybe_s_name
            if self.is_line_right_before_table(line, sentieon_algo=sentieon_algo):
                break

        sample_data = dict()

        # skip to the histogram
        if lines:
            line = lines.pop(0)
            while lines and not self.is_line_right_before_table(line, sentieon_algo=sentieon_algo):
                line = lines.pop(0)

        # check the header
        if lines:
            line = lines.pop(0)
            if headers != line.strip().split("\t"):
                continue

        # slurp the data
        if lines:
            line = lines.pop(0).rstrip()
            while line:
                fields = line.split("\t")
                assert len(fields) == len(headers)
                for i in range(len(fields)):
                    fields[i] = formats[i](fields[i])

                sample_data[fields[0]] = dict(zip(headers, fields))
                if lines:
                    line = lines.pop(0).rstrip()

        # append the data
        if sample_data:
            all_data[s_name] = sample_data

            self.add_data_source(f, s_name, section="Histogram")
            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

    data = self.ignore_samples(all_data)

    # Write data to file
    self.write_data_file(data, f"{self.anchor}_histogram")

    return data
