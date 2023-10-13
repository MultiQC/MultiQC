def read_histogram(self, program_key, headers, formats):
    """
    Reads a Picard HISTOGRAM file.

    Args:
        self: the Picard QC module
        program_key: the key used to find the program (ex. picard/quality_by_cycle)
        headers: the list of expected headers for the histogram
        formats: the list of methods to apply to re-format each field (on a given row)
    """
    all_data = dict()

    assert len(formats) == len(headers)

    # Go through logs and find Metrics
    for f in self.find_log_files(program_key, filehandles=True):
        s_name = f["s_name"]
        lines = iter(f["f"])
        for l in lines:
            maybe_s_name = self.extract_sample_name(l, f)
            if maybe_s_name:
                s_name = maybe_s_name
            if l.startswith("## HISTOGRAM"):
                break

        sample_data = dict()

        try:
            # skip to the histogram
            line = next(lines)
            while not line.startswith("## HISTOGRAM"):
                line = next(lines)

            # check the header
            line = next(lines)
            if headers != line.strip().split("\t"):
                continue

            # slurp the data
            line = next(lines).rstrip()
            while line:
                fields = line.split("\t")
                assert len(fields) == len(headers)
                for i in range(len(fields)):
                    fields[i] = formats[i](fields[i])

                sample_data[fields[0]] = dict(zip(headers, fields))
                line = next(lines).rstrip()

        except StopIteration:
            pass

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
