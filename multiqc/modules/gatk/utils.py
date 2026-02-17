def parse_report(lines, table_names):
    """Parse a GATK report https://software.broadinstitute.org/gatk/documentation/article.php?id=1244

    Only GATTable entries are parsed.  Tables are returned as a dict of tables.
    Each table is a dict of arrays, where names correspond to column names, and arrays
    correspond to column values.

    Args:
        lines (file handle): an iterable over the lines of a GATK report.
        table_names (dict): a dict with keys that are GATK report table names
            (e.g. "#:GATKTable:Quantized:Quality quantization map"), and values that are the
            keys in the returned dict.

    Returns:
        {
            table_1:
                {
                    col_1: [ val_1, val_2, ... ]
                    col_2: [ val_1, val_2, ... ]
                    ...
                }
            table_2:
                ...
        }
    """

    report = dict()
    lines = (line for line in lines)
    for line in lines:
        line = line.rstrip()
        if line in table_names.keys():
            report[table_names[line]] = parse_gatk_report_table(lines)
    return report


def parse_gatk_report_table(lines):
    headers = next(lines).rstrip().split()
    table = {h: [] for h in headers}
    for line in lines:
        line = line.rstrip()

        # testing to see if we have reached the end of a table in a GATKReport
        if line == "":
            break

        for index, value in enumerate(line.split()):
            table[headers[index]].append(value)
    return table
