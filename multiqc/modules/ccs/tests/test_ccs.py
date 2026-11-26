import pytest

PARSABLE_LINES = [
    "",
    "ZMWs input :",
    "ZMWs input (A) :",
    "ZMWs input : 93",
    "ZMWs input (A) : 93",
    "Coefficient of correlation    : 28.78%",
    "ZMWs generating CCS (B)  : 44 (47.31%)",
    "Coefficient of correlation  (A)  : 28.78%",
    "Below min length : 0 (-nan%)",
]

PARSED_RESULTS = [
    {},
    {"name": "ZMWs input"},
    {"name": "ZMWs input", "annotation": "A"},
    {"name": "ZMWs input", "count": 93},
    {"name": "ZMWs input", "annotation": "A", "count": 93},
    {"name": "Coefficient of correlation", "percentage": 28.78},
    {"name": "ZMWs generating CCS", "annotation": "B", "count": 44, "percentage": 47.31},
    {"name": "Coefficient of correlation", "percentage": 28.78, "annotation": "A"},
    {"name": "Below min length", "percentage": 0.0, "count": 0},
]

MARK = zip(PARSABLE_LINES, PARSED_RESULTS)

CCS_LINES = [
    [""],
    ["ZMWs input (A) : 93"],
    ["ZMWs input (A)  : 93", "Exclusive ZMW counts for (A):"],
    ["ZMWs filtered (C)  : 49 (52.69%)", "ZMWs input (A)  : 93", "Exclusive ZMW counts for (C):"],
    ["ZMWs input (A)  : 93", "ZMWs filtered (C)  : 49 (52.69%)", "Exclusive ZMW counts for (C):"],
    [
        "ZMWs input                (A) : 44",
        "ZMWs above all thresholds (B) : 39 (89%)",
        "ZMWs below any threshold  (C) : 5 (11%)",
        "",
        "ZMWs for (B):",
        "With same pair                : 39 (100%)",
        "Coefficient of correlation    : 28.78%",
        "",
        "ZMWs for (A):",
        "Allow diff pair               : 44 (100%)",
        "Allow same pair               : 44 (100%)",
        "",
        "Reads for (B):",
        "Above length                  : 39 (100%)",
        "Below length                  : 0 (0%)",
    ],
]

CCS_DATA = [
    {},
    {"ZMWs input": {"count": 93}},
    {"ZMWs input": {"count": 93, "Exclusive ZMW counts": {}}},
    {"ZMWs input": {"count": 93}, "ZMWs filtered": {"count": 49, "percentage": 52.69, "Exclusive ZMW counts": {}}},
    {"ZMWs input": {"count": 93}, "ZMWs filtered": {"count": 49, "percentage": 52.69, "Exclusive ZMW counts": {}}},
    {
        "ZMWs input": {
            "count": 44,
            "ZMWs": {
                "Allow diff pair": {"count": 44, "percentage": 100.0},
                "Allow same pair": {"count": 44, "percentage": 100.0},
            },
        },
        "ZMWs above all thresholds": {
            "count": 39,
            "percentage": 89.0,
            "ZMWs": {
                "With same pair": {"count": 39, "percentage": 100.0},
                "Coefficient of correlation": {"percentage": 28.78},
            },
            "Reads": {
                "Above length": {"count": 39, "percentage": 100.0},
                "Below length": {"count": 0, "percentage": 0.0},
            },
        },
        "ZMWs below any threshold": {"count": 5, "percentage": 11.0},
    },
]

CCS = zip(CCS_LINES, CCS_DATA)


@pytest.mark.parametrize(["line", "data"], MARK)
def test_parsable_lines(line, data):
    from multiqc.modules.ccs.ccs import parse_line

    parsed_line = parse_line(line)
    assert parsed_line == data


@pytest.mark.parametrize(["lines", "data"], CCS)
def test_parse_pacbio_log_ccs(lines, data):
    from multiqc.modules.ccs.ccs import parse_PacBio_log

    result = parse_PacBio_log(lines)
    assert result == data


def test_split_attributes():
    from multiqc.modules.ccs.ccs import split_attributes

    attributes = [
        {"id": "filtered_poor_snr", "name": "Below SNR threshold", "value": 1},
        {"id": "ccs_processing.zmw_filtered_too_few_passes", "name": "Lacking full passes", "value": 2},
        {"id": "ccs_processing.ds_reads_passed_yield", "name": "DS reads: Passed", "value": 3},
        {"id": "ccs_processing.ss_reads_input", "name": "SS reads: Inputs", "value": 4},
        {"id": "ccs_processing.number_of_ccs_reads_q20", "name": "HiFi Reads", "value": 5},
    ]
    groups, other = split_attributes(attributes)
    assert {g: {k: a["value"] for k, a in attrs.items()} for g, attrs in groups.items()} == {
        "zmw": {"filtered_poor_snr": 1, "filtered_too_few_passes": 2},
        "ds_reads": {"passed_yield": 3},
        "ss_reads": {"input": 4},
    }
    assert {k: a["value"] for k, a in other.items()} == {"number_of_ccs_reads_q20": 5}
