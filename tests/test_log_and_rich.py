import logging

from multiqc.core.log_and_rich import _InfoFormatter


def test_multiline_message_keeps_the_module_column():
    record = logging.makeLogRecord({"msg": "first\n\nthird", "module": "ai", "levelno": logging.ERROR})
    assert _InfoFormatter("%(module)18s | %(message)s").format(record).split("\n") == [
        f"{'ai':>18} | first",
        f"{'':>18} | ",
        f"{'':>18} | third",
    ]
