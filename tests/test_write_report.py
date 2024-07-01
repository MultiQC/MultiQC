import multiqc
from multiqc import report
from multiqc.core.update_config import ClConfig


def test_write_stdout(data_dir, capsys):
    """
    Verify stdout option: stdout contains only HTML, nothing else is written to disk
    """
    report.reset()
    mod_dir = data_dir / "modules" / "kallisto"
    assert mod_dir.exists() and mod_dir.is_dir()

    multiqc.run(mod_dir, cfg=ClConfig(filename="stdout"))
    captured = capsys.readouterr()

    assert captured.out.startswith("""<!DOCTYPE html>\n<html lang="en">\n<head>""")
