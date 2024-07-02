import os

import multiqc
from multiqc import report
from multiqc.core.update_config import ClConfig


def test_write_default(data_dir, tmp_path):
    """
    Verify HTML and data directory with default names are written to current dir
    """
    report.reset()
    mod_dir = data_dir / "modules" / "kallisto"
    assert mod_dir.exists() and mod_dir.is_dir()

    files_before = set(os.listdir(tmp_path))
    os.chdir(tmp_path)

    multiqc.run(mod_dir)

    files_after = set(os.listdir(tmp_path))
    assert files_after - files_before == {"multiqc_report.html", "multiqc_data"}


def test_write_stdout(data_dir, tmp_path, capsys):
    """
    Verify stdout option: stdout contains only HTML, nothing else is written to disk
    """
    report.reset()
    mod_dir = data_dir / "modules" / "kallisto"
    assert mod_dir.exists() and mod_dir.is_dir()

    files_before = set(os.listdir(tmp_path))
    os.chdir(tmp_path)

    multiqc.run(mod_dir, cfg=ClConfig(filename="stdout"))

    captured = capsys.readouterr()
    assert captured.out.startswith("""<!DOCTYPE html>\n<html lang="en">\n<head>""")

    files_after = set(os.listdir(tmp_path))
    assert files_before == files_after


def test_custom_content_with_config(tmp_path, capsys):
    file = tmp_path / "mysample-concordance.txt"
    file.write_text("""Sample	'08021342'	'08027127'\n'08021342'	1.0	0.378""")

    conf_file = tmp_path / "multiqc_config.yaml"
    conf_file.write_text(
        """\
custom_data:
    concordance:
        id: 'concordance'
        section_name: 'Concordance Rates'
        plot_type: 'heatmap'
        pconfig:
            id: 'concordance_heatmap'
        sort_rows: true
sp:
    concordance:
        fn: '*concordance.txt'
"""
    )

    multiqc.run(  # pylint: disable=no-member
        file,
        cfg=ClConfig(
            strict=True,
            force=True,
            config_files=[conf_file],
            filename="stdout",
        ),
    )

    out = capsys.readouterr().out
    assert '<h2 class="mqc-module-title" id="concordance-module">Concordance Rates</h2>' in out
    assert '<div class="mqc-section mqc-section-concordance-module">' in out
    assert 'value="0.378"' in out
