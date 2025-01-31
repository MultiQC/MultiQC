import json

import multiqc
from multiqc.core.update_config import ClConfig


def test_rerun_json(data_dir, tmp_path):
    """Test running MultiQC on intermediate output multiqc_data.json.
    Run MultiQC on a set of inputs.
    Then run MultiQC on the multiqc_data.json from the first run. The reports should be identical.
    It should only work with multiqc_data.json is passed explicitly. Otherwise MultiQC output folder should be ignored.
    """
    # Run 1: Normal run on test data
    run_a_dir = tmp_path / "run_a"
    run_a_dir.mkdir()
    multiqc.run(data_dir / "modules/fastp/SAMPLE.json", cfg=ClConfig(output_dir=run_a_dir))

    # Get the first report contents
    with open(run_a_dir / "multiqc_data" / "multiqc_data.json") as f:
        report1_data = json.load(f)

    # Run 2: Run on the intermediate data from run1
    run_b_dir = tmp_path / "run_b"
    run_b_dir.mkdir()
    multiqc.run(run_a_dir / "multiqc_data" / "multiqc_data.json", cfg=ClConfig(output_dir=run_b_dir))

    # Compare reports
    with open(run_b_dir / "multiqc_data" / "multiqc_data.json") as f:
        report2_data = json.load(f)

    assert report1_data == report2_data


def test_rerun_and_combine(data_dir, tmp_path):
    """Test adding new data to report.
    Run MultiQC on a set of inputs.
    Then run MultiQC on the intermediate output multiqc_data.json of that run, plus add new inputs.
    The output should be the same as if MultiQC was run on the combined inputs directly.
    """
    # Run MultiQC on inputs A
    run_a_dir = tmp_path / "run_a"
    multiqc.run(data_dir / "modules/fastp/SAMPLE.json", cfg=ClConfig(output_dir=run_a_dir))

    # Run MultiQC on outputs of A + inputs B
    run_ab_dir = tmp_path / "run_ab"
    multiqc.run(
        run_a_dir / "multiqc_data" / "multiqc_data.json",
        data_dir / "modules/fastp/single_end",
        cfg=ClConfig(output_dir=run_ab_dir),
    )

    with open(run_ab_dir / "multiqc_data" / "multiqc_data.json") as f:
        report_ab_data = json.load(f)

    # Run MultiQC on inputs A + B directly
    run_ab2_dir = tmp_path / "run_ab2"
    multiqc.run(
        data_dir / "modules/fastp/SAMPLE.json",
        data_dir / "modules/fastp/single_end",
        cfg=ClConfig(output_dir=run_ab2_dir),
    )

    with open(run_ab2_dir / "multiqc_data" / "multiqc_data.json") as f:
        report_ab2_data = json.load(f)

    assert report_ab_data == report_ab2_data
