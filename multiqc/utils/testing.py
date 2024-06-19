from multiqc import config


def data_dir():
    repo_path = config.REPO_DIR / "test-data"
    if not repo_path.exists():
        raise FileNotFoundError(
            f"The test data directory expected to be found at {repo_path}. Please, "
            f"clone the repository with the test data by changing to the MultiQC repo root, "
            f"and running: `git clone https://github.com/MultiQC/test-data`"
        )

    return repo_path / "data"
