import json
import os
from collections import Counter
from pprint import pprint

import requests

from multiqc.utils.software_versions import parse_version


def main():
    stats = {}
    stats |= pypi_stats()
    stats |= bioconda_stats()
    stats |= github_clones_stats()
    stats |= github_releases_stats()
    stats |= dockerhub_stats()
    # stats |= biocontainers_stats()
    pprint(stats)


def get_latest_version() -> str:
    pypi_url = "https://pypi.org/pypi/multiqc/json"
    pypi_response = requests.get(pypi_url)
    if pypi_response.status_code != 200:
        raise RuntimeError("Failed to fetch data from PyPI")
    pypi_data = json.loads(pypi_response.text)
    latest_version = pypi_data["info"]["version"]
    return latest_version


def pypi_stats():
    pepy_url = "https://api.pepy.tech/api/v2/projects/multiqc"

    total_downloads = None
    downloads_by_version = Counter()

    # Fetch download count from pepy.tech
    response = requests.get(pepy_url)
    if response.status_code == 200:
        data = json.loads(response.text)
        total_downloads = data["total_downloads"]
        for date, date_downloads_by_version in data["downloads"].items():
            for version, downloads in date_downloads_by_version.items():
                downloads_by_version[version] += downloads
    else:
        print("Failed to fetch data from pepy.tech")

    downloads_by_version = {parse_version(k): v for k, v in downloads_by_version.items()}
    return {
        "pypi_downloads": total_downloads,
        "pypi_downloads_by_version": downloads_by_version,
    }


def bioconda_stats():
    url = "https://raw.githubusercontent.com/bioconda/bioconda-plots/main/plots/multiqc/versions.json"
    # [{"date":"2023-09-05","total":15949,"delta":18,"version":"1.10.1"}, ...

    downloads_by_version = dict()

    response = requests.get(url)
    if response.status_code == 200:
        data = json.loads(response.text)
        downloads_by_version = Counter()
        for version in data:
            downloads_by_version[version["version"]] = version["total"]
    else:
        print("Failed to fetch data from bioconda")

    downloads_by_version = {parse_version(k): v for k, v in downloads_by_version.items()}
    return {
        "bioconda_downloads": sum(downloads_by_version.values()),
        "bioconda_downloads_by_version": downloads_by_version,
    }


def github_clones_stats():
    github_token = os.environ.get("GITHUB_TOKEN")
    url = "https://api.github.com/repos/ewels/MultiQC/traffic/clones"
    headers = {"Authorization": f"token {github_token}"}

    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        github_data = json.loads(response.text)
        clone_count = github_data.get("count", 0)
    else:
        print(f"Failed to fetch data from GitHub, status code: {response.status_code}, url: {url}")
        clone_count = None

    return {"github_clones": clone_count}


def github_releases_stats():
    github_token = os.environ.get("GITHUB_TOKEN")
    url = "https://api.github.com/repos/ewels/MultiQC/releases"
    headers = {"Authorization": f"token {github_token}"}

    downloads_by_version = Counter()

    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = json.loads(response.text)
        for release in data:
            version = parse_version(release["tag_name"])
            for asset in release["assets"]:
                downloads_by_version[version] += asset["download_count"]
    else:
        print(f"Failed to fetch release data from GitHub, status code: {response.status_code}, url: {url}")

    return {
        "github_releases_downloads": sum(downloads_by_version.values()),
        "github_releases_downloads_by_version": downloads_by_version,
    }


def _fetch_dockerhub_count(repo):
    url = f"https://hub.docker.com/v2/repositories/{repo}"
    response = requests.get(url)
    if response.status_code == 200:
        data = json.loads(response.text)
        return data.get("pull_count", 0)
    else:
        print(f"Failed to fetch data from DockerHub, status code: {response.status_code}, url: {url}")
        return None


def _fetch_quay_count(repo):
    url = f"https://quay.io/api/v1/repository/{repo}"
    response = requests.get(url)
    if response.status_code == 200:
        data = json.loads(response.text)
        return data.get("pull_count", 0)
    else:
        print(f"Failed to fetch data from Quay.io, status code: {response.status_code}, url: {url}")
        return None


def dockerhub_stats():
    return {"dockerhub_pulls": _fetch_dockerhub_count("ewels/multiqc/")}


def biocontainers_stats():
    dockerhub = _fetch_dockerhub_count("biocontainers/multiqc/")
    quay = _fetch_quay_count("biocontainers/multiqc/")
    if dockerhub is not None and quay is not None:
        total = dockerhub + quay
    else:
        total = None
    return {"biocontainers_pulls": total}


if __name__ == "__main__":
    main()
