import json
import os
from collections import Counter

import requests


def main():
    latest_version = get_latest_version()
    print(f"Latest Version: {latest_version}")

    pypi_stats = fetch_pypi_stats()

    github_token = os.environ.get("GITHUB_TOKEN")
    github_clones_stats = fetch_github_clones_stats(github_token)
    github_releases_stats = fetch_github_releases_stats(github_token)

    print(f"PyPI downloads: {pypi_stats.get('total_downloads')}")
    print(f"GutHub clones: {github_clones_stats.get('total_clones')}")
    print(f"GitHub downloads: {github_releases_stats.get('total_downloads')}")

    dockerhub_stats = fetch_dockerhub_stats()
    biocontainers_stats = fetch_biocontainers_stats()

    print(f"DockerHub pulls: {dockerhub_stats.get('total_pulls')}")
    print(f"BioContainers pulls: {biocontainers_stats.get('total_pulls')}")


def get_latest_version() -> str:
    pypi_url = "https://pypi.org/pypi/multiqc/json"
    pypi_response = requests.get(pypi_url)
    if pypi_response.status_code != 200:
        raise RuntimeError("Failed to fetch data from PyPI")
    pypi_data = json.loads(pypi_response.text)
    latest_version = pypi_data["info"]["version"]
    return latest_version


def fetch_pypi_stats():
    pepy_url = "https://api.pepy.tech/api/v2/projects/multiqc"

    total_downloads = None
    downloads_by_version = Counter()

    # Fetch download count from pepy.tech
    response = requests.get(pepy_url)
    if response.status_code == 200:
        data = json.loads(response.text)
        total_downloads = data["total_downloads"]
        for date, downloads_by_version in data["downloads"].items():
            for version, downloads in downloads_by_version.items():
                downloads_by_version[version] += downloads
    else:
        print("Failed to fetch data from pepy.tech")

    return {
        "total_downloads": total_downloads,
        "downloads_by_version": downloads_by_version,
    }


def fetch_github_clones_stats(github_token=os.environ.get("GITHUB_TOKEN")):
    url = "https://api.github.com/repos/ewels/MultiQC/traffic/clones"
    headers = {"Authorization": f"token {github_token}"}

    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        github_data = json.loads(response.text)
        clone_count = github_data.get("count", 0)
    else:
        print(f"Failed to fetch data from GitHub, status code: {response.status_code}, url: {url}")
        clone_count = None

    return {"total_clones": clone_count}


def fetch_github_releases_stats(github_token):
    url = "https://api.github.com/repos/ewels/MultiQC/releases"
    headers = {"Authorization": f"token {github_token}"}
    total_downloads = 0

    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = json.loads(response.text)
        for release in data:
            for asset in release["assets"]:
                total_downloads += asset["download_count"]
    else:
        print(f"Failed to fetch release data from GitHub, status code: {response.status_code}, url: {url}")

    return {"total_downloads": total_downloads}


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


def fetch_dockerhub_stats():
    return {"total_pulls": _fetch_dockerhub_count("ewels/multiqc/")}


def fetch_biocontainers_stats():
    dockerhub = _fetch_dockerhub_count("biocontainers/multiqc/")
    quay = _fetch_quay_count("biocontainers/multiqc/")
    if dockerhub is not None and quay is not None:
        total = dockerhub + quay
    else:
        total = None
    return {"total_pulls": total}


def fetch_galaxy_stats():
    # Fetching statistics from Galaxy can be a bit challenging because it doesn't
    # offer a simple API endpoint for getting tool usage statistics directly.
    # Moreover, Galaxy instances can be independent from each other, and the tool
    # might be hosted on multiple instances.
    pass


if __name__ == "__main__":
    main()
