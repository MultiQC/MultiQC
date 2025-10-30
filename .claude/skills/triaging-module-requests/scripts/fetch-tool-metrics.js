#!/usr/bin/env node

/**
 * Tool Metrics Fetcher
 *
 * Fetches popularity and usage metrics for bioinformatics tools from various sources:
 * - GitHub: stars, forks, last commit date
 * - PyPI: monthly downloads
 * - Conda: total downloads (via anaconda.org)
 *
 * Usage:
 *   node fetch-tool-metrics.js github owner/repo
 *   node fetch-tool-metrics.js pypi package-name
 *   node fetch-tool-metrics.js conda package-name
 */

const https = require("https");
const { execSync } = require("child_process");

/**
 * Make HTTPS request and return parsed JSON
 */
function httpsGet(url, headers = {}) {
  return new Promise((resolve, reject) => {
    const options = {
      headers: {
        "User-Agent": "MultiQC-Module-Triage",
        ...headers,
      },
    };

    https
      .get(url, options, (res) => {
        let data = "";
        res.on("data", (chunk) => (data += chunk));
        res.on("end", () => {
          try {
            resolve({ status: res.statusCode, data: JSON.parse(data) });
          } catch (e) {
            reject(new Error(`Failed to parse JSON: ${e.message}`));
          }
        });
      })
      .on("error", reject);
  });
}

/**
 * Fetch GitHub repository metrics
 */
async function fetchGitHubMetrics(repoPath) {
  try {
    // Try using gh CLI first (respects authentication)
    try {
      const result = execSync(`gh api repos/${repoPath}`, { encoding: "utf-8" });
      const data = JSON.parse(result);

      return {
        source: "github",
        repository: repoPath,
        stars: data.stargazers_count,
        forks: data.forks_count,
        lastCommit: data.pushed_at,
        archived: data.archived,
        url: data.html_url,
      };
    } catch (ghError) {
      // Fallback to public API if gh CLI fails
      console.error("gh CLI failed, trying public API...");
      const { status, data } = await httpsGet(`https://api.github.com/repos/${repoPath}`);

      if (status !== 200) {
        throw new Error(`GitHub API returned ${status}`);
      }

      return {
        source: "github",
        repository: repoPath,
        stars: data.stargazers_count,
        forks: data.forks_count,
        lastCommit: data.pushed_at,
        archived: data.archived,
        url: data.html_url,
      };
    }
  } catch (error) {
    throw new Error(`Failed to fetch GitHub metrics: ${error.message}`);
  }
}

/**
 * Fetch PyPI package metrics
 */
async function fetchPyPIMetrics(packageName) {
  try {
    // Get package info
    const { status, data } = await httpsGet(`https://pypi.org/pypi/${packageName}/json`);

    if (status !== 200) {
      throw new Error(`PyPI API returned ${status}`);
    }

    // Try to get download stats from pypistats.org
    let downloads = null;
    try {
      const statsUrl = `https://pypistats.org/api/packages/${packageName}/recent`;
      const statsResponse = await httpsGet(statsUrl);
      if (statsResponse.status === 200) {
        downloads = {
          lastMonth: statsResponse.data.data.last_month,
          lastWeek: statsResponse.data.data.last_week,
          lastDay: statsResponse.data.data.last_day,
        };
      }
    } catch (e) {
      // Download stats not available
      console.error("PyPI download stats not available");
    }

    return {
      source: "pypi",
      package: packageName,
      downloads: downloads,
      url: `https://pypi.org/project/${packageName}/`,
    };
  } catch (error) {
    throw new Error(`Failed to fetch PyPI metrics: ${error.message}`);
  }
}

/**
 * Fetch Conda package metrics
 */
async function fetchCondaMetrics(packageName) {
  try {
    // Try Anaconda.org API
    const { status, data } = await httpsGet(`https://api.anaconda.org/package/bioconda/${packageName}`);

    if (status !== 200) {
      throw new Error(`Anaconda API returned ${status}`);
    }

    // Sum downloads from all files
    const totalDownloads = data.files ? data.files.reduce((sum, file) => sum + (file.ndownloads || 0), 0) : 0;

    return {
      source: "conda",
      package: packageName,
      downloads: totalDownloads,
      url: `https://anaconda.org/bioconda/${packageName}`,
    };
  } catch (error) {
    throw new Error(`Failed to fetch Conda metrics: ${error.message}`);
  }
}

/**
 * Main function
 */
async function main() {
  const args = process.argv.slice(2);

  if (args.length < 2) {
    console.error("Usage: fetch-tool-metrics.js <source> <identifier>");
    console.error("Sources: github, pypi, conda");
    console.error("Examples:");
    console.error("  fetch-tool-metrics.js github nextflow-io/nextflow");
    console.error("  fetch-tool-metrics.js pypi multiqc");
    console.error("  fetch-tool-metrics.js conda fastqc");
    process.exit(1);
  }

  const [source, identifier] = args;

  try {
    let metrics;

    switch (source.toLowerCase()) {
      case "github":
        metrics = await fetchGitHubMetrics(identifier);
        break;
      case "pypi":
        metrics = await fetchPyPIMetrics(identifier);
        break;
      case "conda":
        metrics = await fetchCondaMetrics(identifier);
        break;
      default:
        throw new Error(`Unknown source: ${source}. Use github, pypi, or conda`);
    }

    console.log(JSON.stringify(metrics, null, 2));
  } catch (error) {
    console.error(`Error: ${error.message}`);
    process.exit(1);
  }
}

// Run if called directly
if (require.main === module) {
  main();
}

module.exports = {
  fetchGitHubMetrics,
  fetchPyPIMetrics,
  fetchCondaMetrics,
};
