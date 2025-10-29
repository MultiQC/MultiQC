#!/usr/bin/env node

/**
 * Tool Metrics Fetcher
 *
 * Fetches popularity and usage metrics for bioinformatics tools from various sources:
 * - GitHub: stars, forks, last commit date, open issues
 * - PyPI: monthly downloads, release info
 * - Conda: package downloads (via anaconda.org)
 * - npm: weekly downloads
 *
 * Usage:
 *   node fetch-tool-metrics.js github owner/repo
 *   node fetch-tool-metrics.js pypi package-name
 *   node fetch-tool-metrics.js conda package-name
 *   node fetch-tool-metrics.js npm package-name
 */

const https = require('https');
const { execSync } = require('child_process');

/**
 * Make HTTPS request and return parsed JSON
 */
function httpsGet(url, headers = {}) {
  return new Promise((resolve, reject) => {
    const options = {
      headers: {
        'User-Agent': 'MultiQC-Module-Triage',
        ...headers
      }
    };

    https.get(url, options, (res) => {
      let data = '';
      res.on('data', (chunk) => data += chunk);
      res.on('end', () => {
        try {
          resolve({ status: res.statusCode, data: JSON.parse(data) });
        } catch (e) {
          reject(new Error(`Failed to parse JSON: ${e.message}`));
        }
      });
    }).on('error', reject);
  });
}

/**
 * Fetch GitHub repository metrics
 */
async function fetchGitHubMetrics(repoPath) {
  try {
    // Try using gh CLI first (respects authentication)
    try {
      const result = execSync(`gh api repos/${repoPath}`, { encoding: 'utf-8' });
      const data = JSON.parse(result);

      return {
        source: 'github',
        repository: repoPath,
        stars: data.stargazers_count,
        forks: data.forks_count,
        openIssues: data.open_issues_count,
        watchers: data.subscribers_count,
        lastCommit: data.pushed_at,
        createdAt: data.created_at,
        language: data.language,
        archived: data.archived,
        license: data.license?.spdx_id,
        defaultBranch: data.default_branch,
        hasWiki: data.has_wiki,
        hasPages: data.has_pages,
        url: data.html_url
      };
    } catch (ghError) {
      // Fallback to public API if gh CLI fails
      console.error('gh CLI failed, trying public API...');
      const { status, data } = await httpsGet(`https://api.github.com/repos/${repoPath}`);

      if (status !== 200) {
        throw new Error(`GitHub API returned ${status}`);
      }

      return {
        source: 'github',
        repository: repoPath,
        stars: data.stargazers_count,
        forks: data.forks_count,
        openIssues: data.open_issues_count,
        watchers: data.subscribers_count,
        lastCommit: data.pushed_at,
        createdAt: data.created_at,
        language: data.language,
        archived: data.archived,
        license: data.license?.spdx_id,
        url: data.html_url
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
          lastDay: statsResponse.data.data.last_day
        };
      }
    } catch (e) {
      // Download stats not available
      console.error('PyPI download stats not available');
    }

    return {
      source: 'pypi',
      package: packageName,
      version: data.info.version,
      summary: data.info.summary,
      author: data.info.author,
      license: data.info.license,
      homePage: data.info.home_page,
      projectUrls: data.info.project_urls,
      requiresPython: data.info.requires_python,
      downloads: downloads,
      uploadTime: data.urls[0]?.upload_time,
      url: `https://pypi.org/project/${packageName}/`
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

    return {
      source: 'conda',
      package: packageName,
      downloads: data.public_downloads || 0,
      latestVersion: data.latest_version,
      owner: data.owner.name,
      description: data.summary,
      url: `https://anaconda.org/bioconda/${packageName}`
    };
  } catch (error) {
    throw new Error(`Failed to fetch Conda metrics: ${error.message}`);
  }
}

/**
 * Fetch npm package metrics
 */
async function fetchNpmMetrics(packageName) {
  try {
    // Get package info
    const { status, data } = await httpsGet(`https://registry.npmjs.org/${packageName}`);

    if (status !== 200) {
      throw new Error(`npm API returned ${status}`);
    }

    // Try to get download stats
    let downloads = null;
    try {
      const statsUrl = `https://api.npmjs.org/downloads/point/last-month/${packageName}`;
      const statsResponse = await httpsGet(statsUrl);
      if (statsResponse.status === 200) {
        downloads = statsResponse.data.downloads;
      }
    } catch (e) {
      console.error('npm download stats not available');
    }

    const latest = data['dist-tags'].latest;
    const versionInfo = data.versions[latest];

    return {
      source: 'npm',
      package: packageName,
      version: latest,
      description: versionInfo.description,
      license: versionInfo.license,
      repository: versionInfo.repository?.url,
      homePage: data.homepage,
      downloads: downloads,
      url: `https://www.npmjs.com/package/${packageName}`
    };
  } catch (error) {
    throw new Error(`Failed to fetch npm metrics: ${error.message}`);
  }
}

/**
 * Main function
 */
async function main() {
  const args = process.argv.slice(2);

  if (args.length < 2) {
    console.error('Usage: fetch-tool-metrics.js <source> <identifier>');
    console.error('Sources: github, pypi, conda, npm');
    console.error('Examples:');
    console.error('  fetch-tool-metrics.js github nextflow-io/nextflow');
    console.error('  fetch-tool-metrics.js pypi multiqc');
    console.error('  fetch-tool-metrics.js conda fastqc');
    console.error('  fetch-tool-metrics.js npm typescript');
    process.exit(1);
  }

  const [source, identifier] = args;

  try {
    let metrics;

    switch (source.toLowerCase()) {
      case 'github':
        metrics = await fetchGitHubMetrics(identifier);
        break;
      case 'pypi':
        metrics = await fetchPyPIMetrics(identifier);
        break;
      case 'conda':
        metrics = await fetchCondaMetrics(identifier);
        break;
      case 'npm':
        metrics = await fetchNpmMetrics(identifier);
        break;
      default:
        throw new Error(`Unknown source: ${source}. Use github, pypi, conda, or npm`);
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
  fetchNpmMetrics
};
