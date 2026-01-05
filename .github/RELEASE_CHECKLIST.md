# Release checklist

This checklist is for my own reference, as I forget the steps every time.

1. Check that everything is up-to-date and ready to go
2. Update version numbers in `pyproject.toml`
3. Generate a new changelog section stub (make sure to export a GitHub token to avoid rate limits):

   ```bash
   GITHUB_TOKEN=$(gh auth token) python scripts/print_changelog.py
   ```

   Then paste it into `CHANGELOG.md` and edit accordingly: group changes if needed, add highlights.

4. Update module documentation markdown files and update the config JSON schema:

   ```bash
   cd test-data && git pull && cd ../
   python scripts/make_module_docs.py
   python scripts/generate_config_schema.py
   ```

5. Install the package again in the `install` mode. Make sure to remove the build directory and egg info that might cache outdated metadata such as entry point names:

   ```bash
   rm -rf build/ *.egg-info
   pip install .
   ```

   This removes the commit hash from the version number when MultiQC runs.

6. Run using test data
   - Check for any command line or javascript errors
   - Check version numbers are printed correctly

7. Create new demo reports for the website
   - Comment out any config in `~/.multiqc_config.yaml`

     ```bash
     mv ~/.multiqc_config.yml ~/.multiqc_config.yml.bkup
     ```

   - Generate reports in the seqeralabs/web repo.

     ```bash
     cd ../../seqera/web/services/website/public/examples
     bash update_examples.sh
     ```

   - Put back homedir config

     ```bash
     mv ~/.multiqc_config.yml.bkup ~/.multiqc_config.yml
     ```

   - Spot any previously unnoticed bugs and fix
   - Upload to the website and push change to Git

8. Commit and push version updates
9. Generate new rich-codex screenshots
   - On github.com, make a new branch using the branch dropdown called `rich-codex` (or whatever).
   - Go to _Actions_ and then [_Docs screenshots_](https://github.com/MultiQC/MultiQC/actions/workflows/screenshots.yml)
   - Click _Run Workflow_ top right, and **select your new branch**
   - Click Run. Wait for the action to complete.
   - Make a PR from this branch to `main` and check that the new screenshot looks ok. Merge if so.
10. Make sure that all tests are passing
11. Make a [release](https://github.com/MultiQC/MultiQC/releases) on GitHub - paste changelog section.
12. Check that [PyPI listing page](https://pypi.python.org/pypi/multiqc/) looks sane
13. Update version numbers to new dev version in `pyproject.toml`
14. Commit and push version bump
15. Run Seqera `#product-releases` [Slack workflow](https://slack.com/shortcuts/Ft06GYSX4UUB/c3733786a0ad2fc1794d1959aed5df19)
16. Look for the automated release PR on `bioconda` and approve / merge
    - IMPORTANT: If any new dependencies added, need a manual PR to add them.
    - Do this quickly, as other people merge the automated PRs really quickly
17. Tweet that new version is released
18. Continue making more awesome :metal:

## Appendix

### Rebuilding Bioconda recipe from scratch

Instructions for complete rebuild of Bioconda:

<details>

```bash
# Update to latest bioconda
cd ../bioconda-recipes
git checkout main
git pull upstream main
git push
git branch -D multiqc
# Build new conda recipe from PyPI to automatically collect new dependencies
git checkout -b multiqc
# Do the conda skeleton to copy the dependencies
cd recipes && mkdir mqctemp && cd mqctemp && code .
conda skeleton pypi multiqc
# Update with new release header - see https://goo.gl/ZfRnmj
cd ../multiqc && code .
# Get the sha256sum of the release
curl -OL https://github.com/MultiQC/MultiQC/archive/v1.5.tar.gz
shasum --algorithm 256 v1.5.tar.gz
# Switch out download for GitHub release and remove all other cruft
# commit changes
cd ../../
git commit -am "MultiQC version 1.23 release"
# Test locally
docker pull bioconda/bioconda-utils-build-env
circleci build
# Push updates
git push -u origin multiqc
# Submit a Pull Request and merge
```

</details>
