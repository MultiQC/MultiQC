# Release checklist

This checklist is for my own reference, as I forget the steps every time.

1. Check that everything is up-to-date and ready to go
2. Update version numbers in `pyproject.toml`
3. Run `python scripts/print_changelog.py` to generate a new changelog section stub, paste it into `CHANGELOG.md` and edit accordingly: group changes if needed, add highlights.
4. Install the package again in `install` mode:

   ```bash
   pip install .
   ```

   - This removes the commit hash from the version number when MultiQC runs
   - If still getting the commit hash in the version, check that the `venv` isn't in a subdirectory of the cloned MultiQC git repo

5. Run using test data
   - Check for any command line or javascript errors
   - Check version numbers are printed correctly
6. Create new demo reports for the website

   - Comment out any config in `~/.multiqc_config.yaml`

     ```bash
     mv ~/.multiqc_config.yml ~/.multiqc_config.yml.bkup
     ```

   - Install `NationalGenomicsInfrastructure/MultiQC_NGI` - **NEEDS Python 3.11**:

     ```bash
     pip install git+https://github.com/NationalGenomicsInfrastructure/MultiQC_NGI@0.7.1
     ```

   - Generate reports in the multiqc/website repo.

     ```bash
     bash update_examples.sh
     ```

   - Put back homedir config

     ```bash
     mv ~/.multiqc_config.yml.bkup ~/.multiqc_config.yml
     ```

   - Spot any previously unnoticed bugs and fix
   - Upload to the website and push change to Git

7. Commit and push version updates
8. Generate new rich-codex screenshots
   - On github.com, make a new branch using the branch dropdown called `rich-codex` (or whatever).
   - Go to _Actions_ and then [_Docs screenshots_](https://github.com/MultiQC/MultiQC/actions/workflows/screenshots.yml)
   - Click _Run Workflow_ top right, and **select your new branch**
   - Click Run. Wait for the action to complete.
   - Make a PR from this branch to `main` and check that the new screenshot looks ok. Merge if so.
9. Make sure that all tests are passing
10. Make a [release](https://github.com/MultiQC/MultiQC/releases) on GitHub - paste changelog section.
11. Check that [PyPI listing page](https://pypi.python.org/pypi/multiqc/) looks sane
12. Update version numbers to new dev version in `setup.py` + a new section in the changelog for the development version
13. Commit and push version bump
14. Run Seqera `#product-releases` [Slack workflow](https://slack.com/shortcuts/Ft06GYSX4UUB/c3733786a0ad2fc1794d1959aed5df19)
15. Look for the automated release PR on `bioconda` and approve / merge
    - IMPORTANT: If any new dependencies added, need a manual PR to add them.
    - Do this quickly, as other people merge the automated PRs really quickly
16. Tweet that new version is released
17. Continue making more awesome :metal:

## Appendix

### Rebuilding BioConda recipe from scratch

Instructions for complete rebuild of BioConda:

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
