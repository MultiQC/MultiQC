# Release checklist

This checklist is for my own reference, as I forget the steps every time.

1. Check that everything is up to date and ready to go
2. Figure out what this release should be called
3. Update version numbers in code: `setup.py`, `CHANGELOG.md`
4. Link the changelog subheading to the as yet non-existant release URL. Add date.
5. Install the package again in `install` mode:

   ```bash
   pip install .
   ```

   - This removes the commit hash from the version number when MultiQC runs

6. Run using test data
   - Check for any command line or javascript errors
   - Check version numbers are printed correctly
7. Create new demo reports for the website

   - Comment out any config in `~/.multiqc_config.yaml`

     ```bash
     mv ~/.multiqc_config.yaml ~/.multiqc_config.yaml.bkup
     ```

     ```bash
     mv ~/.multiqc_config.yaml.bkup ~/.multiqc_config.yaml
     ```

   - Spot any previously unnoticed bugs and fix
   - Upload to the website and push change to Git

8. Commit and push version updates
9. Make a [release](https://github.com/ewels/MultiQC/releases) on GitHub - paste changelog section.
10. Check that [PyPI listing page](https://pypi.python.org/pypi/multiqc/) looks sane
11. Update version numbers to new dev version in `setup.py` + a new section in the changelog for the development version
12. Commit and push version bump
13. Make a new release on `bioconda` (assuming new modules were added):

    ```bash
    # Update to latest bioconda
    cd ../bioconda-recipes
    git checkout master
    git pull upstream master
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
    curl -OL https://github.com/ewels/MultiQC/archive/v1.5.tar.gz
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

14. Tell UPPMAX about the new version and ask for the module system to be updated.
15. Tweet that new version is released
16. Continue making more awesome :metal:
