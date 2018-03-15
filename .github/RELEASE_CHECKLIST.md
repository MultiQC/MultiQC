# Release checklist
This checklist is for my own reference, as I forget the steps every time.

1. Check that everything is up to date and ready to go
2. Figure out what this release should be called
3. Update version numbers in code: `setup.py`, `CHANGELOG.md`
4. Link the changelog subheading to the as yet non-existant release URL. Add date.
5. Install the package again in `install` mode:
    ```bash
    python setup.py install
    ```
6. Run using test data
    * Check for any command line or javascript errors
    * Check version numbers are printed correctly
7. Create new demo reports for the website and upload.
    * Spot any previously unnoticed bugs and fix
8. Release on PyPI:
    ```bash
    rm -rf dist/
    python setup.py sdist bdist_wheel
    twine upload dist/*.tar.gz
    ```
9. Test that it pip installs:
    ```bash
    conda create --name testing --yes python pip && source activate testing
    pip install multiqc
    multiqc .
    source deactivate && conda remove --name testing --all --yes && conda clean --all --yes
    ```
10. Commit and push version updates
11. Make a [release](https://github.com/ewels/MultiQC/releases) on GitHub - paste changelog section.
12. Check that [PyPI listing page](https://pypi.python.org/pypi/multiqc/) looks sane
13. Make a new release on `bioconda`:
    ```bash
    # Update to latest bioconda
    cd ../bioconda-recipes
    git checkout master
    git pull upstream master
    git push
    git branch -d multiqc
    # Rebuild test environment
    ./simulate-travis.py --bootstrap /tmp/miniconda --overwrite
    # Build new conda recipe from PyPI to automatically collect new dependencies
    git checkout -b multiqc
    cd recipes
    rm -r multiqc
    conda skeleton pypi multiqc
    # Switch out download for GitHub release and remove all other cruft
    # Test locally
    cd ../
    ./simulate-travis.py --packages multiqc --force
    # Push updates
    git commit -am "MultiQC update"
    git push origin multiqc
    # Submit a Pull Request and merge
    ```
14. Tell UPPMAX about the new version and ask for the module system to be updated.
15. Describe new release on [SeqAnswers thread](http://seqanswers.com/forums/showthread.php?p=195831#post195831)
16. Tweet that new version is released
17. Update version numbers to new dev version in `setup.py`
18. Add a new section in the changelog for the development version
19. Commit and push. Continue making more awesome :metal:
