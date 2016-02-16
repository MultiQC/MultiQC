# Contributing to MultiQC
Suggestions to the developer are welcome! Contributions to the code are
even more welcome ;)

## Contribution workflow:
1. Create an issue describing the bug / suggestion / improvement / desired module
   [here](https://github.com/ewels/MultiQC/issues).
    1. Please make it clear if you're planning to work on this yourself -
       this avoids someone else duplicating your work by accident.
2. Fork this repository to your GitHub account
3. Make the necessary changes / additions within your forked repository
4. Submit a Pull Request and wait for the code to be reviewed and merged.

## Release checklist
This checklist is for my own reference, as I forget the steps every time.

1. Check that everything is up to date and ready to go
2. Figure out what this release should be called
3. Update version numbers in code: `setup.py`, `CHANGELOG.md`, `docs/README.md`
4. Link the changelog subheading to the as yet non-existant release URL. Add date.
5. Install the package again for local development, but with the new version number:
```
python setup.py develop
```
6. Run using test data
  * Check for any command line or javascript errors
  * Check version numbers are printed correctly
7. Release on PyPI: 
```
python setup.py sdist upload
```
8. Test that it pip installs:
```
conda create --name testing python=2.7
pip install multiqc
multiqc .
source deactivate
conda remove --name testing --all
```
9. Commit and push version updates
10. Make a [release](https://github.com/ewels/MultiQC/releases) on GitHub - paste changelog section.
11. Check that [PyPI listing page](https://pypi.python.org/pypi/multiqc/) looks sane
12. Update version numbers to new dev version in `setup.py`, `CHANGELOG.md`, `docs/README.md`
13. Tell UPPMAX about the new version and ask for the module system to be updated.
14. Add a new section in the changelog for the development version
15. Commit and push. Continue making more awesome :metal:
