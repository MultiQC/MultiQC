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

1. Check that everything is up to date and read to go
2. Figure out what this release should be called
3. Update version numbers in code: `setup.py`, `CHANGELOG.md`, `docs/README.md`
4. Install the package again for local development, but with the new version number:
```
python setup.py develop
```
5. Run using test data
  * Check for any command line or javascript errors
  * Check version numbers are printed correctly
6. Release on PyPI: 
```
python setup.py sdist upload
```
7. Test that it pip installs:
```
conda create --name testing python=2.7
pip install multiqc
multiqc .
source deactivate
conda remove --name testing --all
```
8. Commit and push version updates
9. Make a [release](https://github.com/ewels/MultiQC/releases) on GitHub - paste changelog section.
10. Check that [PyPI listing page](https://pypi.python.org/pypi/multiqc/) looks sane
11. Update version numbers to new dev version in files mentioned above.
12. Commit and push. Continue making more awesome :metal:
