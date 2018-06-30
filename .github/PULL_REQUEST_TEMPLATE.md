Many thanks to contributing to MultiQC!

Please fill in the appropriate checklist below (delete whatever is not relevant). These are the most common things I request on pull requests (PRs).

## If this PR is _not_ a new module
 - [ ] This comment contains a description of changes (with reason)
 - [ ] `CHANGELOG.md` has been updated
 - [ ] (optional but recommended): https://github.com/ewels/MultiQC_TestData contains test data for this change

## If this PR is for a new module
 - [ ] There is example tool output for tools in the https://github.com/ewels/MultiQC_TestData repository
 - [ ] Code is tested and works locally (including with `--lint` flag)
 - [ ] `CHANGELOG.md` is updated
 - [ ] `README.md` is updated
 - [ ] `docs/README.md` is updated with link to below
 - [ ] `docs/modulename.md` is created
 - [ ] Everything that can be represented with a plot instead of a table is a plot
 - [ ] Report sections have a description and help text (with `self.add_section`)
 - [ ] There aren't any huge tables with > 6 columns (explain reasoning if so)
 - [ ] Each table column has a different colour scale to its neighbour, which relates to the data (eg. if high numbers are bad, they're red)
 - [ ] Module does not do any significant computational work
