---
name: FastQ Screen
url: http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
description: Screens a library of sequences in FastQ format against a set of sequence databases to see if the composition of the library matches with what you expect
---

By default, the module creates a plot that emulates the FastQ Screen output
with blue and red stacked bars showing unique and multimapping read counts.
This plot only works for a handful of samples however, so if
`# samples * # organisms >= 160`, a simpler stacked barplot is shown. This
is also shown when generating flat-image plots.

To always show this style of plot, add the following line to a MultiQC config file:

```yaml
fastqscreen_simpleplot: true
```
