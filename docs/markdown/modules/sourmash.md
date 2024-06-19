---
name: sourmash
url: https://github.com/sourmash-bio/sourmash
description: Quickly search, compare, and analyze genomic and metagenomic data sets.
---

The sourmash module produces summary statistics from the
[sourmash](https://github.com/sourmash-bio/sourmash) tool.
The module can summarise data from the following sourmash output files
(descriptions from command line help output):

- `sourmash compare`
  - create a similarity matrix comparing many samples.
- `sourmash gather`
  - search a metagenome signature against databases.

Additional information on sourmash and its outputs is available on
the [sourmash documentation website](https://sourmash.readthedocs.io/en/latest/).

`sourmash gather` is modelled after the Kraken module, and builds a bar graph that
shows the coverage of top-5 genomes covered most by all samples. The number of top
genomes can be customized in the config file:

```yaml
sourmash:
  gather:
    top_n: 5
```
