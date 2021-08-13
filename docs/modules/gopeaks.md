---
Name: gopeaks
URL: https://github.com/maxsonBraunLab/gopeaks
Description: >
  Peak caller for CUT&TAG sequencing data.
---

[Gopeaks](https://github.com/maxsonBraunLab/gopeaks) is designed to call peaks from aligned CUT&TAG sequencing reads. Gopeaks uses a binomial distribution to model the read counts in sliding windows across the genome and calculate peak regions that are enriched over the background.

The gopeaks module parses runtime data output from [gopeaks](https://github.com/maxsonBraunLab/gopeaks) by searching for files ending in "\_gopeaks.json".

By default the module will add data from the gopeaks.json file to the table data, and also include a barchart of the number peaks found for each sample. The following configuration options can be switched to `False` to turn off either or both display options.

```yaml
gopeaks_config:
  show_table_data: True
  show_peak_counts_plot: True
```
