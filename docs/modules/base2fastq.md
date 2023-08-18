---
name: bases2fastq
url: https://docs.elembio.io/docs/bases2fastq/
description: >
  bases2fastq convert raw data from Element AVITI system to fastq and output sequencing run statistics
---

Note that the default maximum file size is 50Mb. Some Runstats.json file will exceed that size if you have many samples or longer reads per run. If samples are missing, please check the size of all json files in the folder and change the maximum file size accordingly.

You can configure the threshold and parse your files by changing the
`log_filesize_limit` config option. For example, to parse files up to 100MB in
size, add the following to your MultiQC config file:

```yaml
log_filesize_limit: 100000000
```
