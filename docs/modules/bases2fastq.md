---
name: Bases2Fastq
url: https://docs.elembio.io/docs/bases2fastq/
description: >
  Bases2Fastq converts raw data from Element AVITI system to FastQ and outputs sequencing run statistics.
---

Note that the default maximum file size is 50Mb. Some Runstats.json file will exceed that size if you have many samples or longer reads per run. If samples are missing, please check the size of all json files in the folder and change the maximum file size accordingly.

You can configure the threshold and parse your files by changing the
`log_filesize_limit` config option.
See the [documentation](https://multiqc.info/docs/usage/troubleshooting/#big-log-files) for more information.
