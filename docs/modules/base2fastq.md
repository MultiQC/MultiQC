---
name: bases2fastq
url: https://docs.elembio.io/docs/bases2fastq/
description: >
  bases2fastq convert raw data from Element AVITI system to fastq and output sequencing run statistics
---

Note that the default maximum file size is 10Mb. Some Runstats.json file will exceed that size if you have many samples (>=3) or long reads per run. If some samples/runs do not show up in the final report, please check the size of all json files in the folder and change the maximum file size accordingly.