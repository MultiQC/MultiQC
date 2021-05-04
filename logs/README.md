Test output files for sambamba markdup were generated from snakemake, but basic commands would be:

sambamba markdup -r -t {threads} --tmpdir=data/markd --io-buffer-size=512 {input} {output} > {log} 2>&1
