---
Name: VerifyBAMID
URL: https://genome.sph.umich.edu/wiki/VerifyBamID
Description: >
    VerifyBamID checks whether reads match known genotypes or are
    contaminated as a mixture of two samples.
---
A key step in any genetic analysis is to verify whether data being generated matches expectations. verifyBamID checks whether reads in a BAM file match previous genotypes for a specific sample. In addition, it detects possible sample mixture from population allele frequency only, which can be particularly useful when the genotype data is not available.

Using a mathematical model that relates observed sequence reads to an hypothetical true genotype, verifyBamID tries to decide whether sequence reads match a particular individual or are more likely to be contaminated (including a small proportion of foreign DNA), derived from a closely related individual, or derived from a completely different individual.

This module currently only imports data from the `.selfSM` output.
The chipmix and freemix columns are imported into the general statistics table.
A verifyBAMID section is then added, with a table containing the entire selfSM file.

If no chip data was parsed, these columns will not be added to the MultiQC report.

Should you wish to remove one of these columns from the general statistics table add the below lines to the table_columns_visible section of your config file

    table_columns_visible:
        verifyBAMID:
            CHIPMIX: False
            FREEMIX: False



This was designed to work with verifyBamID 1.1.3 January 2018
