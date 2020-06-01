---
Name: HOPS
URL: https://www.github.com/rhubler/HOPS
Description: >
    This tool performs screening of output from the ancient DNA optimised 
    BLAST-replacement tool MALT, to identify taxa that have expected ancient 
    DNA characteristics.

---

This module takes the JSON output of the HOPS postprocessing R script (Version 
\>= 0.34). to recreate the possible positives heatmap, with the heat intensity
representing the number of 'ancient DNA characteristics' categories (small
edit distance, damage, both edit distance and aDNA damage) that a particular 
taxon has.