---
Name: BBMap
URL: http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/
Description: >
	BBMap is a suite of pre-processing, assembly, alignment, and statistics
	tools for DNA/RNA sequencing reads.
---

The BBMap module produces summary statistics from the 
[BBMap](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/) suite of tools. 
The module can summarize data from the following BBMap output files
(desriptions from `bbmap.sh` help output).
Coverage:
    covstats=<file>     Per-scaffold coverage info.
    rpkm=<file>         Per-scaffold RPKM/FPKM counts.
    covhist=<file>      Histogram of # occurrences of each depth level.
    basecov=<file>      Coverage per base location.
    bincov=<file>       Print binned coverage per location (one line per X bases).
Histogram and statistics:
    scafstats=<file>    Statistics on how many reads mapped to which scaffold.
    refstats=<file>     Statistics on how many reads mapped to which reference
                        file; only for BBSplit.
    bhist=<file>        Base composition histogram by position.
    qhist=<file>        Quality histogram by position.
    qchist=<file>       Count of bases with each quality value.
    aqhist=<file>       Histogram of average read quality.
    bqhist=<file>       Quality histogram designed for box plots.
    lhist=<file>        Read length histogram.
    gchist=<file>       Read GC content histogram.
    indelhist=<file>    Indel length histogram.
    mhist=<file>        Histogram of match, sub, del, and ins rates by
                        read location.
    statsfile=<file>    Mapping statistics are printed here.

Additional information on the BBMap tools is available on 
[SeqAnswers](http://seqanswers.com/forums/showthread.php?t=41057)
