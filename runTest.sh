#!/bin/bash
python setup.py install
multiqc /Users/alexanderpeltzer/IDEA/MultiQC_TestData/data/modules/ivar/ --lint -f -m ivar 
