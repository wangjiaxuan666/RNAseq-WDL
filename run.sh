#!/bin/bash

#---------Wang Jiaxuan 2022.03.08------------
# if run_RNAseq run failured but don't know why
# can use the code below this line
#--------------------------------------------


# java \
# -jar WDL/cromwell-58.jar \
# run WDL/rna_seq.wdl \
# -i WDL/input.RNAseq.json \
# -m WDL/output.RNAseq.json

./run_RNAseq -i WDL/input.RNAseq.json