#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Main Snakefile executing metaSeq pipeline
# Author:            Chao | fangchao@genomics.cn
# Version:           V0.1
# Last modified:     02 Feb 2019 (since 30 Jan 2019)
# ===================================================================

# Init
## Read config
configfile: "config.yaml"

#  Main
rule all:
    input:
        expand("beadPool/{sample}.B.filter.dist", sample=config["samples"])

# Module #01: Beadbarcods detection
include: "rules/BBprep.smk"

# Module #02: cluster beads assembly solution
include: "rules/beadsCluster.smk"

# Module #03: Assemble Draft
#include: src + "/rules/assemble.smk"

# Module #04: binning method
include: "rules/athena.smk"
