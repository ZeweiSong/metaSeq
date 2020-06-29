#!/bin/bash

# Copyright 2013-2019, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Modified by Chao Fang <fangchao@genomics.cn> to fit stLFR needs
# Build a deeply customed Kraken2 database according to a merged taxononmy tree
# from SSU(138),LSU(132) Silva database and UNITE(04.04.2020) database

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

KRAKEN2_DB_NAME=$1
KRAKEN2_THREAD_CT=$2

mkdir -p "$KRAKEN2_DB_NAME"
pushd "$KRAKEN2_DB_NAME"
mkdir -p data taxonomy library
pushd data

cat mergeTaxon.txt added.txt|sort > paths.txt
build_silva_taxonomy.pl paths.txt
popd
mv data/names.dmp data/nodes.dmp data/paths.txt taxonomy/
cp data/added.txt.acc_taxid seqid2taxid.map
popd

kraken2-build --db $KRAKEN2_DB_NAME --build --threads $KRAKEN2_THREAD_CT
