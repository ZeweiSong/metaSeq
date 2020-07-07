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
if [ -e hash.k2d ]; then
  rename k2d bak.k2d *.k2d
fi
mkdir -p data taxonomy library
pushd data
if [ -e mergeTaxon.txt ]; then
	cat mergeTaxon.txt added.txt|sort > paths.txt
else
	sort added.txt > paths.txt
fi
build_silva_taxonomy.pl paths.txt
popd
cp data/names.dmp data/nodes.dmp data/paths.txt taxonomy/
if [ -e data/seqid2taxid.map ];then
	cat data/{added.acc_taxid,seqid2taxid.map} > seqid2taxid.map
else
	cp data/added.acc_taxid seqid2taxid.map
fi
popd

kraken2-build --db $KRAKEN2_DB_NAME --build --threads $KRAKEN2_THREAD_CT

bracken-build -d $KRAKEN2_DB_NAME -t $KRAKEN2_THREAD_CT
