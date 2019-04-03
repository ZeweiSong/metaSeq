# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Purpose:           align reads to assembly and assembly to ref
# Parameters:        1
#                    ---------------------------
#                    sample name (the directory name)
#                    cluster number
#                    continue from break point if defined
#                    ---------------------------
# Author:            fangchao@genomics.cn
# Version:           V0.11
# Last modified:     19 Feb 2019 (since 19 Jan 2019)
# ===================================================================
#
if [ ! $1 ];
then
  echo -e "Usage:\n  reAlign.sh [threads] <R1> <R2> <scaffolds> <reference> <output dir> [force]" && exit
fi
thread=$1
r1=$2
r2=$3
scaf=$4
refDB=$5
outDir=$6
force=$7

echo [BC] Post aligment pipeline start:
echo [BC] info: Scaffolds:     $scaf
echo [BC] info: REF Database:  $refDB

mkdir -p $outDir
echo [BC] info: output dir :      $outDir

echo
sam="$outDir/read2scaf.sam"
if [[ -f $sam && -z $force ]];
then
  echo "[BC] BC alignemnt results exists. Skiped (add \$3 to force re-run)"
else
  echo "[BC] BC alignemnt start. Indexing scaffolds.fasta"
  #awk -F "_" '($0~/^>/){if($6>=50&&$4>=100){s=1}else{s=0}};(s==1){print}' \
  awk -F "_" '($0~/^>/){if($4>=100){s=1;S=S+1}else{s=0}};(s==1 && S<=30){print}' \
   $scaf > $outDir/scaffolds.fasta
  bwa index $scaf -p $outDir/scaffolds.bwa
  echo "[BC] bwa mem aligning"
  bwa mem -a -v 3 -t $thread $outDir/scaffolds.bwa $r1 $r2 -o $sam
fi


### BLAST for scaffolds to ref ###
echo
scaf6="$outDir/scaf2ref.blast6"
if [[ -f $scaf6 && -z $force ]];
then
  echo "[BC] BLAST results exists. Skiped (add \$6 to force re-run)"
else
  echo "[BC] BLAST start"
  echo "[CMD] blastn -num_threads $thread -db $refDB -query $scaf -out $sam6 -outfmt 6"
  blastn -num_threads $thread -db $refDB -query $scaf -out $scaf6 -outfmt 6
fi


echo "[BC] Done!"
