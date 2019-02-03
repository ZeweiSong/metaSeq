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
# Version:           V0.1
# Last modified:     19 Jan 2019 (since 19 Jan 2019)
# ===================================================================
#
thread=$1
BCDIR=$2
refDB=$3
force=$4

outDir="$BCDIR/align"

echo [BC] Post aligment pipeline start:
echo [BC] info: sample directory: $BCDIR
echo [BC] info: REF Database:     $refDB

mkdir -p $BCDIR/align
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
   $BCDIR/scaffolds.fasta > $outDir/scaffolds.fasta
  bwa index $outDir/scaffolds.fasta -p $outDir/scaffolds.bwa
  echo "[BC] bwa mem aligning"
  bwa mem -a -v 3 -t $thread $outDir/scaffolds.bwa $BCDIR/sort.{1,2}.fq -o $sam
fi


### BLAST for scaffolds to ref ###
echo
scaf6="$outDir/scaf2ref.blast6"
if [[ -f $scaf6 && -z $force ]];
then
  echo "[BC] BLAST results exists. Skiped (add \$6 to force re-run)"
else
  echo "[BC] BLAST start"
  echo "[CMD] blastn -num_threads $thread -db $refDB -query $BCDIR/scaffolds.fasta -out $sam6 -outfmt 6"
  blastn -num_threads $thread -db $refDB -query $BCDIR/scaffolds.fasta -out $scaf6 -outfmt 6
fi


echo "[BC] Done!"
