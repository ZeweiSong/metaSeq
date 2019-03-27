# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Purpose:           assemble and annotate clustered beads
# Parameters:        8
#                    ---------------------------
#                    mode [F|B]
#                    sample directory
#                    cluster level
#                    cluster ID
#                    tag
#                    reference 1
#                    reference 2
#                    force run
#                    ---------------------------
# Author:            fangchao@genomics.cn
# Version:           V0.2
# Last modified:     04 Jan 2019 (since 04 Jan 2019)
# ===================================================================
#
mode=$1 #[F|B]
samDir=$2
level=$3
cluster=$4
tag=$5
ref1=$6
ref2=$7
force=$8

# if [ $mode == "F" ];
# then
   refDB=$6
   refBT=$7
# else
#   refSSU=$6
#   refLSU=$7
# fi

clusterFmt=`printf "%06d" $cluster`
if [ $level -gt 0 ];
then
  ASB="Assemble_Lv$level"
else
  ASB="Assemble_$level"
fi

subDir=$tag$clusterFmt
rdir=""
scaf=""
if [ $mode == "spades" ]; then
  rdir="spades";
  scaf="$samDir/$ASB/$subDir/scaffolds.fasta"
elif [ $mode == "idba" ]; then
  rdir="idba";
  scaf="$samDir/$ASB/$subDir/$rdir/scaffold.fa"
elif [ $mode == "megahit" ]; then
  rdir="megahit";
  scaf="$samDir/$ASB/$subDir/$rdir/final.contigs.fa"
fi


echo [BC] Barcode cluster assembly pipeline start:
echo [BC] info: sample directory: $samDir/$ASB/$subDir
echo [BC] info: REF Database:     $refDB

echo
if [[ -f $scaf && -z $force ]];
then
  echo "[BC] Assemble results exists. Skiped (add \$6 to force re-run)"
else
  if [ $mode == "spades" ]; then
    spaLog="$samDir/$ASB/$subDir/spades.log"
    if [[ $force || ! -f $spaLog ]];
    then
      spades.py --meta -t 8 -o $samDir/$ASB/$subDir \
      -1 $samDir/$ASB/$subDir/sort.1.fq -2 $samDir/$ASB/$subDir/sort.2.fq
    else
      spades.py --continue -o $samDir/$ASB/$subDir
      #spades.py --meta -t 48 -o $samDir/$ASB/$subDir -1 $samDir/$ASB/$subDir/sort.1.fq -2 $samDir/$ASB/$subDir/sort.2.fq
    fi
  elif [ $mode == "idba" ]; then
    fq2fa --merge --filter $samDir/$ASB/$subDir/sort.1.fq \
    $samDir/$ASB/$subDir/sort.2.fq $samDir/$ASB/$subDir/sort.pair.fa
    cmd="idba_ud -o $samDir/$ASB/$subDir/idba -r $samDir/$ASB/$subDir/sort.pair.fa \
    --mink 11 --step 22 --maxk 121 --min_contig 999 --num_threads 16"
	echo $cmd;
	$cmd
  elif [ $mode == "megahit" ]; then
    echo "[BC] megahit selected."
    megahit --k-min 21 --k-step 22 --k-max 121  \
    -1 $samDir/$ASB/$subDir/sort.1.fq -2 $samDir/$ASB/$subDir/sort.2.fq \
    -o $samDir/$ASB/$subDir/$rdir
  fi
fi


### BLAST for scaffolds ###
echo
scaf6="$samDir/$ASB/$subDir/$rdir/scaffolds.$mode.BLAST.tax.blast6"
if [[ -f $scaf6 && -z $force ]];
then
  echo "[BC] Scaffolds BLAST results exists. Skiped (add \$6 to force re-run)"
else
  echo "[BC] Scaffolds BLAST start"

  blastn -num_threads 4 -db $refBT -query $scaf \
  -out $scaf6 -outfmt 6

  echo "[BC] Adding taxonomy info for BLAST"
  metabbq anno.pl $refBT.ids $scaf6 > $scaf6.anno
  metabbq binWrite best -u -m -i $scaf6.anno
fi

echo "[BC] Done!"
#awk '$4>50' $ASB/$subDir/scaffolds.$mode.tax.blast6| perl ../src/refScore.pl -r \
# ../REF/ee_its_database.ITSx.anno -o $ASB/$subDir/scaffolds.$mode.tax.stat
