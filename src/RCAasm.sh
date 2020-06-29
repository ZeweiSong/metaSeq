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
#                    force run
#                    ---------------------------
# Author:            fangchao@genomics.cn
# Version:           V0.0
# Last modified:     03 Jun 2020 (since 05 Jun 2019)
# ===================================================================
#
mode=$1 #[F|B]
samDir=$2
level=$3
cluster=$4
tag=$5
mem=$6
cpu=$7
config=$8
force=$9


clusterFmt=`printf "%08.0f" $cluster`
if [[ $level -gt 0 ]];
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
  scaf="$samDir/$ASB/$subDir/$rdir/scaffolds.fasta"
elif [ $mode == "idba" ]; then
  rdir="idba";
  scaf="$samDir/$ASB/$subDir/$rdir/scaffold.fa"
elif [ $mode == "megahit" ]; then
  rdir="megahit";
  scaf="$samDir/$ASB/$subDir/$rdir/final.contigs.fa"
fi
clip="$samDir/$ASB/$subDir/$rdir/RCAclip.fa"

echo [BC] Barcode cluster assembly pipeline start:
echo [BC] info: sample directory: $samDir/$ASB/$subDir

echo
echo "[BC] list beads contained in this cluster"

fq1="$samDir/$ASB/$subDir/sort.1.fq"
fq2="$samDir/$ASB/$subDir/sort.2.fq"
#fo1="$samDir/$ASB/$subDir/noAd.1.fq"
#fo2="$samDir/$ASB/$subDir/noAd.2.fq"
fo1=$fq1
fo2=$fq2

#perl -e 'open R,"<'$fq2'";open O,">'$fo2'";while(<>){
#$i2=<R>;$s2=<R>;<R>;$q2=<R>;$id=$_;$s=<>;
#if($s=~/GCCAGCAACCGCGGTAA|TTACCGCGGTTGCTGGC|GGTCCGTGTTTCAAGACG|CGTCTTGAAACACGGACC/ || $s2=~/GCCAGCAACCGCGGTAA|TTACCGCGGTTGCTGGC|GGTCCGTGTTTCAAGACG|CGTCTTGAAACACGGACC/){<>;<>}else{<>;$q=<>;print "$id"."$s+\n$q";print O "$i2"."$s2+\n$q2"}}' < $fq1 > $fo1

#mkdir $ASB/$subDir
if [ $level == "BC" ];
then
  spadesMode="--meta"
elif [ $level == "BI" ];
then
  spadesMode="--sc"
else
  fq="$samDir/$ASB/$subDir/sort.1.fq"
  metabbq beadsWrite3.pl -x --r1 $fq
  awk '$1!~/0000/{print ($3-$2+1)/4}' $fq.idx|sort|uniq -c|sort -k2,2nr|awk '{b=b+$1;print $0"\t"b}' > $fq.stat
fi
#perl ../src/beadsWrite3.pl --r1 ./clean/fastp.sort.1.fq --r2 ./clean/fastp.sort.2.fq -b $ASB/$subDir/beads.lst -o $ASB/$subDir -p sort -f fq -v

echo
if [[ -f $scaf && -z $force ]];
then
  echo "[BC] Assemble results exists. Skiped (add \$9 to force re-run)"
else
  if [ $mode == "spades" ]; then
    spaLog="$samDir/$ASB/$subDir/spades/spades.log"
	mkdir -p $samDir/$ASB/$subDir/spades
    if [[ $force || ! -f $spaLog ]];
    then
      spades.py $spadesMode -t 8 -o $samDir/$ASB/$subDir/spades \
      -1 $samDir/$ASB/$subDir/sort.1.fq -2 $samDir/$ASB/$subDir/sort.2.fq
    else
      spades.py --continue -o $samDir/$ASB/$subDir/spades
      #spades.py --meta -t 48 -o $samDir/$ASB/$subDir -1 $samDir/$ASB/$subDir/sort.1.fq -2 $samDir/$ASB/$subDir/sort.2.fq
    fi
  elif [ $mode == "idba" ]; then
    fq2fa --merge --filter $samDir/$ASB/$subDir/sort.1.fq \
    $samDir/$ASB/$subDir/sort.2.fq $samDir/$ASB/$subDir/sort.pair.fa
    cmd="idba_ud -o $samDir/$ASB/$subDir/idba -r $samDir/$ASB/$subDir/sort.pair.fa \
    --mink 11 --step 22 --maxk 121 --min_contig 500 --num_threads $cpu"
	echo $cmd;
	$cmd
  elif [ $mode == "megahit" ]; then
    echo "[BC] megahit selected."
    if [[ -d $samDir/$ASB/$subDir/$rdir ]];
    then
      cmd="rm -rf $samDir/$ASB/$subDir/$rdir.old && rename $rdir $rdir.old $samDir/$ASB/$subDir/$rdir"
      echo "[BC] $cmd";
      rm -rf "$samDir/$ASB/$subDir/$rdir.old"
      rename $rdir $rdir.old $samDir/$ASB/$subDir/$rdir
    fi
#    megahit --k-min 21 --k-step 22 --k-max 121  -m $mem -t $cpu \
    cmd="megahit --k-min 21 --k-step 20 --prune-level 0 --min-count 1 -m $mem -t $cpu \
    -1 $samDir/$ASB/$subDir/sort.1.fq -2 $samDir/$ASB/$subDir/sort.2.fq \
    -o $samDir/$ASB/$subDir/$rdir &> $samDir/$ASB/$subDir/megahit.log"
    echo $cmd;
    $cmd
  fi
fi

# RUN RCA CLIP
echo
if [[ -f $clip && -z $force ]];
then
  echo "[BC] RCACLIP results exists. Skiped (add \$9 to force re-run)"
else
  #read primers from config
  while read line;do eval "$line"; done < $config
  echo
  echo "Run RCACLIP"
  cmd="metabbq RCACLIPS -i $scaf -o $clip -fwd $FWD -rev $REV -a $Ad -v 2> ${clip/fa/log}"
  echo "exec: $cmd"
  $cmd
  echo "DONE"
fi
