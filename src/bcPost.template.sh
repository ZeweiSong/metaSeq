# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Purpose:           assemble and annotate clustered beads
# Parameters:        1
#                    ---------------------------
#                    sample name (the directory name)
#                    cluster number
#                    continue from break point if defined
#                    ---------------------------
# Author:            fangchao@genomics.cn
# Version:           V0.1
# Last modified:     04 Jan 2019 (since 04 Jan 2019)
# ===================================================================
#
mode=$1 #[F|B]
samDir=$2
level=$3
cluster=$4
ident=$5
force=$6

if [ $mode == "F" ];
then
  refDB="REF/ee_its_database.fasta"
  refBT="REF/ee_its_database.blast"
else
  refSSU="REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta"
  refLSU="REF/silva132/LSU/SILVA_132_LSURef_tax_RNA.fasta"
fi

clusterFmt=`printf "%05d" $cluster`
if [ $level -gt 0 ];
then
  ASB="Assemble_Lv$level"
else
  ASB="Assemble"
fi

subDir=BC$clusterFmt

if [ -z $ident ];then
  ident=97
fi

echo [BC] Barcode cluster assembly pipeline start:
echo [BC] info: sample directory: $samDir/$ASB/$subDir
echo [BC] info: REF Database:     $refDB

echo
echo "[BC] list beads contained in this cluster"
#mkdir $ASB/$subDir
awk -v c=$cluster '$1==c{print $2}' $samDir/VSEARCH/read.merge.derep.2T.bc.cluster|sort > $samDir/$ASB/$subDir/beads.lst
#perl ../src/beadsWrite3.pl --r1 ./clean/fastp.sort.1.fq --r2 ./clean/fastp.sort.2.fq -b $ASB/$subDir/beads.lst -o $ASB/$subDir -p sort -f fq -v

echo
scaf="$samDir/$ASB/$subDir/scaffolds.fasta"
if [[ -f $scaf && -z $force ]];
then
  echo "[BC] BLAST SSU alignemnt results exists. Skiped (add \$6 to force re-run)"
else
  spaLog="$samDir/$ASB/$subDir/spades.log"
  if [[ $force || ! -f $spaLog ]];
  then
    spades.py --careful --cov-cutoff auto -t 8 -o $samDir/$ASB/$subDir -1 $samDir/$ASB/$subDir/sort.1.fq -2 $samDir/$ASB/$subDir/sort.2.fq
  else
    spades.py --continue -o $samDir/$ASB/$subDir
    #spades.py --meta -t 48 -o $samDir/$ASB/$subDir -1 $samDir/$ASB/$subDir/sort.1.fq -2 $samDir/$ASB/$subDir/sort.2.fq
  fi
fi

if [ $mode == 'F' ];
then
  #### ITSx ###
  echo
  itsFa="$samDir/$ASB/$subDir/scaffolds.$mode.ITSm.fasta"
  if [[ -f $itsFa && -z $force ]];
  then
    echo "[BC] ITSx results exists. Skiped (add \$6 to force re-run)"
  else
    echo "[BC] Predicting ITS positions in scaffolds"
    ITSx -i $samDir/$ASB/$subDir/scaffolds.fasta -o $samDir/$ASB/$subDir/scaffolds --cpu 8 --save_regions all --detailed_results T
    cat $samDir/$ASB/$subDir/scaffolds.{SSU,ITS1,5_8S,ITS2,LSU}.fasta | \
    awk '{if($0~/^>/){id=$0}else{len=length($0);if(len>=50){print id"\n"$0}}}' > $itsFa
  fi
  #### barrnap ###
  echo
  barrnapFa="$samDir/$ASB/$subDir/scaffolds.$mode.barrnap.fasta"
  if [[ -f $barrnapFa && -z $force ]];
  then
    echo "[BC] BLAST SSU alignemnt results exists. Skiped (add \$6 to force re-run)"
  else
    echo "[BC] Annotating SSU by BLAST"
    barrnap -kingdom euk --reject 0.1 < $samDir/$ASB/$subDir/scaffolds.fasta \
     > $samDir/$ASB/$subDir/scaffolds.barrnap.positions.txt --outseq $barrnapFa
    sed -i 's/^>\(.*\)::\(.*\):\(.*\)$/>\2|N|\1 (\3)/' $barrnapFa
  fi
  # #### usearch annotation ###
  # echo
  # usearchFile="$samDir/$ASB/$subDir/scaffolds.$mode.tax\_$ident.blast6"
  # if [[ -f $usearchFile && -z $force ]];
  # then
  #   echo "[BC] VSEARCH alignemnt results exists. Skiped (add \$6 to force re-run)"
  # else
  #   echo "[BC] Annotating by VSEARCH"
  #   vsearch --threads 8 --usearch_global $samDir/$ASB/$subDir/scaffolds.$mode.fasta --db \
  #   $refDB --id 0.$ident --maxaccepts 0 --maxrejects 0 \
  #   --blast6out $usearchFile \
  #   --samout $samDir/$ASB/$subDir/scaffolds.$mode.tax\_$ident.sam
  #
  #   echo "[BC] Adding taxonomy info for VSERACH"
  #   perl src/anno.pl $refDB.ids $usearchFile \
  #   > $samDir/$ASB/$subDir/scaffolds.$mode.tax\_$ident.blast6.anno
  # fi

  ### BLAST for predicted units ###
  echo
  units6="$samDir/$ASB/$subDir/scaffolds.$mode.BLAST.units.blast6"
  if [[ -f $units6 && -z $force ]];
  then
    echo "[BC] units BLAST results exists. Skiped (add \$6 to force re-run)"
  else
    echo "[BC] units BLAST start"
    cat $itsFa $barrnapFa > $samDir/$ASB/$subDir/scaffolds.$mode.units.fasta

    blastn -num_threads 8 -db $refBT -out $units6 -outfmt 6 \
    -query $samDir/$ASB/$subDir/scaffolds.$mode.units.fasta

    echo "[BC] units annotating"
    perl src/anno.pl $refBT.ids $units6 > $units6.anno
  fi
  ### BLAST for scaffolds ###
  echo
  scaf6="$samDir/$ASB/$subDir/scaffolds.$mode.BLAST.tax_$ident.blast6"
  if [[ -f $scaf6 && -z $force ]];
  then
    echo "[BC] Scaffolds BLAST results exists. Skiped (add \$6 to force re-run)"
  else
    echo "[BC] Scaffolds BLAST start"

    blastn -num_threads 8 -db $refBT -query $samDir/$ASB/$subDir/scaffolds.fasta \
    -out $scaf6 -outfmt 6

    echo "[BC] Adding taxonomy info for BLAST"
    perl src/anno.pl $refBT.ids $scaf6 > $scaf6.anno
  fi
else
  ###   BAC   ###
  ### barrnap ###
  echo
  barrnapFa="$samDir/$ASB/$subDir/scaffolds.$mode.barrnap.fasta"
  if [[ -f $barrnapFa && -z $force ]];
  then
    echo "[BC] barrnap prediction results exists. Skiped (add \$6 to force re-run)"
  else
    echo "[BC] barrnap predicting by BLAST"
    barrnap -kingdom bac --reject 0.1 < $samDir/$ASB/$subDir/scaffolds.fasta \
     > $samDir/$ASB/$subDir/scaffolds.barrnap.positions.txt --outseq $barrnapFa
    sed -i 's/^>\(.*\)::\(.*\):\(.*\)$/>\2|N|\1 (\3)/' $barrnapFa
  fi

  ### BLAST for predicted units ###
  echo
  units6="$samDir/$ASB/$subDir/scaffolds.$mode.BLAST.units.blast6"
  if [[ -f $units6 && -z $force ]];
  then
    echo "[BC] units BLAST results exists. Skiped (add \$6 to force re-run)"
  else
    echo "[BC] units BLAST start"
    blastn -num_threads 8 -db $refSSU -query $barrnapFa -out ${units6/units/SSU} -outfmt 6
    blastn -num_threads 8 -db $refLSU -query $barrnapFa -out ${units6/units/LSU} -outfmt 6
    cat ${units6/units/SSU} ${units6/units/LSU} > $units6
    echo "[BC] units BLAST annotating"
    perl src/anno.pl $refSSU.ids ${units6/units/SSU} > ${units6/units/SSU}.anno
    perl src/anno.pl $refLSU.ids ${units6/units/LSU} > ${units6/units/LSU}.anno
    cat ${units6/units/SSU}.anno ${units6/units/LSU}.anno > $units6.anno
  fi
  ### BLAST for scaffolds ###
  echo
  scaf6="$samDir/$ASB/$subDir/scaffolds.$mode.BLAST.tax_$ident.blast6"
  if [[ -f $scaf6 && -z $force ]];
  then
    echo "[BC] Scaffolds BLAST results exists. Skiped (add \$6 to force re-run)"
  else
    echo "[BC] Scaffolds BLAST start"
    blastn -num_threads 8 -db $refSSU -out ${scaf6/tax/SSU} -outfmt 6 \
    -query $samDir/$ASB/$subDir/scaffolds.fasta
    blastn -num_threads 8 -db $refLSU -out ${scaf6/tax/LSU} -outfmt 6 \
    -query $samDir/$ASB/$subDir/scaffolds.fasta
    cat ${scaf6/tax/SSU} ${scaf6/tax/LSU} > $scaf6
  fi
  if [[ -f $scaf6.anno && -z $force ]];
  then
    echo "[BC] Scaffolds BLAST results exists. Skiped (add \$6 to force re-run)"
  else
    echo "[BC] Adding taxonomy info for BLAST"
    perl src/anno.pl $refSSU.ids ${scaf6/tax/SSU} > ${scaf6/tax/SSU}.anno
    perl src/anno.pl $refLSU.ids ${scaf6/tax/LSU} > ${scaf6/tax/LSU}.anno
    cat ${scaf6/tax/SSU}.anno ${scaf6/tax/LSU}.anno > $scaf6.anno
  fi
fi

echo "[BC] Done!"
#awk '$4>50' $ASB/$subDir/scaffolds.$mode.tax.blast6| perl ../src/refScore.pl -r \
# ../REF/ee_its_database.ITSx.anno -o $ASB/$subDir/scaffolds.$mode.tax.stat
