# practice refere to :https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
threads=$1
fq=$2
pfx=$4 #contained dir and prefix
force=$5

if [[ -z $1 ]];
then
  echo "metabbq template.vsearch.preCluster.sh threads fq1 fq2 pfx force" && exit
fi
echo VSEARCH start:
echo info: reference: "$REF"
echo "info: Sample:    $fq"
outDir=`dirname $pfx`
baseN=`basename $pfx`
echo info: output dir: "$outDir"
echo info: output pfx: "$baseN"
mkdir -p "$outDir"

echo
if [[ -f $pfx.merge.fasta && -z $force ]];
then
  echo "Merge fasta exists. Skiped (add \$5 to force re-run)"
else
  echo Merge seq
  vsearch --threads $threads --fastq_mergepairs $fq1 --reverse $fq2 \
  --fastaout $pfx.merge.fasta \
  --fastqout_notmerged_fwd $pfx.unmergee.F.fasta \
  --fastqout_notmerged_rev $pfx.unmergee.R.fasta
fi

echo
if [[ -f $pfx.merge.derep.fasta && -z $force ]];
then
  echo "Dereplicate file exists. Skiped (add \$5 to force re-run)"
else
  echo "Dereplicate (non-singleton)"
  vsearch --threads $threads --derep_fulllength $pfx.merge.fasta  \
  --output $pfx.merge.derep.fasta \
  --sizeout --uc $pfx.merge.derep.uc --fasta_width 0
fi

echo
if [[ -f $pfx.merge.derep.preclustered.uc && -z $force ]];
then
  echo "Precluster result exists. Skiped (add \$5 to force re-run)"
else
  echo Precluster at 97% before chimera detection
  vsearch --threads $threads --cluster_size $pfx.merge.derep.fasta \
  --id 0.97 --strand plus --sizein --sizeout --fasta_width 0 \
  --uc $pfx.merge.derep.preclustered.uc --minuniquesize 2 \
  --centroids $pfx.merge.derep.preclustered.fasta
  echo "Unique (non-singleton) sequences:" $(grep -c "^>" $pfx.merge.derep.preclustered.fasta)
fi


echo
if [[ -f $pfx.merge.derep.uc.bc && -z $force ]];
then
  echo "Beads Cluster from duplication result found. Skiped (add \$5 to force re-run)"
else
  echo "Beads Cluster from duplication"
  metabbq beadStat.pl -match 100 -m uc -i $pfx.merge.derep.uc -o $pfx.merge.derep.uc.bc -v
fi

echo
if [[ -f $pfx.merge.derep.preclustered.uc.bc && -z $force ]];
then
  echo "Beads Cluster from preCluster result found. Skiped (add \$5 to force re-run)"
else
  echo "Beads Cluster from preCluster"
  metabbq beadStat.pl -m uc -ident 100 -match 100 -i $pfx.merge.derep.preclustered.uc -o \
  $pfx.merge.derep.preclustered.uc.bc -v
fi

echo
if [[ -f $pfx.merge.derep.2.bc && -z $force ]];
then
  echo "Beads Cluster seems already merged. Skiped (add \$5 to force re-run)"
else
  echo "Beads Cluster merging"
  metabbq beadStat.pl -m merge -i $pfx.merge.derep.uc.bc,$pfx.merge.derep.preclustered.uc.bc \
  -o $pfx.merge.derep.2.bc -v
  awk '$4>1{print $0}' $pfx.merge.derep.2.bc | sort -nrk4,4 > $pfx.merge.derep.2T.bc.s4
fi

awk -F '\t' '$4>1{print $2"\t"$3"\t"100/$4/$5}' $pfx.merge.derep.2T.bc.s4 > $pfx.merge.derep.2T.bc.dist

echo
echo "Find communities"
metabbq filterMASHoutput.pl -i $pfx.merge.derep.2T.bc.dist -o $pfx.merge.derep.2T.bc.reNum.dist -M $pfx.merge.derep.2T.bc.map
convert -i $pfx.merge.derep.2T.bc.reNum.dist -w $pfx.merge.derep.2T.bc.w -o $pfx.merge.derep.2T.bc.bin
community $pfx.merge.derep.2T.bc.bin -w $pfx.merge.derep.2T.bc.w -l -1 > $pfx.merge.derep.2T.bc.graph.tree

echo
echo "Extracting Beads Cluster Reads according to the max level"
#maxLv=`hierarchy $pfx.merge.derep.2T.bc.graph.tree -n|head -1|sed 's/Number of levels: //'`
#((maxLv=$maxLv -1 ))
maxLv=1

hierarchy $pfx.merge.derep.2T.bc.graph.tree -l $maxLv > $pfx.merge.derep.2T.bc.node\_Lv$maxLv.lst
metabbq change.id.pl -n $pfx.merge.derep.2T.bc.node\_Lv$maxLv.lst \
-m $pfx.merge.derep.2T.bc.map -v -o $pfx.merge.derep.2T.bc.cluster\_Lv$maxLv

cut -d " " -f2 $pfx.merge.derep.2T.bc.node\_Lv$maxLv.lst|sort|uniq -c|sort -nr > $pfx.merge.derep.2T.bc.node\_Lv$maxLv.rank

#idxFile=${fq1/gz/gz.idx}
idxFile=$fq1.idx
echo
if [[ -f $idxFile && -z $force ]];
then
  echo "Beads barcode index found."
else
  echo "Beads barcode index haven't found. Generate now..."
  metabbq beadsWrite3.pl -x $fq1
fi

perl -e 'open IN,"<'$idxFile'";while(<IN>){chomp;@a=split;$num=($a[2]-$a[1]+1)/4;$HASH{$a[0]}=$num};while(<>){chomp;@a=split;$HB{$a[0]}{R}+=$HASH{$a[1]};$HB{$a[0]}{C}++};foreach my $c (sort {$a<=>$b} keys %HB){print "$c\t$HB{$c}{C}\t$HB{$c}{R}\n"}' <  $pfx.merge.derep.2T.bc.cluster\_Lv$maxLv >  \
$pfx.merge.derep.2T.bc.cluster\_Lv$maxLv.count

#sort -k3,3nr $pfx.merge.derep.2T.bc.cluster\_Lv$maxLv.count| \
#awk 'FNR==1{b=$2;r=$3}($2>b/10 && $3>r/10 && FNR<300){b=$2;r=$3;print $0}'|sort -k1,1n \

#awk 'BEGIN{b=$2;r=$3}($2<=50&&$3>=50){b=$2;r=$3;print $0}'|sort -k1,1n \
sort -k3,3nr $pfx.merge.derep.2T.bc.cluster\_Lv$maxLv.count| \
awk 'BEGIN{b=$2;r=$3}($3>=50){b=$2;r=$3;print $0}'|sort -k1,1n \
> $pfx.merge.derep.2T.bc.cluster\_Lv$maxLv.count.main
perl -e 'open IN,"sort -k1,1n '$pfx'.merge.derep.2T.bc.cluster\_Lv'$maxLv'.count.main|";while(<IN>){@a=split(/\t/,$_);push @ids, $a[0]}; $tag= shift @ids; while(<>){@a=split /\t/, $_; if($a[0] < $tag){next}elsif($a[0] == $tag){print $_}else{$tag= shift @ids;print $_ if $a[0] == $tag } }' < $pfx.merge.derep.2T.bc.cluster\_Lv$maxLv |sort -k1,1n > $pfx.merge.derep.2T.bc.cluster\_Lv$maxLv.main

cp $pfx.merge.derep.2T.bc.cluster\_Lv$maxLv.main $pfx.merge.derep.2T.bc.cluster_LvMax.main



#Find uniq reads in each bead
perl -e 'while(<>){chomp;@a=split;$num=($a[2]-$a[1]+1)/4;print "$a[0]\t$num\n"}' \
< $idxFile > $pfx.in.beads.tsv
metabbq beadStat.pl -m bc -i $pfx.merge.derep.uc -o $pfx.merge.derep.bead.uc \
|sort -k1,1 > $pfx.merge.derep.bead.dup

awk 'BEGIN{bc=0;d=0;c=0}(bc!=$1&&FNR>1){print bc"\t"c"\t"d;d=0;c=0}{c=c+1;d=d+$2;bc=$1}END{print bc"\t"c"\t"d}' $pfx.merge.derep.bead.dup > $pfx.merge.derep.bead.dup.uniq

metabbq binWrite uniq -t $pfx.in.beads.tsv -c $pfx.merge.derep.2T.bc.cluster_Lv$maxLv.main \
-u $pfx.merge.derep.bead.dup.uniq -o $pfx.in.beads.uniq.tsv -v

awk 'BEGIN{c=0}($1!~/0000/ && $3>=50){print c"\t"$0;c=c+1}' $pfx.in.beads.uniq.tsv >  $pfx.individual.beads.list



# Stop at here
