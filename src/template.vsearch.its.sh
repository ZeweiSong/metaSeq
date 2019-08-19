# practice refere to :https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
threads=$1
REF=$2
fq1=$3
fq2=$4
pfx=$5 #contained dir and prefix
pct=$6 #Percentage of cluster and alignment

echo VSEARCH start:
echo info: reference: "$REF"
echo "info: Sample:
    $fq1
    $fq2"
outDir=`dirname $pfx`
baseN=`basename $pfx`
echo info: output dir: "$outDir"
echo info: output pfx: "$baseN"
mkdir -p "$outDir"

echo
echo Merge seq
vsearch --threads $threads --fastq_mergepairs $fq1 --reverse $fq2 \
--fastaout $pfx.merge.fasta \
--fastqout_notmerged_fwd $pfx.unmergee.F.fasta \
--fastqout_notmerged_rev $pfx.unmergee.R.fasta

echo
echo "Dereplicate (non-singleton)"
vsearch --threads $threads --derep_fulllength $pfx.merge.fasta  \
--output $pfx.merge.derep.fasta \
--sizeout --uc $pfx.merge.derep.uc --fasta_width 0

echo
echo Precluster at 98% before chimera detection

vsearch --threads $threads --cluster_size $pfx.merge.derep.fasta \
--id 0.98 --strand plus --sizein --sizeout --fasta_width 0 \
--uc $pfx.merge.derep.preclustered.uc --minuniquesize 2 \
--centroids $pfx.merge.derep.preclustered.fasta

echo "Unique (non-singleton) sequences:" $(grep -c "^>" $pfx.merge.derep.preclustered.fasta)

echo
echo De novo chimera detection
vsearch --threads $threads --uchime_denovo $pfx.merge.derep.preclustered.fasta \
--sizein --sizeout --fasta_width 0 \
--nonchimeras $pfx.merge.derep.denovo.nonchimeras.fasta

echo Unique sequences after de novo chimera detection: $(grep -c "^>" $pfx.merge.derep.denovo.nonchimeras.fasta)

echo
echo Reference chimera detection
vsearch --threads $threads --uchime_ref $pfx.merge.derep.denovo.nonchimeras.fasta \
--db $REF   --sizein --sizeout --fasta_width 0 \
--nonchimeras $pfx.merge.derep.ref.nonchimeras.fasta

echo Unique sequences after reference-based chimera detection: $(grep -c "^>" $pfx.merge.derep.ref.nonchimeras.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

perl ./src/map.pl $pfx.merge.derep.fasta $pfx.merge.derep.preclustered.uc \
$pfx.merge.derep.ref.nonchimeras.fasta > $pfx.merge.derep.nonchimeras.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" $pfx.merge.derep.nonchimeras.fasta)


echo
echo Cluster nonchimeras at $pct% and relabel with OTU_n, generate OTU table
vsearch --threads $threads --cluster_size  $pfx.merge.derep.nonchimeras.fasta \
--id 0.$pct --strand plus --sizein --sizeout --fasta_width 0 --relabel OTU_ \
--centroids $pfx.merge.derep.nonchimeras.otus.$pct.fasta \
--biomout $pfx.merge.derep.nonchimeras.otutab.$pct.biom

echo
echo nonchimeras annotation at $pct%
vsearch --threads $threads --usearch_global $pfx.merge.derep.nonchimeras.otus.$pct.fasta --db $REF --id 0.$pct \
--maxaccepts 0 --maxrejects 0 --blast6out $pfx.merge.derep.nonchimeras.otus.tax.$pct.blast6 \
--samout $pfx.merge.derep.nonchimeras.otus.tax.$pct.sam

# echo
# echo stat occurrence
# awk 'FNR>1{a=0;for(i=2;i<=NF;i++){if($i>0){a=a+1}};print $1"\t"a/(NF-1)}' $pfx.merge.derep.nonchimeras.otutab.txt|\
# sort -nrk2,2 > $pfx.merge.derep.nonchimeras.otutab.rank

perl -e 'while(<>){if($_=~/"id":"(OTU_\d+)"/){$r++;$HASH{$r}{OTU}=$1;$HASH{$r}{val}=0;$HASH{$r}{abun}=0};
if($_=~/\[(\d+),(\d+),(\d+)\]/){$HASH{$1}{val}++;$HASH{$1}{abun}+=$3}};
foreach my $i (keys %HASH){print "$HASH{$i}{OTU}\t$HASH{$i}{val}\t".$HASH{$i}{val}/$r."\t".$HASH{$1}{abun}/$HASH{$i}{val}."\n"}'\
< $pfx.merge.derep.nonchimeras.otutab.$pct.biom |\
sort -nrk2,2 > $pfx.merge.derep.nonchimeras.otutab.$pct.rank

echo
# echo stat annotation OTU with occurrence > 6%
# for i in `awk '($3>0.06){print $1}' $pfx.merge.derep.nonchimeras.otutab.rank`;do
#   grep "$i;" $pfx.merge.derep.nonchimeras.otus.tax.blast6;
# done| perl src/refScore.pl -r REF/ee_its_database.ITSx.anno \
# -o $pfx.merge.derep.nonchimeras.otus.pct6.tax.stat

echo All done!
