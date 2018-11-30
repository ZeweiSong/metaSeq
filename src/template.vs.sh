# practice refere to :https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
threads=$1
REF=$2
sam=$3
community=$4
ph="_"
echo VSEARCH start:
echo info: Sample: "$sam | community = $community | REF = $REF "
echo info: output dir: "$sam/community$ph$community/"

pfx="$sam/community$ph$community/scaffolds"

echo
echo "Dereplicate (non-singleton)"
vsearch --threads $threads --derep_fulllength $pfx.fasta  \
--output $pfx.derep.fasta \
--sizeout --uc $pfx.derep.uc --fasta_width 0

echo
echo Precluster at 98% before chimera detection

vsearch --threads $threads --cluster_size $pfx.derep.fasta \
--id 0.98 --strand plus --sizein --sizeout --fasta_width 0 \
--uc $pfx.derep.preclustered.uc --minuniquesize 2 \
--centroids $pfx.derep.preclustered.fasta

echo "Unique (non-singleton) sequences:" $(grep -c "^>" $pfx.derep.preclustered.fasta)

echo
echo De novo chimera detection
vsearch --threads $threads --uchime_denovo $pfx.derep.preclustered.fasta \
--sizein --sizeout --fasta_width 0 \
--nonchimeras $pfx.derep.denovo.nonchimeras.fasta

echo Unique sequences after de novo chimera detection: $(grep -c "^>" $pfx.derep.denovo.nonchimeras.fasta)

echo
echo Reference chimera detection
vsearch --threads $threads --uchime_ref $pfx.derep.denovo.nonchimeras.fasta \
--db $REF   --sizein --sizeout --fasta_width 0 \
--nonchimeras $pfx.derep.ref.nonchimeras.fasta

echo Unique sequences after reference-based chimera detection: $(grep -c "^>" $pfx.derep.ref.nonchimeras.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

perl ../src/map.pl $pfx.derep.fasta $pfx.derep.preclustered.uc \
$pfx.derep.ref.nonchimeras.fasta > $pfx.derep.nonchimeras.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" $pfx.derep.nonchimeras.fasta)


echo
echo Cluster nonchimeras at 97% and relabel with OTU_n, generate OTU table
vsearch --threads $threads --cluster_size  $pfx.derep.nonchimeras.fasta \
--id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --relabel OTU_ \
--centroids $pfx.derep.nonchimeras.otus.fasta \
--otutabout $pfx.derep.nonchimeras.otutab.txt

echo
echo nonchimeras annotation at 97%
vsearch --threads $threads --usearch_global $pfx.derep.nonchimeras.otus.fasta --db $REF --id 0.97 \
--maxaccepts 0 --maxrejects 0 --blast6out $pfx.derep.nonchimeras.otus.tax.blast6 \
--samout $pfx.derep.nonchimeras.otus.tax.sam
