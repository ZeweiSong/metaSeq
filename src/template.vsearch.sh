# practice refere to :https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
threads=$1
REF=$2
sam=$3
bead=$4
tag=$5

b1=${bead:0:4}

echo VSEARCH start:
echo info: Sample: "$sam | bead = $bead | REF = $REF | tag = $tag"
echo info: output dir: "VSEARCH/$sam/$b1/$bead/"

pfx="VSEARCH/$sam/$b1/$bead/$tag"

echo
echo init and decmpress fastq files
mkdir -p VSEARCH/$sam/$b1/$bead/
pigz -dc beadPool/$sam/$b1/$bead.1.fq.gz > beadPool/$sam/$b1/$bead.1.fq
pigz -dc beadPool/$sam/$b1/$bead.2.fq.gz > beadPool/$sam/$b1/$bead.2.fq

echo
echo Merge seq
vsearch --threads $threads --fastq_mergepairs beadPool/$sam/$b1/$bead.1.fq \
--reverse beadPool/$sam/$b1/$bead.2.fq \
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
echo Cluster nonchimeras at 97% and relabel with OTU_n, generate OTU table
vsearch --threads $threads --cluster_size  $pfx.merge.derep.nonchimeras.fasta \
--id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --relabel OTU_ \
--centroids $pfx.merge.derep.nonchimeras.otus.fasta \
--otutabout $pfx.merge.derep.nonchimeras.otutab.txt

echo
echo nonchimeras annotation at 97%
vsearch --threads $threads --usearch_global $pfx.merge.derep.nonchimeras.otus.fasta --db $REF --id 0.97 \
--maxaccepts 0 --maxrejects 0 --blast6out $pfx.merge.derep.nonchimeras.otus.tax.blast6 \
--samout $pfx.merge.derep.nonchimeras.otus.tax.sam

echo "####EXTRA####"
echo "nonchimeras annotation to GENOME (test mode) at 97%"
vsearch --threads $threads --usearch_global $pfx.merge.derep.nonchimeras.otus.fasta --db REF/ncbi.5Ref.fasta --id 0.97 \
--maxaccepts 0 --maxrejects 0 --blast6out $pfx.merge.derep.nonchimeras.otus.tax.GENOME.blast6
