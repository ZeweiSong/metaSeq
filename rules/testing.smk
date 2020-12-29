#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:    Testing methods temporarily saved here
# Call from:      ../Snakefile
# Author:         Chao | fangchao@genomics.cn
# ===================================================================

configfile: "config.yaml"

rule TEST_mash_freq:
	input: "{sample}/mash/BI.msh"
	output: "{sample}/mash/BI.msh.freq"
	threads: 36
	shell:
		"mash info -c {input} > {output}"

rule TEST_mash_freqCombine:
	input:
		f ="{sample}/mash/BI.msh.freq",
		t = "{sample}/mash/BI.msh.tsv"
	output: "{sample}/mash/BI.msh.freq.tsv"
	threads: 1
	shell:
		"metabbq beadStat mfreq -i {input.f} -t {input.t} -o {output} -v"

rule TEST_mash_freqStat:
	input: "{sample}/mash/BI.msh.freq.tsv",
	output: "{sample}/mash/BI.msh.freq.stat"
	threads: 1
	shell:
		"awk '{{b=sprintf(\"%.0f\",$2/10)*10;c=sprintf(\"%.0f\",$3/100)*100;print b\"\\t\"c\"\\t\"$4}}' < {input} |sort|uniq -c > {output}"

rule TEST_clean_merge:
	input:
		f1 = "{sample}/clean/fastp.sort.1.fq",
		f2 = "{sample}/clean/fastp.sort.2.fq"
	output:
		fa = "{sample}/clean/fastp.sort.merged.fa",
		log ="{sample}/clean/fastp.sort.merged.fa.log"
	threads: 24
	shell:
		"vsearch --threads {threads} --fastq_mergepairs {input.f1} --reverse {input.f2} "
		"--fasta_width 0 --fastaout {output.fa} &> {output.log}"

rule TEST_clean_clust:
	input:
		fa = "{sample}/clean/fastp.sort.merged.fa"
	output:
		fa0="{sample}/clean/merged.derep.fasta",
		uc0="{sample}/clean/merged.derep.uc",
		fa1="{sample}/clean/merged.preclust.fasta",
		uc1="{sample}/clean/merged.preclust.uc",
		fa2="{sample}/clean/merged.nonchime.fasta",
		uc2="{sample}/clean/merged.nonchime.uc"
	params:
		pct = config['p_VS_clust_Sid'],
	threads: 24
	shell:
		"vsearch --threads {threads} --derep_fulllength {input} "
		"--id .99 --strand both --fasta_width 0 --minuniquesize 1 "
		"--output {output.fa0} -uc {output.uc0}\n"
		"vsearch --threads {threads} --cluster_fast {output.fa0} "
		"--id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
		"--centroids {output.fa1} -uc {output.uc1}\n"
		"vsearch --threads {threads} --uchime_denovo {output.fa1} "
		"--id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
		"--relabel SSU_ --relabel_keep "
		"--nonchimeras {output.fa2} -uc {output.uc2}"


rule TEST_clean_merge2nt:
	input:
		fa = "{sample}/clean/fastp.sort.merged.fa"
	output:
		m6 = "{sample}/clean/fastp.sort.merged.nt.m6"
	log:	"{sample}/clean/fastp.sort.merged.log"
	threads: 16
	shell:
		"blastn -num_threads {threads} -db $MDB/NCBI/blast_nt/nt -query {input.fa} "
		"-outfmt '6 std qlen staxid ssciname' -out {output.m6}"


rule TEST_asm_clip2nt:
	input:
		fa = "{sample}/summary.BI.megahit.clip.fasta"
	output:
		m6 = "{sample}/summary.BI.megahit.clip.nt.m6"
	threads: 32
	shell:
		"blastn -num_threads {threads} -db $MDB/NCBI/blast_nt/nt -query {input.fa} "
		"-outfmt '6 std qlen staxid ssciname' -out {output.m6}"

rule TEST_asm_clip_heads:
	input:
		fa = "{sample}/summary.BI.megahit.clip.fasta"
	output: "{sample}/summary.BI.megahit.clip.heads"
	threads: 16
	shell:
		"grep '>' {input}|sed 's/>//;s/_/\\t/g;s/multi=//;s/len=//' > {output}"

rule TEST_asm_clip2zymo:
	input:
		fa = "{sample}/summary.BI.megahit.clip.all.fasta"
	output: "{sample}/summary.BI.megahit.clip.zymo.bam"
	threads: 16
	shell:
		"bwa mem -t {threads} ../../Source/REF/zymo/D6305.genomes.bwa {input} "
		"|samtools view -b - > {output}\n"

rule TEST_asm_clip2zymo_check:
	input: "{sample}/summary.BI.megahit.clip.zymo.bam"
	output: "{sample}/summary.BI.megahit.clip.zymo.region"
	shell:
		"perl ../04.test_vsearch/checkRegion.pl ../../Results/ZYMO/D6305.ssrRNA2genome.bed {input} {output}"

rule TEST_asm_clip_kmerStat:
	input: "{sample}/summary.BI.megahit.clip.fasta"
	output: "{sample}/summary.BI.megahit.clip.kmerStat"
	shell:
		"metabbq kmerStat -n 30 -k 19 -i {input} -o {output} "


rule TEST_asm_rna_heads:
	input:
		fa = "{sample}/summary.BI.megahit.rRNA.fasta"
	output: "{sample}/summary.BI.megahit.rRNA.heads"
	threads: 16
	shell:
		"grep '>' {input}|sed 's/>//;s/[:_]/\t/g;s/multi=//;s/len=//'|cut -f 1,4-10,12-13 > {output}"


rule TEST_asm_rRNA2nt:
	input:
		fa = "{sample}/summary.BI.megahit.rRNA.fasta"
	output:
		m6 = "{sample}/summary.BI.megahit.rRNA.nt.m6"
	threads: 16
	shell:
		"blastn -num_threads {threads} -db $MDB/NCBI/blast_nt/nt -query {input.fa} "
		"-outfmt '6 std qlen staxid ssciname' -out {output.m6}"

rule TEST_asm_rRNA2slv:
	input:
		fa = "{sample}/summary.BI.megahit.rRNA.fasta"
	output:
		m6 = "{sample}/summary.BI.megahit.rRNA.nt.m6"
	threads: 16
	shell:
		"blastn -num_threads {threads} -db $MDB/NCBI/blast_nt/nt -query {input.fa} "
		"-outfmt '6 std qlen staxid ssciname' -out {output.m6}"


rule TEST_clip_cluster:
	input: "{sample}/summary.BI.megahit.clip.fasta"
	output:
		fa1="{sample}/VSEARCH/clip.preclust.fasta",
		uc1="{sample}/VSEARCH/clip.preclust.uc",
		fa2="{sample}/VSEARCH/clip.nonchime.fasta"
	params:
		pct = config['p_VS_clust_Sid'],
	threads: 24
	shell:
		"vsearch --threads {threads} --cluster_fast {input} "
		"--id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
		"--centroids {output.fa1} -uc {output.uc1}\n"
		"vsearch --threads {threads} --uchime_denovo {output.fa1} "
		"--iddef 0 --id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
		"--relabel LFR_ --relabel_keep "
		"--nonchimeras {output.fa2}"


rule TEST_clip_kmer_bwa:
	input: "{sample}/KMER/clip.kmer.fa"
	output:
		m6="{sample}/KMER/clip.kmer.zymo.genome.m6",
		ss="{sample}/KMER/clip.kmer.zymo.genome.region"
	threads: 24
	shell:
		"blastn -num_threads {threads} -query {input} -db ../../Source/REF/zymo/D6305.genomes.blast "
		"-out {output} -outfmt 6 -word_size 7 -evalue 10\n"
		"perl ../04.test_vsearch/checkRegion.pl ../../Results/ZYMO/D6305.ssrRNA2genome.bed {output.m6} {output.stat}"

rule TEST_clip_kmer_blast:
	input: "{sample}/KMER/clip.kmer.fa"
	output:
		m6="{sample}/KMER/clip.kmer.zymo.genome.bam"
	threads: 24
	shell:
		"bwa mem -t {threads} -a -db ../../Source/REF/zymo/D6305.genomes.bwa {input} "
		"|samtools view -b - > {output}"


rule TEST_clip_rRNA_combineCluster:
	input:
		clip = "{sample}/summary.BI.megahit.clip.fasta",
		rrna = "{sample}/summary.BI.megahit.rRNA.fasta"
	output:
		fa0="{sample}/VSEARCH/comb.clip_rRNA.fa",
		fa1="{sample}/VSEARCH/comb.clip_rRNA.preclust.fa",
		uc1="{sample}/VSEARCH/comb.clip_rRNA.preclust.uc",
		fa2="{sample}/VSEARCH/comb.clip_rRNA.nonchime.fa"
	params:
		pct = config['p_VS_clust_Sid'],
	threads: 24
	shell:
		"metabbq comb2files.pl {input.rrna} {input.clip} {output.fa0}\n"
		"vsearch --threads {threads} --cluster_fast {output.fa0} "
		"--id .95 --strand both --fasta_width 0 --minuniquesize 1 "
		"--centroids {output.fa1} -uc {output.uc1}\n"
		"vsearch --threads {threads} --uchime_denovo {output.fa1} "
		"--id .95 --strand both --fasta_width 0 --minuniquesize 1 "
		"--relabel LFR_ --relabel_keep "
		"--nonchimeras {output.fa2}"


rule TEST_combine_c5:
	input: "{sample}/VSEARCH/comb.clip_rRNA.preclust.fa"
	output:"{sample}/VSEARCH/comb.clip_rRNA.nonchime.c5.uc"
	threads: 24
	shell:
		"vsearch --threads {threads} --cluster_fast {input} "
		"--id .5 --strand both --minuniquesize 10 "
		"--uc {output}\n"

rule TEST_combine_mafft:
	input: "{sample}/VSEARCH/comb.clip_rRNA.nonchime.fa"
	output:"{sample}/VSEARCH/comb.clip_rRNA.nonchime.mafft.fa"
	threads: 24
	shell:
		"mafft --thread {threads} --reorder {input} > {output}"


rule TEST_LFR_mafft:
	input: "{sample}/VSEARCH/barrnap.LFRs.fasta"
	output:
		ssu = "{sample}/VSEARCH/barrnap.LFRs.S.fa",
		lsu = "{sample}/VSEARCH/barrnap.LFRs.L.fa",
		sma = "{sample}/VSEARCH/barrnap.LFRs.S.mafft.fa",
		lma = "{sample}/VSEARCH/barrnap.LFRs.L.mafft.fa"
	threads: 24
	shell:
		"grep -A1 16S --no-group-separator {input} > {output.ssu}\n"
		"grep -A1 23S --no-group-separator {input} > {output.lsu}\n"
		"mafft --thread {threads} --reorder {output.ssu} > {output.sma}\n"
		"mafft --thread {threads} --reorder {output.lsu} > {output.lma}\n"

rule TEST_clipS_silva:
	input:
		ssu = "{sample}/summary.BI.megahit.clip.fasta",
	output:
		ssu = "{sample}/summary.BI.megahit.clip.SSU.m6",
		ann = "{sample}/summary.BI.megahit.clip.SSU.m6.more",
	threads: 16
	shell:
		"blastn -num_threads {threads} -perc_identity 90 "
		"-db $STL/Source/REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta "
		"-query {input.ssu} -outfmt '6 std qlen' -out {output.ssu}\n"
		"metabbq anno.pl $STL/Source/REF/silva132/SSU/taxmap_embl_ssu_ref_132.tax {output.ssu} > {output.ann}\n"

rule TEST_clipL_silva:
	input:
		lsu = "{sample}/summary.BI.megahit.clip.fasta"
	output:
		ann = "{sample}/summary.BI.megahit.clip.LSU.m6.more",
		lsu = "{sample}/summary.BI.megahit.clip.LSU.m6",
	threads: 16
	shell:
		"blastn -num_threads {threads} -perc_identity 90 "
		"-db $STL/Source/REF/silva132/LSU/SILVA_132_LSURef_tax_RNA.fasta "
		"-query {input.lsu} -outfmt '6 std qlen' -out {output.lsu} \n"
		"metabbq anno.pl $STL/Source/REF/silva132/LSU/taxmap_embl_lsu_ref_132.tax {output.lsu} > {output.lsu}.more\n"

rule TEST_ssu_silva:
	input:
		ssu = "{sample}/VSEARCH/barrnap.SSU.fasta",
	output:
		ssu = "{sample}/VSEARCH/barrnap.SSU.m6",
		ann = "{sample}/VSEARCH/barrnap.SSU.m6.more",
	threads: 8
	shell:
		"blastn -num_threads {threads} -perc_identity 90 "
		"-db $STL/Source/REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta "
		"-query {input.ssu} -outfmt '6 std qlen' -out {output.ssu}\n"
		"metabbq anno.pl $STL/Source/REF/silva132/SSU/*txt {output.ssu} > {output.ann}\n"

rule TEST_lsu_silva:
	input:
		lsu = "{sample}/VSEARCH/barrnap.LSU.fasta"
	output:
		ann = "{sample}/VSEARCH/barrnap.LSU.m6.more",
		lsu = "{sample}/VSEARCH/barrnap.LSU.m6",
	threads: 8
	shell:
		"blastn -num_threads {threads} -perc_identity 90 "
		"-db $STL/Source/REF/silva132/LSU/SILVA_132_LSURef_tax_RNA.fasta "
		"-query {input.lsu} -outfmt '6 std qlen' -out {output.lsu} \n"
		"metabbq anno.pl $STL/Source/REF/silva132/LSU/*txt {output.lsu} > {output.lsu}.more\n"
rule TEST_cdhitS_silva:
	input:
		ssu = "{sample}/VSEARCH/barrnap.cdhitS.fasta",
	output:
		ssu = "{sample}/VSEARCH/barrnap.cdhitS.m6",
		ann = "{sample}/VSEARCH/barrnap.cdhitS.m6.more",
	threads: 8
	shell:
		"blastn -num_threads {threads} -perc_identity 90 "
		"-db $STL/Source/REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta "
		"-query {input.ssu} -outfmt '6 std qlen' -out {output.ssu}\n"
		"metabbq anno.pl $STL/Source/REF/silva132/SSU/*txt {output.ssu} > {output.ann}\n"

rule TEST_cdhitL_silva:
	input:
		lsu = "{sample}/VSEARCH/barrnap.cdhitL.fasta"
	output:
		ann = "{sample}/VSEARCH/barrnap.cdhitL.m6.more",
		lsu = "{sample}/VSEARCH/barrnap.cdhitL.m6",
	threads: 8
	shell:
		"blastn -num_threads {threads} -perc_identity 90 "
		"-db $STL/Source/REF/silva132/LSU/SILVA_132_LSURef_tax_RNA.fasta "
		"-query {input.lsu} -outfmt '6 std qlen' -out {output.lsu} \n"
		"metabbq anno.pl $STL/Source/REF/silva132/LSU/*txt {output.lsu} > {output.lsu}.more\n"


rule TEST_BB_hybridize:
	input: "{sample}/summary.BI.megahit.clip.fasta"
	output: "{sample}/summary.BI.megahit.clip.num"
	shell:
		"grep '^>' {input} > {ouput}"

rule TEST_BB_mashstat:
	input: "{sample}/mash/bMin2.msh"
	output: "{sample}/mash/bMin2.msh.tsv"
	shell:
		"mash info -t {input} > {output}"


### stat
rule STAT_barrnap_fa:
	input:
		ssu = "{sample}/VSEARCH/barrnap.SSU.fasta",
		lsu = "{sample}/VSEARCH/barrnap.LSU.fasta"
	output: "{sample}/stat/barrnap.RSU.stat"
	threads: 1
	shell:
		"metabbq stat RSU -l {input.lsu} -s {input.ssu} -o {output}"

### stat
rule STAT_callSNV_clip_2_mock_ref:
	input: "{sample}/summary.BI.megahit.clip.all.fasta"
	output:
		bam = "{sample}/CLIP/clips2REF.bam",
		bcf = "{sample}/CLIP/clips2REF.bcf",
		tsv = "{sample}/CLIP/clips2REF.bcf.tsv"
	threads: 8
	params: config['MockRef']['bwa_db']
	shell:
		"awk 'FNR%2==1{{if($0~/flag=[123]/){{p=1}}else{{p=0}}}}(p==1){{print}}' {input}|bwa mem -k 141 -t {threads} {params} - | samtools sort -O bam -o {output.bam} - && "
		"bcftools mpileup --threads {threads} -Ob -o {output.bcf} -f {params} {output.bam} && "
		"bcftools view {output.bcf}|metabbq readBCF2tsv > {output.tsv}"

rule STAT_clip_2_mock_ref:
	input: "{sample}/summary.BI.megahit.clip.all.fasta"
	output:
		m6 = "{sample}/ANNO/clip2MockRef.m6"
	threads: 8
	params: config['MockRef']['blastdb']
	shell:
		"blastn -num_threads {threads} -perc_identity 95 -word_size 77 -db {params} "
		"-query {input} -outfmt '6 std qlen' -out {output.m6}\n"

rule STAT_clip_2_known_bb_stat:
	input: "{sample}/ANNO/clip2MockRef.m6"
	output:
		b= "{sample}/ANNO/clip2MockRef.bead.anno",
		c= "{sample}/ANNO/clip2MockRef.clip.anno",
		g= "{sample}/ANNO/clip2MockRef.ref.cov"
	params: config['MockRef']['bed']
	shell:
		"metabbq checkBeadRegion -i {input} -r {params} -o {output.b} -c {output.c} -g {output.g} \n"


### stat LOTU
rule STAT_lotu_2_mock_ref:
	input: "{sample}/CLIP/LOTU.fa"
	output:
		m6 = "{sample}/ANNO/LOTU.to.MockRef.m6",
		m0 = "{sample}/ANNO/LOTU.to.MockRef.m0"
	threads: 4
	params: config['MockRef']['blastdb']
	shell:
		"blastn -num_threads {threads} -perc_identity 95 -word_size 141 -db {params} "
		"-query {input} -outfmt '6 std qlen slen' -out {output.m6} &\n"
		"blastn -num_threads {threads} -perc_identity 95 -word_size 141 -db {params} "
		"-query {input} -out {output.m0}\n"

rule STAT_lotu_2_known_bb_stat:
	input: "{sample}/ANNO/LOTU.to.MockRef.m6"
	output:
		b= "{sample}/ANNO/LOTU.to.MockRef.bead.anno",
		c= "{sample}/ANNO/LOTU.to.MockRef.clip.anno",
		g= "{sample}/ANNO/LOTU.to.MockRef.ref.cov"
	params: config['MockRef']['bed']
	shell:
		"metabbq checkBeadRegion -i {input} -r {params} -o {output.b} -c {output.c} -g {output.g} \n"

rule STAT_lotu_2_mock_quast:
	input: "{sample}/CLIP/LOTU.fa"
	output:
		m6 = "{sample}/CLIP/LOTU.to.MockRef.quast"
	threads: 8
	params: config['MockRef']['blastdb']
	shell:
		"quast.py -t {threads} {input} -r {params} -o {output}"

###################### Temp ####################################################
rule STAT_clip_2_known_fix_ref:
	input: "{sample}/summary.BI.megahit.clip.fasta"
	output:
		m6 = "{sample}/summary.BI.megahit.clip2MockRef.fix.m6"
	threads: 8
	params: "$LFR/Source/REF/zymo/D6305.fix.rRNA.fa"
	shell:
		"blastn -num_threads {threads} -perc_identity 90 -db {params} "
		"-query {input} -outfmt '6 std qlen' -out {output.m6}\n"

rule STAT_clip_2_known_fix_bb_stat:
	input: "{sample}/summary.BI.megahit.clip2MockRef.fix.m6"
	output:
		b= "{sample}/CLIP/clip2MockRef.fix.bead.anno",
		c= "{sample}/CLIP/clip2MockRef.fix.clip.anno",
		g= "{sample}/CLIP/clip2MockRef.fix.ref.cov"
	params: "$LFR/Source/REF/zymo/D6305.fix.rRNA.barrnap.bed"
	shell:
		"metabbq checkBeadRegion -i {input} -r {params} -o {output.b} -c {output.c} -g {output.g} \n"
###################### Temp ####################################################

rule TEST_clip_2_unite_blast:
	input: "{sample}/summary.BI.megahit.clip.fasta",
	output:
		unite = "{sample}/summary.BI.megahit.clip.UNITE.m6",
		ann = "{sample}/summary.BI.megahit.clip.UNITE.m6.more",
	threads: 8
	shell:
		"blastn -num_threads {threads} -perc_identity 90 "
		"-db $LFR/Source/REF/UNITE/sh_general_release_all_04.02.2020/sh_general_release_dynamic_all "
		"-query {input} -outfmt '6 std qlen' -out {output.unite}\n"
		"metabbq anno.pl $LFR/Source/REF/UNITE/sh_general_release_all_04.02.2020/sh_general_release_dynamic_all.tax {output.unite} > {output.ann}\n"


rule STAT_clip_2_slv_blast:
	input:
		ssu = "{sample}/summary.BI.megahit.clip.SSU.m6.more",
		lsu = "{sample}/summary.BI.megahit.clip.LSU.m6.more"
	output:
		clip = "{sample}/stat/summary.BI.megahit.clip.slv.stat.clips",
		bead = "{sample}/stat/summary.BI.megahit.clip.slv.stat.beads"
	threads: 1
	shell:
		"metabbq stat slv -l {input.lsu} -s {input.ssu} -o {output.clip} -g {output.bead}\n"

rule STAT_clip_2_slv_getTaxRank:
	input:
		clip = "{sample}/stat/summary.BI.megahit.clip.slv.stat.clips",
		bead = "{sample}/stat/summary.BI.megahit.clip.slv.stat.beads"
	output:
		clip = "{sample}/stat/summary.BI.megahit.clip.slv.stat.clips.tax",
		bead = "{sample}/stat/summary.BI.megahit.clip.slv.stat.beads.tax"
	threads: 4
	shell:
		"awk -F \"\\t\" '{{gsub(/ <.*>/,\"\",$4);print}}' {input.bead} "
		"| taxonkit name2taxid -r -i 4 > {output.bead} &\n"
		"awk -F \"\\t\" '{{gsub(/ <.*>/,\"\",$7);print $1\"\\t\"$7}}' {input.clip} "
		"| taxonkit name2taxid -r -i 2 > {output.clip}"

# From clip
rule CLIP_1_cluster:
	input: "{sample}/summary.BI.megahit.clip.fasta"
	output:
		fa1="{sample}/CLIP/preclust.fa",
		uc1="{sample}/CLIP/preclust.uc",
		log="{sample}/CLIP/cluster.log"
	params:
		pct = 0.99,
	threads: config['thread']['vsearch']
	shell:
		"vsearch --threads {threads} --cluster_fast {input} --iddef 4 --id {params.pct} "
		"--strand both --fasta_width 0 --centroids {output.fa1} -uc {output.uc1} --sizeout &> {output.log} \n"

rule CLIP_2_filter:
	input: "{sample}/CLIP/preclust.fa"
	output:
		fa1="{sample}/CLIP/preclust_1k.fa",
	params:
		len = 999,
		minSize= 2,
	threads: config['thread']['vsearch']
	shell:
		"awk -F '_|;size=' '/>/{{if($9>{params.len}&&$10>{params.minSize}){{p=1}}else{{p=0}}}} (p==1){{print}}' "
		"{input} > {output}\n"

if config["sampleType"] == "F":
	rule CLIP_3_predict_barranp:
		input: outXac
		output:
			fa2="{sample}/CLIP/all.barrnap.fa",
			gff="{sample}/CLIP/all.barrnap.gff",
			log="{sample}/CLIP/all.barrnap.log",
		params:
			pfx = "{sample}/CLIP/all.ITS",
		threads: 8
		shell:
			"barrnap --threads {threads} --kingdom euk --outseq {output.fa2} < {input} > {output.gff} 2> {output.log}\n"
	rule CLIP_4_sumgff:
		input: "{sample}/CLIP/all.barrnap.gff"
		output: "{sample}/CLIP/all.barrnap.sum.tsv"
		threads: 1
		shell:
			"metabbq IO sumgff -d euk -i {input} -o {output}"
	rule CLIP_3_predict_ITSx:
		input: outXac
		output:
			log="{sample}/CLIP/all.ITS.log",
			pos="{sample}/CLIP/all.ITS.positions.txt",
		params:
			pfx = "{sample}/CLIP/all.ITS",
		threads: 8
		shell:
			"ITSx --cpu {threads} --save_regions all -i {input} -o {params.pfx} &> {output.log}\n"
	rule CLIP_4_sumITSxPos:
		input: "{sample}/CLIP/all.ITS.positions.txt"
		output: "{sample}/CLIP/all.ITS.positions.sum.tsv"
		threads: 1
		shell:
			"metabbq IO itsx -i {input} -o {output}"
	rule PRED_0_getSSU:
		input: "{sample}/CLIP/all.barrnap.fa"
		output:"{sample}/PREDICT/barrnap.ssu.fa"
		shell:
			"metabbq IO printUniqSubunit -i {input} -t 18S_rRNA -o {output} -v"
	rule PRED_0_getLSU:
		input: "{sample}/CLIP/all.barrnap.fa"
		output:"{sample}/PREDICT/barrnap.lsu.fa"
		shell:
			"metabbq IO printUniqSubunit -i {input} -t 28S_rRNA -o {output} -v"
	rule PRED_0_getITSxSubunits:
		input:
			inf = outXcafa,
			pos = "{sample}/CLIP/all.ITS.positions.txt"
		output:
			SSU = "{sample}/CLIP/all.ITS.SSU.fasta",
			LSU = "{sample}/CLIP/all.ITS.LSU.fasta"
		params: "{sample}/CLIP/all.ITS"
		shell:
			"metabbq IO printITSxSubunit -i {input.inf} -p {input.pos} -o {params} -v"

else:
	rule CLIP_3_predict:
		input: outXac
		output:
			fa="{sample}/CLIP/all.barrnap.fa",
			gff="{sample}/CLIP/all.barrnap.gff",
			log="{sample}/CLIP/all.barrnap.log"
		threads: 8
		shell:
			"barrnap --threads {threads} --kingdom bac --outseq {output.fa} < {input} > {output.gff} 2> {output.log}\n"
	rule CLIP_4_sumgff:
		input: "{sample}/CLIP/all.barrnap.gff"
		output: "{sample}/CLIP/all.barrnap.sum.tsv"
		threads: 1
		shell:
			"metabbq IO sumgff -i {input} -o {output}"
	rule PRED_0_getSSU:
		input: "{sample}/CLIP/all.barrnap.fa"
		output:"{sample}/PREDICT/barrnap.ssu.fa"
		shell:
			"metabbq IO printUniqSubunit -i {input} -t 16S_rRNA -o {output} -v"
	rule PRED_1_ClustSSU:
		input:"{sample}/PREDICT/barrnap.ssu.fa"
		output:"{sample}/PREDICT/barrnap.ssu.clade.uc_0.995"
		params:"{sample}/PREDICT/barrnap.ssu.clade.uc"
		shell:
			"metabbq IO makeclade -i {input} -o {params} -v"
	rule PRED_0_getLSU:
		input: "{sample}/CLIP/all.barrnap.fa"
		output:"{sample}/PREDICT/barrnap.lsu.fa"
		shell:
			"metabbq IO printUniqSubunit -i {input} -t 23S_rRNA -o {output} -v"
	rule PRED_1_ClustLSU:
		input:"{sample}/PREDICT/barrnap.lsu.fa"
		output:"{sample}/PREDICT/barrnap.lsu.clade.uc_0.995"
		params:"{sample}/PREDICT/barrnap.lsu.clade.uc"
		shell:
			"metabbq IO makeclade -i {input} -o {params} -v"

#test SOTU mapping to refs
rule zymo_sotu_2_MockRef:
    input:
        SSU = "{sample}/OTU/pred.SSU.clust.fa",
        LSU = "{sample}/OTU/pred.LSU.clust.fa"
    output:
        SSU = "{sample}/OTU/pred.SSU.clust.zymo.m6",
        LSU = "{sample}/OTU/pred.LSU.clust.zymo.m6"
    params:
        LSU = "$LFR/Source/REF/zymo/D6305.rRNA.fa",
        SSU = "$LFR/Source/REF/zymo/D6305.ssrRNA.fasta"
    threads: 8
    shell:
        "blastn -num_threads {threads} -db {params.LSU} -query {input.LSU} -outfmt '6 std qlen staxid ssciname' -out {output.LSU} &\n"
        "blastn -num_threads {threads} -db {params.SSU} -query {input.SSU} -outfmt '6 std qlen staxid ssciname' -out {output.SSU} \n"



##LOTU test mapping to close Ref

if config["sampleType"] == "F":
    kingdom="euk"
    LSU="28S"
    SSU="18S"
else:
    kingdom="bac"
    LSUDB="$STL/Source/REF/silva132/LSU/SILVA_132_LSURef_tax_RNA.fasta"
    LSUAN="$STL/Source/REF/silva132/LSU/taxmap_embl_lsu_ref_132.tax"
    SSU="16S"


rule LOTU_closeDB_blastn:
    input: "{sample}/CLIP/LOTU.fa"
    output: "{sample}/ANNO/LOTU.map.MockRef.m6"
    params: config["MockRef"]["bwa_db"]
    threads: 8
    shell:
        "blastn -num_threads {threads} -perc_identity 95 -db {params} "
        " -query {input} -outfmt '6 std qlen slen' -out {output}\n"

rule LOTU_closeDB_bwa:
    input: "{sample}/CLIP/LOTU.fa"
    output: "{sample}/ANNO/LOTU.map.MockRef.bam"
    params: config["MockRef"]["bwa_db"]
    threads: 8
    shell:
        "bwa mem -t {threads} {params} {input} "
        "|samtools view -b - > {output}\n"


## CLEAN: randomly downsize:
rule downsize_50_10:
    input:
        fq1="{sample}/clean/fastp.sort.1.fq.gz",
        fq2="{sample}/clean/fastp.sort.2.fq.gz"
    output:
        fq1="{sample}/clean/fastp.sort.1_00.fq.gz",
        fq2="{sample}/clean/fastp.sort.2_00.fq.gz"
    params:
        pfx1 = "{sample}/clean/fastp.sort.1",
        pfx2 = "{sample}/clean/fastp.sort.2",
        p = config["clean_ds_pieces"],
        k  = config["clean_ds_keep_files"]
    threads: 4
    shell:
        "metabbq randomlyPick.pl -i {input.fq1} -o {params.pfx1} -n {params.p} -m {params.keeps} -s {params.p} -t 2 &\n"
        "metabbq randomlyPick.pl -i {input.fq2} -o {params.pfx2} -n {params.p} -m {params.keeps} -s {params.p} -t 2 -v"


rule test_LOTU_i95:
    input:"{sample}/CLIP/id95def4.clust.uc"
    output:"{sample}/CLIP/id95def4.clust.H.BID"
    shell:
        """
awk '(/^H/&&$3>999){{split($9,a,"_");split($10,b,"_");print a[1];print b[1]}}' {input}|sort|uniq > {output}
        """


###############################################################################
# Avoid hybridize beads by clustering clips
###############################################################################

rule testKK2_2_listClustBID:
    input:
        fa = outXac,
        uc = "{sample}/CLIP/id95def4.clust.uc",
        anno="{sample}/ANNO/clip2MockRef.bead.anno"
    output: "{sample}/CLIP/mock.id95def4.clust.fa"
    threads: 1
    shell:
        "metabbq IO uch -a {input.anno} -i {input.fa} -u {input.uc} -o {output} -c 999"

rule testKK2_3_ClipClust995:
    input: "{sample}/CLIP/mock.id95def4.clust.fa"
    output:
        uc = "{sample}/CLIP/mock.id995def2.clust.uc",
        ofa = "{sample}/CLIP/mock.id995def2.clust.fa",
        cfa = "{sample}/KRAKEN/mock/library/added.fna"
    threads: 8
    shell:
        "vsearch --threads {threads} --cluster_fast {input} --strand both --iddef 2 --id .995 "
        " --fasta_width 0 --relabel LOTU_ --uc {output.uc}  --centroids {output.ofa} \n"
        "cp -r {output.ofa} {output.cfa}"



if config["sampleType"] == "F":
    rule testKK2_4_makeTaxonomyTree:
        input:
            uc = "{sample}/CLIP/mock.id995def2.clust.uc",
            anno="{sample}/ANNO/clip2MockRef.bead.anno"
        output: "{sample}/KRAKEN/mock/data/added.txt"
        params: "Eukaryota"
        shell:
            "metabbq IO maketaxon -d {params} -u {input.uc} -a {input.anno} -t $LFR/Source/REF/KRAKEN2/rDNA/data/mergeTaxon.txt -o {output}"
else:
    rule testKK2_4_makeTaxonomyTree:
        input:
            uc = "{sample}/CLIP/mock.id995def2.clust.uc",
            anno="{sample}/ANNO/clip2MockRef.bead.anno"
        output: "{sample}/KRAKEN/mock/data/added.txt"
        params: "Bacteria"
        shell:
            "metabbq IO maketaxon -u {input.uc} -a {input.anno} -t $LFR/Source/REF/KRAKEN2/rDNA/data/mergeTaxon.txt -o {output}"


rule testKK2_5_makedb:
    input:
        fa = "{sample}/KRAKEN/mock/library/added.fna",
        ad = "{sample}/KRAKEN/mock/data/added.txt"
    output: "{sample}/KRAKEN/mock/hash.k2d"
    params: "{sample}/KRAKEN/mock"
    threads: 4
    shell:
        """
if [[ -e {params}/data/mergeTaxon.txt ]]; then rm {params}/data/mergeTaxon.txt ;fi
cp -r $LFR/Source/REF/KRAKEN2/rDNA/data/mergeTaxon.txt {params}/data/mergeTaxon.txt
metabbq buildKraken2db.sh {params} {threads}
        """

rule testKK2_6_mapping:
    input:
        fa1 = "{sample}/clean/fastp.sort.1.fq.gz",
        fa2 = "{sample}/clean/fastp.sort.2.fq.gz",
        db  = "{sample}/KRAKEN/mock/hash.k2d"
    output:
        raw = "{sample}/KRAKEN/reads2mock.out",
        prf = "{sample}/KRAKEN/reads2mock.prof"
    params: "{sample}/KRAKEN/mock"
    threads: 4
    shell:
        "kraken2 --db {params} --threads {threads} --paired {input.fa1} {input.fa2} "
        " | tee {output.raw} | metabbq IO kk2prf -d {params} -i - -o {output.prf}"

# rule testKK2_6_profling:
#     input:
#         kk2 = "{sample}/testKK2/reads2db.out",
#         db  = "{sample}/testKK2/db/hash.k2d"
#     output: "{sample}/testKK2/reads2db.prof"
#     params: "{sample}/testKK2/db"
#     threads: 1
#     shell:
#         "metabbq IO kk2prf -d {params} -i {input.kk2} -o {output} -v"

rule testKK2_7_sum2species:
    input:"{sample}/KRAKEN/reads2mock.prof"
    output: "{sample}/KRAKEN/reads2mock.sp.prf"
    threads: 1
    shell:
        "metabbq IO sumprf -r species -i {input} -o {output}"
