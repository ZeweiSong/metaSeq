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
		fa = "{sample}/summary.BI.megahit.clip.fasta"
	output: "{sample}/summary.BI.megahit.clip.zymo.bam"
	threads: 16
	shell:
		"bwa mem -t {threads} -db ../../Source/REF/zymo/D6305.genomes.bwa {input} "
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
rule STAT_clip_2_known_ref:
	input: "{sample}/summary.BI.megahit.clip.fasta"
	output:
		m6 = "{sample}/summary.BI.megahit.clip2CloseRef.m6"
	threads: 8
	params: config['closeRef']['blastdb']
	shell:
		"blastn -num_threads {threads} -perc_identity 90 -db {params} "
		"-query {input} -outfmt '6 std qlen' -out {output.m6}\n"

rule STAT_clip_2_known_bb_stat:
	input: "{sample}/summary.BI.megahit.clip2CloseRef.m6"
	output:
		b= "{sample}/CLIP/clip2CloseRef.bead.anno",
		c= "{sample}/CLIP/clip2CloseRef.clip.anno",
		g= "{sample}/CLIP/clip2CloseRef.ref.cov"
	params: config['closeRef']['bed']
	shell:
		"metabbq checkBeadRegion -i {input} -r {params} -o {output.b} -c {output.c} -g {output.g} \n"

###################### Temp ####################################################
rule STAT_clip_2_known_fix_ref:
	input: "{sample}/summary.BI.megahit.clip.fasta"
	output:
		m6 = "{sample}/summary.BI.megahit.clip2CloseRef.fix.m6"
	threads: 8
	params: "$LFR/Source/REF/zymo/D6305.fix.rRNA.fa"
	shell:
		"blastn -num_threads {threads} -perc_identity 90 -db {params} "
		"-query {input} -outfmt '6 std qlen' -out {output.m6}\n"

rule STAT_clip_2_known_fix_bb_stat:
	input: "{sample}/summary.BI.megahit.clip2CloseRef.fix.m6"
	output:
		b= "{sample}/CLIP/clip2CloseRef.fix.bead.anno",
		c= "{sample}/CLIP/clip2CloseRef.fix.clip.anno",
		g= "{sample}/CLIP/clip2CloseRef.fix.ref.cov"
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
	rule CLIP_3_predict:
		input: "{sample}/CLIP/preclust_1k.fa"
		output:
			fa1="{sample}/CLIP/preclust_1k.ITS.positions.txt",
			fa2="{sample}/CLIP/preclust_1k.barrnap.fa",
			gff="{sample}/CLIP/preclust_1k.barrnap.gff",
			log="{sample}/CLIP/preclust_1k.barrnap.log",
		params:
			pfx = "{sample}/CLIP/preclust_1k.ITS",
		threads: 2
		shell:
			"ITSx -t F --cpu {threads} -i {input} -o {params.pfx} &> {output.fa1} &\n"
			"barrnap --threads {threads} --kingdom euk --outseq {output.fa2} < {input} > {output.gff} 2> {output.log}\n"

else:
	rule CLIP_3_predict:
		input: "{sample}/CLIP/preclust_1k.fa"
		output:
			fa="{sample}/CLIP/preclust_1k.barrnap.fa",
			gff="{sample}/CLIP/preclust_1k.barrnap.gff",
			log="{sample}/CLIP/preclust_1k.barrnap.log"
		threads: 2
		shell:
			"barrnap --threads {threads} --kingdom bac --outseq {output.fa} < {input} > {output.gff} 2> {output.log}\n"


#test SOTU mapping to refs
rule zymo_sotu_2_closeRef:
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


rule LOTU_closeDB_mapping:
    input: "{sample}/CLIP/all.LOTU.fa"
    output: "{sample}/CLIP/all.LOTU.map.closeRef.m6"
    params: config["closeRef"]["bwa_db"]
    threads: 8
    shell:
        "blastn -num_threads {threads} -perc_identity 95 -db {params} "
        " -query {input} -outfmt '6 std qlen slen' -out {output}\n"



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
        pfx2 = "{sample}/clean/fastp.sort.2"
        p = config["clean_ds_pieces"],
        k  = config["clean_ds_keep_files"]
    threads: 4
    shell:
        "metabbq randomlyPick.pl -i {input.fq1} -o {params.pfx1} -n {params.p} -m {params.keeps} -s {params.p} -t 2 &\n"
        "metabbq randomlyPick.pl -i {input.fq2} -o {params.pfx2} -n {params.p} -m {params.keeps} -s {params.p} -t 2 -v"
