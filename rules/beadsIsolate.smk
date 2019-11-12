#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beads Isolate Assembly Solution
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

configfile: "config.yaml"

rule BIA_0_cutoff:
    input: "{sample}/clean/fastp.sort.1.fq.idx"
    output: "{sample}/Assemble_BI/ID.lst"
    params:
        p_rpb  = config["p_rpb_min"]
    shell:
        "awk '(!/0000/&&(($3-$2+1)/4)>{params.p_rpb}){{print FNR\"\\t\"$1\"\\t\"($3-$2+1)/4}}' {input} > {output}"

rule BIA_1_write:
    input:
        s1 = "{sample}/clean/fastp.sort.1.fq",
        s2 = "{sample}/clean/fastp.sort.2.fq",
        bi = "{sample}/Assemble_BI/ID.lst"
    output: "{sample}/Assemble_BI/log"
    params:
        outDir = "{sample}/Assemble_BI",
        threads = config["threads"]
    shell:
        "metabbq beadsWrite3.pl -d {input.bi} -f fq -t 4 -o {params.outDir} -v "
        "--r1 {input.s1}  --r2 {input.s2} &> {output}\n"

rule BIA_2_initASMsh:
    input:
        log = "{sample}/Assemble_BI/log",
        bi = "{sample}/Assemble_BI/ID.lst",
    output: "{sample}/batch.assemble.BI.sh"
    params:
        samDir = "{sample}",
        outDir = "{sample}/Assemble_BI",
        mode   = config["method"]["assemble"]["mode"],
        mem    = config["p_asm_mem"],
        cpu    = config["p_asm_cpu"]
    shell:
        "for i in `sort -k2,2nr {input.bi} | cut -f1`; do "
        "echo metabbq RCAasm.sh {params.mode} {params.samDir} BI $i BI {params.mem} {params.cpu} ; done > {output}\n"

outXs = "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.tsv"
outXa= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.fasta"
outXn= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.anno"

rule BIA_3_assemble:
    input: "{sample}/batch.assemble.BI.sh"
    output: "{sample}/log/batch.assemble.BI.done"
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        dir  = "{sample}"
    shell:
        """
# 0. Run assembly one by one
dl=$[$(wc -l {input}|tr ' ' '\\n' |head -1) / {params.divide}]
split -d -l $dl {input} {params.dpfx}
for i in {params.dpfx}??;do sh $i & done > {params.dir}/running.BI.asm.log && wait
mv {params.dir}/running.BI.asm.log {params.dir}/log/batch.assemble.BI.done
        """

rule BIA_4_summary1:
    input: "{sample}/log/batch.assemble.BI.done"
    output:
        fa   = outXa
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        outDir = "{sample}/Assemble_BI",
        mAsm   = config["p_asm_min"],
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
# 1. pick robust assemble fasta
for i in `ls {params.outDir}|grep BI`;do
  metabbq clusterHelper asmPick -l {params.mAsm} -b $i -i {params.outDir}/$i/{params.mode}/final.contigs.fa
done > {output.fa}
        """

rule BIA_4_summary2:
    input: outXa
    output:
        stat = outXs
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        outDir = "{sample}/Assemble_BI",
        mAsm   = config["p_asm_min"],
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
# 2. stat assemble results
for i in `ls {params.outDir}|grep BI`;do
  awk -v bc=$i '/^>/{{print bc"\\t"$0}}' {params.outDir}/$i/{params.mode}/final.contigs.fa | sed 's/>//;s/ flag=/\\t/;s/ multi=/\\t/;s/ len=/\\t/'
done > {output.stat}
        """

outPfx="{sample}/summary.BI." + str(config["method"]["assemble"]["mode"])
outXd= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.fasta"
outXc= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".clip.fasta"
rule BIA_5_RCACLIP:
    input: outXa
    output: outXc
    params:
        A = config['AdRCA'],
        F = config['pmrFwd'],
        R = config['pmrRev']
    shell:
        "metabbq RCACLIPS -a {params.A} -fwd {params.F} -rev {params.R} -i {input} -o {output} -v 2> {output}.log"

outRc= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".rRNA.fasta"
if config["sampleType"] == "F":
    kingdom="euk"
    LSU="28S"
    SSU="18S"
else:
    kingdom="bac"
    LSU="23S"
    SSU="16S"

rule BIA_6_rrnadetect:
    input: outXc
    output: outRc
    log: outRc + ".barrnap"
    params:
        k = kingdom,
    threads: config['thread']['vsearch']
    shell:
        "barrnap --kingdom {params.k} --threads {threads} --reject 0.1 {input} --outseq {output} &> {log}"

rule BIA_7_cut_LSU_and_SSU:
    input: outRc
    output:
        S= "{sample}/VSEARCH/barrnap.SSU.fasta",
        L= "{sample}/VSEARCH/barrnap.LSU.fasta"
    params:
        S = SSU,
        L = LSU
    threads: 2
    shell:
        "grep -A1 {params.S} {input} > {output.S} & "
        "grep -A1 {params.L} {input} > {output.L} & "
        "wait"

rule VSEARCH_1_SSU_cluster:
    input: "{sample}/VSEARCH/barrnap.SSU.fasta"
    output:
        fa1="{sample}/VSEARCH/barrnap.preclustS.fasta",
        uc1="{sample}/VSEARCH/barrnap.preclustS.uc",
        fa2="{sample}/VSEARCH/barrnap.cdhitS.fasta",
        uc2="{sample}/VSEARCH/barrnap.cdhitS.fasta.uc"
    params:
        pct = config['p_VS_clust_Sid'],
    threads: config['thread']['vsearch']
    shell:
        "vsearch --threads {threads} --cluster_fast {input} "
        "--id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
        "--centroids {output.fa1} -uc {output.uc1}\n"
        "vsearch --threads {threads} --cluster_fast {output.fa1} "
        "--iddef 0 --id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
        "--relabel SSU_ --relabel_keep "
        "--centroids {output.fa2} -uc {output.uc2}"

rule VSEARCH_1_LSU_cluster:
    input: "{sample}/VSEARCH/barrnap.LSU.fasta"
    output:
        fa1="{sample}/VSEARCH/barrnap.preclustL.fasta",
        uc1="{sample}/VSEARCH/barrnap.preclustL.uc",
        fa2="{sample}/VSEARCH/barrnap.cdhitL.fasta",
        uc2="{sample}/VSEARCH/barrnap.cdhitL.fasta.uc"
    params:
        pct = config['p_VS_clust_Lid']
    threads: config['thread']['vsearch']
    shell:
        "vsearch --threads {threads} --cluster_fast {input} "
        "--id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
        "--centroids {output.fa1} -uc {output.uc1}\n"
        "vsearch --threads {threads} --cluster_fast {output.fa1} "
        "--iddef 0 --id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
        "--relabel LSU_ --relabel_keep "
        "--centroids {output.fa2} -uc {output.uc2}"

rule VSEARCH_2_merge:
    input:
        S = "{sample}/VSEARCH/barrnap.cdhitS.fasta",
        L = "{sample}/VSEARCH/barrnap.cdhitL.fasta"
    output:
        "{sample}/VSEARCH/barrnap.LFRs.fasta"
    shell:
        "cat {input.S} {input.L} > {output}"

rule CloseRef_1_makeIndex:
    input: "{sample}/VSEARCH/barrnap.LFRs.fasta"
    output:"{sample}/VSEARCH/contig.LFRs.fasta.index.sa"
    params:"{sample}/VSEARCH/contig.LFRs.fasta.index"
    shell: "bwa index {input} -p {params}"

rule CloseRef_2_mapping:
    input:
        R1 = "{sample}/clean/fastp.sort.1.fq",
        R2 = "{sample}/clean/fastp.sort.2.fq",
        db = "{sample}/VSEARCH/contig.LFRs.fasta.index.sa"
    output: "{sample}/VSEARCH/contig.LFRs.bwa.bam"
    params: "{sample}/VSEARCH/contig.LFRs.fasta.index"
    threads: config['thread']['blastn']
    shell:
        "bwa mem -t {threads} {params} {input.R1} {input.R2} "
        "| samtools view -b -@ {threads} > {output}"

rule Quatification:
    input: "{sample}/VSEARCH/contig.LFRs.bwa.bam"
    output:
        stat = "{sample}/VSEARCH/contig.LFRs.bwa.stat",
        sbb1 = "{sample}/VSEARCH/contig.LFRs.bwa.sbb1",
        prop = "{sample}/VSEARCH/contig.LFRs.bwa.prop"
    shell:
        "samtools view {input} | metabbq beadStat sam -i - -o {output.stat} -v\n"
        "metabbq beadStat sam2b -i {output.stat} -o {output.sbb1} -v\n"
        "awk 'FNR>1{{print $5}}' {output.sbb1}|sort|uniq -c > {output.prop}\n"
