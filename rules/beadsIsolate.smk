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

rule VSEARCH_1_Dereplicate:
    input: outXc
    output:
        fa="{sample}/VSEARCH/contig.derep.full.fasta",
        uc="{sample}/VSEARCH/contig.derep.full.uc"
    threads: config['thread']['vsearch']
    shell:
        "vsearch --threads {threads} --derep_fulllength {input} "
        "--output {output.fa} -uc {output.uc} --fasta_width 0"

rule VSEARCH_2_Precluster:
    input: "{sample}/VSEARCH/contig.derep.full.fasta"
    output:
        fa="{sample}/VSEARCH/contig.preclust.fasta",
        uc="{sample}/VSEARCH/contig.preclust.uc"
    params:
        pct = config['p_VS_preClust']
    threads: config['thread']['vsearch']
    shell:
        "vsearch --threads {threads} --cluster_fast {input} "
        "--id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
        "--centroids {output.fa} -uc {output.uc}"

rule VSEARCH_3_chimeraDetection:
    input: "{sample}/VSEARCH/contig.preclust.fasta"
    output:"{sample}/VSEARCH/contig.nonchimeras.fasta"
    threads: config['thread']['vsearch']
    shell:
        "vsearch --threads {threads} --uchime_denovo {input} "
        "--fasta_width 0 --nonchimeras {output}"

rule VSEARCH_4_postCluster:
    input: "{sample}/VSEARCH/contig.nonchimeras.fasta"
    output:
        fa="{sample}/VSEARCH/contig.LFRs.fasta",
        uc="{sample}/VSEARCH/contig.LFRs.uc"
    params:
        pct = config['p_VS_postClust']
    threads: config['thread']['vsearch']
    shell:
        "vsearch --threads {threads} --cluster_fast {input} "
        "--id {params.pct} --strand both --fasta_width 0 --minuniquesize 1 "
        "--relabel LFR_ --relabel_keep "
        "--centroids {output.fa} -uc {output.uc}"

rule CloseRef_1_makeIndex:
    input: "{sample}/VSEARCH/contig.LFRs.fasta"
    output:"{sample}/VSEARCH/contig.LFRs.fasta.bwt"
    shell: "bwa index {input}"

rule CloseRef_2_mapping:
    input:
        R1 = "{sample}/clean/fastp.sort.1.fq",
        R2 = "{sample}/clean/fastp.sort.2.fq",
        db = "{sample}/VSEARCH/contig.LFRs.fasta.bwt"
    output: "{sample}/VSEARCH/contig.LFRs.bwa"
    params: "{sample}/VSEARCH/contig.LFRs.fasta"
    threads: config['thread']['blastn']
    shell:
        "bwa mem -t {threads} {params} {input.R1} {input.R2} > {output}"

rule Quatification:
    input: "{sample}/VSEARCH/contig.LFRs.bwa"
    output:
        stat = "{sample}/VSEARCH/contig.LFRs.bwa.stat",
        sbb1 = "{sample}/VSEARCH/contig.LFRs.bwa.sbb1",
        prop = "{sample}/VSEARCH/contig.LFRs.bwa.prop"
    shell:
        "metabbq beadStat sam -i {input} -o {output.stat} -v\n"
        "metabbq beadStat sam2b -i {output.stat} -o {output.sbb1} -v\n"
        "awk '$6>$9&&$6/$2>0.5{{print $5}}' {output.sbb1}|sort|uniq -c > {output.prop}\n"
