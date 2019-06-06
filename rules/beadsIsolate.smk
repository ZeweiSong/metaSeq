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
    output: "{sample}/Assemble_BI/BI.log"
    params:
        outDir = "{sample}/Assemble_BI",
        threads = config["threads"]
    shell:
        "metabbq beadsWrite3.pl -d {input.bi} -f fq -t 4 -o {params.outDir} -v "
        "--r1 {input.s1}  --r2 {input.s2} > {output}\n"

rule BIA_2_initASMsh:
    input:
        log = "{sample}/Assemble_BI/BI.log",
        bi = "{sample}/Assemble_BI/ID.lst",
    output: "{sample}/batch.assemble.BI.sh"
    params:
        samDir = "{sample}",
        outDir = "{sample}/Assemble_BI",
        mode   = config["method"]["assemble"]["mode"]
    shell:
        "for i in `sort -k2,2nr {input.bi} | cut -f1`; do "
        "echo metabbq RCAasm.sh {params.mode} {params.samDir} BI $i BI ; done > {output}\n"

outXs = "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.tsv"
outXa= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.fasta"
outXn= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.anno"

rule BIA_X_summary:
    input: "{sample}/batch.assemble.BI.sh"
    output:
        stat = outXs,
        fa   = outXa,
        anno = outXn
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        outDir = "{sample}/Assemble_BI",
        mAsm   = config["p_asm_min"],
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
# 1. Run assembly one by one
dl=$[$(wc -l {input}|tr ' ' '\\n' |head -1) / {params.divide}]
split -d -l $dl {input} {params.dpfx}
for i in {params.dpfx}??;do sh $i & done && wait

# 2. pick robust assemble fasta
for i in `ls {params.outDir}`;do
  metabbq clusterHelper asmPick -l {params.mAsm} -b $i -i {params.outDir}/$i/{params.mode}/final.contigs.fa
done > {output.fa}

# 3. stat assemble results
for i in `ls {params.outDir}`;do
  awk -v bc=$i '/^>/{{print bc"\\t"$0}}' {params.outDir}/$i/{params.mode}/final.contigs.fa | sed 's/>//;s/ flag=/\\t/;s/ multi=/\\t/;s/ len=/\\t/'
done > {output.stat}

# 4. stat annotation results
for i in `ls {params.outDir}`;do
  awk -v bc=$i '{{print bc"\\t"$0}}' {params.outDir}/$i/{params.mode}/scaffolds.megahit.BLAST.tax.blast6.anno.more
done > {output.anno}
        """
