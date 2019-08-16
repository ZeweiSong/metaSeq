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


outXd= "{sample}/blast/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.fasta"
outXb= "{sample}/blast/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.b6"
rule BIA_5_primerBlast:
    input: outXa
    output: outXb
    params:
        dbPfx   = outXd,
        primer  = config["A_primer_fa_to_find"],
        threads = config["thread"]["spades"]
    shell:
        "makeblastdb -in {input} -input_type fasta -dbtype nucl "
        "-title BI.contig -parse_seqids -out {params.dbPfx}\n"
        "blastn -num_threads {params.threads} -query {params.primer} -db {params.dbPfx} "
        "-out {output} -outfmt 6 -word_size 7 -evalue 10"

outXnt= "{sample}/blast/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.nt.b6"
rule BIA_5_ntBlast:
    input: outXa
    output: outXnt
    params:
        dbPfx   = config["REF_nt"],
        threads = config["thread"]["spades"]
    shell:
        "blastn -num_threads {params.threads} -query {input} -db {params.dbPfx} "
        "-out {output} -outfmt '6 std staxid ssciname'"

outXsilva= "{sample}/blast/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.silva.b6"
rule BIA_6_silvaBlast:
    input: outXa
    output: outXsilva
    params:
        dbPfx   = config["REF_silva"],
        threads = config["thread"]["blastn"]
    shell:
        "blastn -num_threads {params.threads} -query {input} -db {params.dbPfx} "
        "-out {output} -outfmt 6"
