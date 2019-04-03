#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beadbarcode detection and preparation
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================


rule BB_1_stLFR:
    input:
        r1 = "{sample}/input/rawSeq_1.fq.gz",
        r2 = "{sample}/input/rawSeq_2.fq.gz",
        bfile = config["BB_LIST"]
    params:
        title = "{sample}"
    output:
        r1 = "{sample}/clean/fastp.1.fq",
        r2 = "{sample}/clean/fastp.2.fq",
        json = "{sample}/clean/fastp.json",
        html = "{sample}/clean/fastp.html"
    log:
        "{sample}/clean/fastp.log"
    benchmark:
        "benchmarks/{sample}.fp.benchmark.txt"
    threads: config["thread"]["fastp"]
    shell:
        "fastp --stLFR_barcode_file {input.bfile} --stLFR_pos3 132 "
        "--in1 {input.r1} --in2 {input.r2} "
        "--adapter_sequence AGAGTTTGATCATGGCTCAG --adapter_sequence_r2 AAGGAGGTGATCCAGCCGCA "
        "--out1 {output.r1} --out2 {output.r2} --json {output.json} --html {output.html} "
        "--disable_trim_poly_g --report_title {params.title} "
        "-w {threads} -V &> {log}\n"

rule BB_2_sortR1:
    input:
        r1 = "{sample}/clean/fastp.1.fq"
    output:
        s1 = "{sample}/clean/fastp.sort.1.fq"
    params:
        tmp = config["tmp"]
    threads: config["thread"]["pigz"]
    shell:
        "cat {input.r1} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" > {output.s1}"

rule BB_3_idxR1:
    input:  "{sample}/clean/fastp.sort.1.fq"
    output:
        idx = "{sample}/clean/fastp.sort.1.fq.idx",
        stat= "{sample}/clean/BB.stat"
    shell:
        "metabbq beadsWrite3.pl -x --r1 {input} -v \n"
        "awk '$1!~/0000/{{print ($3-$2+1)/4}}' {output.idx}|sort|uniq -c|sort -k2,2nr|awk '{{b=b+$1;print $0\"\\t\"b}}' > {output.stat}"

rule BB_2_sortR2:
    input:
        r2 = "{sample}/clean/fastp.2.fq"
    output:
        s2 = "{sample}/clean/fastp.sort.2.fq"
    params:
        tmp = config["tmp"]
    threads: config["thread"]["pigz"]
    shell:
        "cat {input.r2} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" > {output.s2}"

rule BB_4_split:
    input:
        idx= "{sample}/clean/fastp.sort.1.fq.idx",
        s1 = "{sample}/clean/fastp.sort.1.fq",
        s2 = "{sample}/clean/fastp.sort.2.fq"
    output:
        x1 = "{sample}/clean/bxxxx.sort.1.fq",
        x2 = "{sample}/clean/bxxxx.sort.2.fq"
    shell:
        "export pointer=`grep -v 0000 {input.idx}|head -1 |cut -f2` \n"
        "awk -v p=$pointer 'FNR>=p{{print}}' {input.s1} > {output.x1} & "
        "awk -v p=$pointer 'FNR>=p{{print}}' {input.s2} > {output.x2} & "
        "wait"
