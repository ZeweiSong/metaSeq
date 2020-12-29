#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beadbarcode detection and preparation
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

if config["method"]["QC"]["cutOA"]:
    rule BB_1_stLFR:
        input:
            r1 = "{sample}/input/rawSeq_1.fq.gz",
            r2 = "{sample}/input/rawSeq_2.fq.gz",
            bfile = config["BB_LIST"]
        params:
            title = "{sample}",
            process=config["fastp_process"],
            pos1  = config["stLFR_pos1"],
            pos2  = config["stLFR_pos2"],
            pos3  = config["stLFR_pos3"]
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
            "fastp -L --stLFR_barcode_file {input.bfile} --reads_to_process {params.process} "
            "--stLFR_pos1 {params.pos1} --stLFR_pos2 {params.pos2} --stLFR_pos3 {params.pos3} "
            "--in1 {input.r1} --in2 {input.r2} "
            "--out1 {output.r1} --out2 {output.r2} --json {output.json} --html {output.html} "
            "--disable_trim_poly_g --report_title {params.title} "
            "-w {threads} -V -B -W 30 &> {log}\n"
else:
    rule BB_1_stLFR:
        input:
            r1 = "{sample}/input/rawSeq_1.fq.gz",
            r2 = "{sample}/input/rawSeq_2.fq.gz",
            bfile = config["BB_LIST"]
        params:
            title = "{sample}",
            process=config["fastp_process"],
            pos1  = config["stLFR_pos1"],
            pos2  = config["stLFR_pos2"],
            pos3  = config["stLFR_pos3"]
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
            "fastp -L --stLFR_barcode_file {input.bfile} --reads_to_process {params.process} "
            "--stLFR_pos1 {params.pos1} --stLFR_pos2 {params.pos2} --stLFR_pos3 {params.pos3} "
            "--in1 {input.r1} --in2 {input.r2} "
            "--out1 {output.r1} --out2 {output.r2} --json {output.json} --html {output.html} "
            "--disable_trim_poly_g --report_title {params.title} "
            "-w {threads} -V &> {log}\n"

detectedFq1="{sample}/clean/fastp.1.fq"
detectedFq2="{sample}/clean/fastp.2.fq"

rule BB_2_sortR1:
    input:
        r1 = detectedFq1
    output:
        s1 = "{sample}/clean/fastp.sort.1.fq"
    params:
        tmp = "{sample}/.tmp.1"
    threads: config["thread"]["pigz"]
    shell:
        "mkdir -p {params.tmp}\n"
        "cat {input.r1} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" > {output.s1}"

rule BB_2_sortR2:
    input:
        r2 = detectedFq2
    output:
        s2 = "{sample}/clean/fastp.sort.2.fq"
    params:
        tmp = "{sample}/.tmp.2"
    threads: config["thread"]["pigz"]
    shell:
        "mkdir -p {params.tmp}\n"
        "cat {input.r2} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" > {output.s2}"

rule BB_3_idx:
    input:
        f1 = "{sample}/clean/fastp.sort.1.fq",
        f2 = "{sample}/clean/fastp.sort.2.fq"
    output:
        idx = "{sample}/clean/fastp.sort.1.fq.idx",
        stat= "{sample}/clean/BB.stat"
    shell:
        "metabbq beadsWrite3.pl -x --r1 {input.f1} -v \n"
        "awk '$1!~/0000/{{print ($3-$2+1)/4}}' {output.idx}|sort|uniq -c|sort -k2,2nr|awk '{{b=b+$1;print $0\"\\t\"b}}' > {output.stat}"

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

rule BB_5_gzip:
    input:
        s1 = "{sample}/clean/fastp.sort.1.fq",
        s2 = "{sample}/clean/fastp.sort.2.fq"
    output:
        x1 = "{sample}/clean/fastp.sort.1.fq.gz",
        x2 = "{sample}/clean/fastp.sort.2.fq.gz"
    shell:
        "gzip {input.s1} {input.s2}"

rule BB_6_removeAd:
    input:
        s1 = "{sample}/clean/fastp.sort.1.fq.gz",
        s2 = "{sample}/clean/fastp.sort.2.fq.gz"
    output:
        x1 = "{sample}/clean/fastp.sort.-ad.1.fq.gz",
        x2 = "{sample}/clean/fastp.sort.-ad.2.fq.gz"
    shell:
        "metabbq delAd.pl RCAad.seq {input.s1} {input.s2} {output.s1} {output.s2}"
