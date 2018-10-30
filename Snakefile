configfile: "config.yaml"

thead4fastp = 1
if config["threads"] > 6:
    thead4fastp = 6

rule all:
    input:
        expand("beadPool/{sample}.{sfx}.filter.dist",
        sample=config["samples"],sfx=["1","2"])

rule s1_stLFR_fastp:
    input:
        r1 = "rawSeq/{sample}_1.fq.gz",
        r2 = "rawSeq/{sample}_2.fq.gz",
        bfile = "database/barcode.list"
    params:
        title = "{sample}"
    output:
        r1 = "clean/{sample}.fp.1.fq.gz",
        r2 = "clean/{sample}.fp.2.fq.gz",
        json = "clean/{sample}.fp.json",
        html = "clean/{sample}.fp.html"
    log:
        "clean/{sample}.fp.log"
    benchmark:
        "benchmarks/{sample}.fp.benchmark.txt"
    threads: thead4fastp
    shell:
        "fastp --stLFR_barcode_file {input.bfile} "
        "--in1 {input.r1} --in2 {input.r2} "#"--disable_adapter_trimming "
        "--adapter_sequence CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGATCACCAAGGATCGCCATAGTCCATGCTAAAGGACGTCAGGAAGGGCGATCTCAGG "
        "--adapter_sequence_r2 TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGGCGACGGCCACGAAGCTAACAGCCAATCTGCGTAACAGCCAAACCTGAGATCGCCC "
        "--out1 {output.r1} --out2 {output.r2} --json {output.json} --html {output.html} "
        "--disable_trim_poly_g --report_title {params.title} "
        "-w {threads} &> {log}"

rule s1_stLFR_sort:
    input:
        json = "clean/{sample}.fp.json"
    output:
        sort = "clean/{sample}.stlfr.lst"
    shell:
        "perl src/jsonBarcode2txt.pl < {input.json} |sort -nrk2 > {output.sort}"

rule s2_saveBeads1:
    input:
        fq = "clean/{sample}.fp.1.fq.gz",
        list = "clean/{sample}.stlfr.lst"
    params:
        dir = "beadPool/{sample}",
        cut = "1000",
        log = "beadPool/{sample}.info",
        slog = "beadPool/{sample}.1.info"
    output: "beadPool/{sample}/small.1.fa"
    benchmark:
        "benchmarks/{sample}.saveBeads1.benchmark.txt"
    shell:
        #"python src/stlfr_json2individualFASTA.py -i {input} -d {params.dir} -z {params.b0d};"
        "perl src/beadsWrite.pl {input.fq} {input.list} {params.cut} 1 {params.dir} |"
        "paste - - | sort -k2,2 -t \"/\" | tr \"\\t\" \"\\n\" > {output}"

rule s2_saveBeads2:
    input:
        fq = "clean/{sample}.fp.2.fq.gz",
        list = "clean/{sample}.stlfr.lst"
    params:
        dir = "beadPool/{sample}",
        cut = "1000",
        log = "beadPool/{sample}.info",
        slog = "beadPool/{sample}.2.info"
    output: "beadPool/{sample}/small.2.fa"
    benchmark:
        "benchmarks/{sample}.saveBeads2.benchmark.txt"
    shell:
        #"python src/stlfr_json2individualFASTA.py -i {input} -d {params.dir} -z {params.b0d};"
        "perl src/beadsWrite.pl {input.fq} {input.list} {params.cut} 2 {params.dir} | "
        "paste - - | sort -k2,2 -t \"/\" | tr \"\\t\" \"\\n\" > {output}"

rule s3_makeMash1:
    input: "beadPool/{sample}.1.log"
    params:
        dir = "beadPool/{sample}/",
    output: "beadPool/{sample}.1.msh"
    benchmark:
        "benchmarks/{sample}.1.mash.benchmark.txt"
    shell:
        "mash sketch {params.dir}*/*.1.fa -o {output}"

rule s3_makeMash2:
    input: "beadPool/{sample}.2.log"
    params:
        dir = "beadPool/{sample}/",
    output: "beadPool/{sample}.2.msh"
    benchmark:
        "benchmarks/{sample}.2.mash.benchmark.txt"
    shell:
        "mash sketch {params.dir}*/*.2.fa -o {output}"

rule s4_filterMash1:
    input: "beadPool/{sample}.1.msh"
    params:
        threshold = "0.04"
    output: "beadPool/{sample}.1.filter.dist"
    benchmark:
        "benchmarks/{sample}.1.compareMash.benchmark.txt"
    shell:
        #"ls {params.dir}/*.fa|grep -v 0000|xargs -n 1 mash dist {input} > {output}"
        #"ls {params.dir}/*.fa|xargs -n 1 mash dist {input} > {output}"
        "mash dist {input} {input} | perl src/filterMASHoutput.pl - {params.threshold} {output};"
        "sed 's/\t/,/g' {output} > {output}.csv"

rule s4_filterMash2:
    input: "beadPool/{sample}.2.msh"
    params:
        threshold = "0.04"
    output: "beadPool/{sample}.2.filter.dist"
    benchmark:
        "benchmarks/{sample}.2.compareMash.benchmark.txt"
    shell:
        #"ls {params.dir}/*.fa|grep -v 0000|xargs -n 1 mash dist {input} > {output}"
        #"ls {params.dir}/*.fa|xargs -n 1 mash dist {input} > {output}"
        "mash dist {input} {input} | perl src/filterMASHoutput.pl - {params.threshold} {output};"
        "sed 's/\t/,/g' {output} > {output}.csv"
