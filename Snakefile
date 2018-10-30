configfile: "config.yaml"

thead4fastp = 1
if config["threads"] > 6:
    thead4fastp = 6

rule all:
    input:
        expand("beadPool/{sample}.{sfx}.filter.dist",
        sample=config["samples"],sfx=["1","2"])

rule s1_stLFR:
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
        html = "clean/{sample}.fp.html",
        freq = "clean/{sample}.code.freq.tsv"
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
        "-w {threads} &> {log};\n"
        "perl src/jsonBarcode2txt.pl < {output.json} | sort -nrk2 > {output.freq}"

rule s2_sortR1:
    input:
        r1 = "clean/{sample}.fp.1.fq.gz"
    output:
        s1 = "clean/{sample}.fps.1.fq.gz"
    params:
        tmp = config["tmp"]
    shell:
        "gunzip -c {input.r1} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" |gzip > {output.s1}"

rule s2_sortR2:
    input:
        r2 = "clean/{sample}.fp.2.fq.gz"
    output:
        s2 = "clean/{sample}.fps.2.fq.gz"
    params:
        tmp = config["tmp"]
    shell:
        "gunzip -c {input.r2} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" |gzip > {output.s2}"

rule s3_saveBeadR1:
    input:
        fq = "clean/{sample}.fps.1.fq.gz",
        freq = "clean/{sample}.code.freq.tsv"
    params:
        dir = "beadPool/{sample}",
        cut = config["beadSelect"],
        log = "beadPool/{sample}.info"
    output: "beadPool/{sample}.1.info"
    benchmark:
        "benchmarks/{sample}.saveBeads1.benchmark.txt"
    shell:
        #"python src/stlfr_json2individualFASTA.py -i {input} -d {params.dir} -z {params.b0d};"
        "perl src/beadsWrite.pl {input.fq} {input.freq} {params.cut} 1 {params.dir} > {output}"

rule s3_saveBeadR2:
    input:
        fq = "clean/{sample}.fps.2.fq.gz",
        freq = "clean/{sample}.code.freq.tsv"
    params:
        dir = "beadPool/{sample}",
        cut = config["beadSelect"],
        log = "beadPool/{sample}.info"
    output: "beadPool/{sample}.2.info"
    benchmark:
        "benchmarks/{sample}.saveBeads2.benchmark.txt"
    shell:
        #"python src/stlfr_json2individualFASTA.py -i {input} -d {params.dir} -z {params.b0d};"
        "perl src/beadsWrite.pl {input.fq} {input.freq} {params.cut} 2 {params.dir} > {output}"

rule s4_makeMashR1:
    input: "beadPool/{sample}.1.info"
    params:
        dir = "beadPool/{sample}/",
        threshold = config["distCutoff"]
    output:
        msh = "beadPool/{sample}.1.msh",
        dist = "beadPool/{sample}.1.dist",
        filter = "beadPool/{sample}.1.filter.dist"
    benchmark:
        "benchmarks/{sample}.1.mash.benchmark.txt"
    shell:
        "mash sketch {params.dir}*/*.1.fa -o {output.msh}\n"
        "mash dist {output.msh} {output.msh} > {output.dist}\n"
        "perl src/filterMASHoutput.pl {output.dist} {params.threshold} {output.filter}\n"
        "sed 's/\t/,/g' {output.filter} > {output.filter}.csv"

rule s4_makeMashDistR2:
    input: "beadPool/{sample}.2.info"
    params:
        dir = "beadPool/{sample}/",
        threshold = config["distCutoff"]
    output:
        msh = "beadPool/{sample}.2.msh",
        dist = "beadPool/{sample}.2.dist",
        filter = "beadPool/{sample}.2.filter.dist"
    benchmark:
        "benchmarks/{sample}.2.mash.benchmark.txt"
    shell:
        "mash sketch {params.dir}*/*.2.fa -o {output.msh}\n"
        "mash dist {output.msh} {output.msh} > {output.dist}\n"
        "perl src/filterMASHoutput.pl {output.dist} {params.threshold} {output.filter}\n"
        "sed 's/\t/,/g' {output.filter} > {output.filter}.csv"

rule s4_makeMashDistRB:
    input:
        R1 = "beadPool/{sample}.1.msh",
        R2 = "beadPool/{sample}.2.msh"
    output:
        dist = "beadPool/{sample}.B.dist",
        filter = "beadPool/{sample}.B.filter.dist"
    params:
        threshold = config["distCutoff"]
    shell:
        "mash dist {input.R1} {input.R2} > {output.dist}\n"
        "perl src/filterMASHoutput.pl {output.dist} {params.threshold} {output.filter}\n"
