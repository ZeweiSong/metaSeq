configfile: "config.yaml"

if config["threads"] > 4:
    thread4pigz = 4
if config["threads"] > 6:
    thread4fastp = 6
if config["threads"] > 8:
    thread4SPAdes = config["threads"]
else:
    thread4fastp = config["threads"]
    thread4pigz = config["threads"]
    thread4SPAdes  = config["threads"]

divideNum = config["divideSH"]
parallel = ["{:02d}".format(item) for item in range(divideNum)]

rule all:
    input:
        expand("beadPool/{sample}.B.filter.dist", sample=config["samples"])


rule s1_stLFR:
    input:
        r1 = "rawSeq/{sample}_1.fq.gz",
        r2 = "rawSeq/{sample}_2.fq.gz",
        bfile = "database/barcode.list"
    params:
        title = "{sample}"
    output:
        r1 = "{sample}/clean/fastp.1.fq.gz",
        r2 = "{sample}/clean/fastp.2.fq.gz",
        json = "{sample}/clean/fastp.json",
        html = "{sample}/clean/fastp.html"
    log:
        "{sample}/clean/fastp.log"
    benchmark:
        "benchmarks/{sample}.fp.benchmark.txt"
    threads: thread4fastp
    shell:
        "fastp --stLFR_barcode_file {input.bfile} --stLFR_pos3 132 "
        "--in1 {input.r1} --in2 {input.r2} "
        "--adapter_sequence AGAGTTTGATCATGGCTCAG --adapter_sequence_r2 AAGGAGGTGATCCAGCCGCA "
        "--out1 {output.r1} --out2 {output.r2} --json {output.json} --html {output.html} "
        "--disable_trim_poly_g --report_title {params.title} "
        "-w {threads} -V &> {log}\n"

rule s2_sortR1:
    input:
        r1 = "{sample}/clean/fastp.1.fq.gz"
    output:
        s1 = "{sample}/clean/fastp.sort.1.fq.gz"
    params:
        tmp = config["tmp"]
    threads: thread4pigz
    shell:
        "pigz -p {threads} -dc {input.r1} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" |pigz -p {threads} > {output.s1}"

rule s3_idxR1:
    input:  "{sample}/clean/fastp.sort.1.fq.gz"
    output: "{sample}/clean/fastp.sort.1.fq.gz.idx"
    shell:
        "metabbq beadsWrite3.pl -x --r1 {input} -v "

rule s2_sortR2:
    input:
        r2 = "{sample}/clean/fastp.2.fq.gz"
    output:
        s2 = "{sample}/clean/fastp.sort.2.fq.gz"
    params:
        tmp = config["tmp"]
    threads: thread4pigz
    shell:
        "pigz -p {threads} -dc {input.r2} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" |pigz -p {threads} > {output.s2}"

rule vsearch:
    input:
        s1 = "{sample}/clean/fastp.sort.1.fq.gz",
        s2 = "{sample}/clean/fastp.sort.2.fq.gz",
        idx= "{sample}/clean/fastp.sort.1.fq.gz.idx"
    output: "{sample}/VSEARCH/read.merge.derep.2T.bc.graph.tree"
    params:
        pfx= "{sample}/VSEARCH/read",
        threads = config["threads"]
    shell:
        "metabbq template.vsearch.preCluster.sh {params.threads} {input.s1} {input.s2} {params.pfx}\n"
