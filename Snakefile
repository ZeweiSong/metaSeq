configfile: "config.yaml"


if config["threads"] > 4:
    thead4fastp = 4
    thread4pigz = 4
    thread4SPAdes = 4
else:
    thead4fastp = config["threads"]
    thread4pigz = config["threads"]
    thread4SPAdes  = config["threads"]

rule all:
    input:
        expand("beadPool/{sample}.B.filter.dist",
        sample=config["samples"])

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
    threads: thread4pigz
    shell:
        "pigz -p {threads} -dc {input.r1} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" |pigz -p {threads} > {output.s1}"

rule s2_sortR2:
    input:
        r2 = "clean/{sample}.fp.2.fq.gz"
    output:
        s2 = "clean/{sample}.fps.2.fq.gz"
    params:
        tmp = config["tmp"]
    threads: thread4pigz
    shell:
        "pigz -p {threads} -dc {input.r2} | paste - - - - | sort -T {params.tmp} -k2,2 -t \"/\" | "
        "tr \"\\t\" \"\\n\" |pigz -p {threads} > {output.s2}"

rule s3_saveBeadR1:
    input:
        fq = "clean/{sample}.fps.1.fq.gz",
        freq = "clean/{sample}.code.freq.tsv"
    params:
        dir = "beadPool/{sample}",
        cut = config["beadSelect"],
        log = "beadPool/{sample}.info"
    threads: thread4pigz
    output: "beadPool/{sample}.1.info"
    benchmark:
        "benchmarks/{sample}.saveBeads1.benchmark.txt"
    shell:
        "perl src/beadsWrite2.pl --r1 {input.fq} -L {input.freq} -c {params.cut} -p {threads} -s 1 -f fafq -o {params.dir} > {output}"

rule s3_saveBeadR2:
    input:
        fq = "clean/{sample}.fps.2.fq.gz",
        freq = "clean/{sample}.code.freq.tsv"
    params:
        dir = "beadPool/{sample}",
        cut = config["beadSelect"],
        log = "beadPool/{sample}.info"
    threads: thread4pigz
    output: "beadPool/{sample}.2.info"
    benchmark:
        "benchmarks/{sample}.saveBeads2.benchmark.txt"
    shell:
        "perl src/beadsWrite2.pl --r1 {input.fq} -L {input.freq} -c {params.cut} -p {threads} -s 2 -f fafq -o {params.dir} > {output}"

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

rule s5_SPAdes_shell:
    input:
        inf1 = "beadPool/{sample}.1.info",
        inf2 = "beadPool/{sample}.2.info",
    output:
        sh = "Assemble/{sample}.sh"
    params:
        iDir = "beadPool/{sample}",
        oDir = "Assemble/{sample}",
        threshold = config["beadSelect"],
        threads = thread4SPAdes
    shell:
        "awk '($2>={params.threshold}){{b1=substr($1,1,4);print \""
        "spades.py -t {params.threads} -o {params.oDir}/\"b1\"/\"$1\" "
        "-1 {params.iDir}/\"b1\"/\"$1\".1.fq.gz -2 {params.iDir}/\"b1\"/\"$1\".2.fq.gz\"}}' "
        "{input.inf1} > {output.sh}"

rule s5_bwa_each_bead:
    input: "Assemble/{sample}.sh"
    output: "Assemble/bwa.{sample}.sh"
    params:
        threads = thread4SPAdes,
        REF = config["REF_GEN"]
    shell:
        "awk '{{print \"bwa mem -t {params.threads} {params.REF} \"$7,$9\" -o \"$5\"/bwa.reads.sam\";"
        "print \"bwa bwasw -t {params.threads} {params.REF} \"$5\"/scaffolds.fasta > \"$5\"/bwa.scaffolds.sam\"}}' "
        " {input} > {output} "

rule s6_vsearch_EE_each_bead:
    input: "Assemble/{sample}.sh"
    output: "VSEARCH/{sample}.sh"
    params:
        threads = thread4SPAdes,
        REF = config["REF_ITS"]
    shell:
        "awk '{{split($5,tag,\"/\");print \"sh src/template.vsearch.sh {params.threads} "
        "{params.REF} \"tag[2],tag[4]\" EEITS\"}}' {input} > {output} "
