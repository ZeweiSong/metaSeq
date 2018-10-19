configfile: "config.yaml"

rule all:
    input:
        "seqJSON/clean.json"
    threads: 24

rule splitBarcode_cOMG:
    input:
        r1 = "rawSeq/{sample}_1.fq.gz",
        r2 = "rawSeq/{sample}_2.fq.gz"
    params:
        pfx= "OA",
        Qsys = "33",
        minLen = "30",
        Scut = "20",
        Qcut = "10",
        stLFR = "2",
        bfile = "database/barcode.list",
        debug = "2"
    output:
        r1 = "clean/{sample}.OA.clean.1.fq.gz",
        r2 = "clean/{sample}.OA.clean.2.fq.gz"
    benchmark:
        "benchmarks/{sample}.OA.benchmark.txt"
    shell:
        "cOMG OAs1 {input.r1},{input.r2} "
        "clean/{wildcards.sample}.OA {params.Qsys} {params.minLen} "
        "{params.Scut} {params.Qcut} {params.stLFR} "
        "144 {params.bfile} {params.debug}"

rule splitBarcode_fastp:
    input:
        r1 = "rawSeq/{sample}_1.fq.gz",
        r2 = "rawSeq/{sample}_2.fq.gz",
        bfile = "database/barcode.list"
    params:
        pfx= "OA",
        Qsys = "33",
        minLen = "30",
        Scut = "20",
        Qcut = "10",
        stLFR = "2",
        debug = "2"
    output:
        r1 = "clean/{sample}.fp.1.fq.gz",
        r2 = "clean/{sample}.fp.2.fq.gz"
    log:
        "clean/{sample}.fp.log"
    benchmark:
        "benchmarks/{sample}.fp.benchmark.txt"
    threads: 6
    shell:
        "fastp --stLFR_barcode_file {input.bfile} "
        "--in1 {input.r1} --in2 {input.r2} --disable_adapter_trimming "
        "--out1 {output.r1} --out2 {output.r2} --disable_quality_filtering "
        "-w {threads} &> {log}"

rule splitBarcode_py:
    input:
        r1 = "rawSeq/{sample}_1.fq.gz",
        r2 = "rawSeq/{sample}_2.fq.gz",
        bfile = "database/barcode.list"
    params:
        pfx= "OA",
        Qsys = "33",
        minLen = "30",
        Scut = "20",
        Qcut = "10",
        stLFR = "2",
        debug = "2"
    output:
        r1 = "clean/{sample}.py.1.fq.gz",
        r2 = "clean/{sample}.py.2.fq.gz"
    log:
        "clean/{sample}.py.log"
    benchmark:
        "benchmarks/{sample}.py.benchmark.txt"
    shell:
        "python src/stlfr_split.py -r1 {input.r1} -r2 {input.r2} "
        "-b {input.bfile} -bl 54 -o clean/{wildcards.sample}.py -fastq &> {log}"

rule cutadapt:
    input:
        r1 = "instance/tmp/clean.clean.1.fq.gz",
        r2 = "instance/tmp/clean.clean.2.fq.gz"
    params:
        a = "CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGATCACCAAGGATCGCCATAGTCCATGCTAAAGGACGTCAGGAAGGGCGATCTCAGG",
        A = "TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGGCGACGGCCACGAAGCTAACAGCCAATCTGCGTAACAGCCAAACCTGAGATCGCCC",
        m = "60",
        e = "0.1"
    threads: 24
    output:
        r1 = "instance/tmp/clean.r1.cutPair.fq",
        r2 = "instance/tmp/clean.r2.cutPair.fq"
    shell:
        '''
        cutadapt {input.r1} {input.r2} -a {params.a}\
        -A {params.A} -o {output.r1} -p {output.r2}\
        -m {params.m} -e {params.e} -j {threads}
        '''

rule pear:
    input:
        r1 = "instance/tmp/clean.r1.cutPair.fq",
        r2 = "instance/tmp/clean.r2.cutPair.fq"
    threads: 24
    output:
        assembledFQ = "instance/tmp/clean.assembled.fastq",
        discardedFQ = "instance/tmp/clean.discarded.fastq",
        unassForFQ = "instance/tmp/clean.unassembled.forward.fastq",
        unassRevFQ = "instance/tmp/clean.unassembled.reverse.fastq"
    shell:
        '''
        pear -f {input.r1} -r {input.r2} -o {rules.all.params.pfx} -k -j {threads}
        '''

rule makeJSON:
    input:
        assembledFQ = "instance/tmp/clean.assembled.fastq",
        unassForFQ = "instance/tmp/clean.unassembled.forward.fastq",
        unassRevFQ = "instance/tmp/clean.unassembled.reverse.fastq"
    output: "seqJSON/clean.json"
    shell:
        '''
        python ../stlfr_mergepairs2json.py -a {input.assembledFQ} -f {input.unassForFQ} -r {input.unassRevFQ} -o {output}
        '''
