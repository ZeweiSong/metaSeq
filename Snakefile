configfile: "config.yaml"

thead4fastp = 1
if config["threads"] > 6:
    thead4fastp = 6

rule all:
    input:
        "benchmarks/stLFR_summary.rst",
        expand("beadPool/{sample}.filter.dist",sample=config["samples"])

#Skip
rule s1_stLFR_cOMG:
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
        "cOMG_dev OAs1 {input.r1},{input.r2} "
        "clean/{wildcards.sample}.OA {params.Qsys} {params.minLen} "
        "{params.Scut} {params.Qcut} {params.stLFR} "
        "144 {params.bfile} {params.debug}"

rule s1_stLFR_fastp:
    input:
        r1 = "rawSeq/{sample}_1.fq.gz",
        r2 = "rawSeq/{sample}_2.fq.gz",
        bfile = "database/barcode.list"

    output:
        r1 = "clean/{sample}.fp.1.fq.gz",
        r2 = "clean/{sample}.fp.2.fq.gz"
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
        "--out1 {output.r1} --out2 {output.r2} "
        "--disable_trim_poly_g "
        #"--disable_quality_filtering "
        "-w {threads} &> {log}"

#Skip
rule s1_stLFR_py:
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
        r1 = "clean/{sample}.py.r1_split.fq",
        r2 = "clean/{sample}.py.r2_split.fq"
    log:
        "clean/{sample}.py.log"
    benchmark:
        "benchmarks/{sample}.py.benchmark.txt"
    shell:
        "python src/stlfr_split.py -r1 {input.r1} -r2 {input.r2} "
        "-b {input.bfile} -bl 54 -o clean/{wildcards.sample}.py -fastq &> {log};"

rule s1_stat:
    output:
        txt="benchmarks/stLFR_summary.txt",
        rst=report("benchmarks/stLFR_summary.rst",caption="benchmarks/stLFR_summary.rst",category="Step1")
    input:
        perl = expand("benchmarks/{sample}.OA.benchmark.txt",sample=config["samples"]),
        cpp = expand("benchmarks/{sample}.fp.benchmark.txt",sample=config["samples"]),
        python = expand("benchmarks/{sample}.py.benchmark.txt",sample=config["samples"])
    shell:
        '''
echo -e "==========  ========  ========  ========  ========  ========  ========  ========  ========  ========
 program        s       h:m:s   max_rss   max_vms   max_uss   max_pss    io_in     io_out   mean_load
==========  ========  ========  ========  ========  ========  ========  ========  ========  ========" > {output.txt};
for s in {input.perl} {input.cpp} {input.python} ;do
    awk 'FNR==2{{print FILENAME"  "$0}}' $s;
done | sed 's/benchmarks\///;s/.benchmark.txt/  /' >> {output.txt}
echo -e "==========  ========  ========  ========  ========  ========  ========  ========  ========  ========">> {output.txt};
cat {output.txt} | column -t > {output.rst}
        '''

#Skip
rule s2_cutadapt:
    input:
        r1 = "clean/{sample}.fp.1.fq.gz",
        r2 = "clean/{sample}.fp.2.fq.gz"
    params:
        a = "CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGATCACCAAGGATCGCCATAGTCCATGCTAAAGGACGTCAGGAAGGGCGATCTCAGG",
        A = "TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGGCGACGGCCACGAAGCTAACAGCCAATCTGCGTAACAGCCAAACCTGAGATCGCCC",
        m = "60",
        e = "0.1"
    threads: config["threads"]
    output:
        r1 = "clean/{sample}.cutadapt.1.fq.gz",
        r2 = "clean/{sample}.cutadapt.2.fq.gz"
    benchmark:
        "benchmarks/{sample}.cutadapt.benchmark.txt"
    shell:
        "cutadapt {input.r1} {input.r2} -a {params.a} "
        "-A {params.A} -o {output.r1} -p {output.r2} "
        "-m {params.m} -e {params.e} -j {threads}"

#Skip
rule s3_pear:
    input:
        r1 = "clean/{sample}.fp.1.fq.gz",
        r2 = "clean/{sample}.fp.2.fq.gz"
    threads: config["threads"]
    output:
        assembledFQ = "clean/{sample}.pear.assembled.fastq",
        discardedFQ = "clean/{sample}.pear.discarded.fastq",
        unassForFQ = "clean/{sample}.pear.unassembled.forward.fastq",
        unassRevFQ = "clean/{sample}.pear.unassembled.reverse.fastq"
    benchmark:
        "benchmarks/{sample}.pear.benchmark.txt"
    shell:
        "pear -f {input.r1} -r {input.r2} -o clean/{wildcards.sample}.pear -k -j {threads}"

#Skip
rule s4_makeJSON:
    input:
        assembledFQ = "clean/{sample}.pear.assembled.fastq",
        unassForFQ = "clean/{sample}.pear.unassembled.forward.fastq",
        unassRevFQ = "clean/{sample}.pear.unassembled.reverse.fastq"
    output: "clean/{sample}.merge.json"
    benchmark:
        "benchmarks/{sample}.makeJSON.benchmark.txt"
    shell:
        "python src/stlfr_mergepairs2json.py -a {input.assembledFQ} "
        "-f {input.unassForFQ} -r {input.unassRevFQ} -o {output}"

#Skip
rule s5_QC:
    input: "clean/{sample}.merge.json"
    output: "clean/{sample}.clean.json"
    benchmark:
        "benchmarks/{sample}.QC.benchmark.txt"
    shell:
        "python src/stlfr_qc_derep.py -i {input} -maxee 1 -o {output}"

#Skip
rule s6_countBeads:
    input: "clean/{sample}.clean.json"
    output:
        count = "clean/{sample}.beads.count",
        dist = "clean/{sample}.beads.dist"
    benchmark:
        "benchmarks/{sample}.countBeads.benchmark.txt"
    shell:
        "python src/stlfr_count.py -i {input} -o {output.count} -d {output.dist}"

rule s7_saveBeads1:
    input: "clean/{sample}.fp.1.fq.gz"
    params:
        dir = "beadPool/{sample}",
        cut = "5",
        log = "beadPool/{sample}.info",
        slog = "beadPool/{sample}.1.info"
    output: "beadPool/{sample}.1.log"
    benchmark:
        "benchmarks/{sample}.saveBeads1.benchmark.txt"
    shell:
        #"python src/stlfr_json2individualFASTA.py -i {input} -d {params.dir} -z {params.b0d};"
        "perl src/beadsWrite.pl {input} {params.cut} 1 {params.dir} > {output};"
        "sort -nrk2 {params.log} > {params.slog}"

rule s7_saveBeads2:
    input: "clean/{sample}.fp.2.fq.gz"
    params:
        dir = "beadPool/{sample}",
        cut = "5",
        log = "beadPool/{sample}.info",
        slog = "beadPool/{sample}.2.info"
    output: "beadPool/{sample}.2.log"
    benchmark:
        "benchmarks/{sample}.saveBeads2.benchmark.txt"
    shell:
        #"python src/stlfr_json2individualFASTA.py -i {input} -d {params.dir} -z {params.b0d};"
        "perl src/beadsWrite.pl {input} {params.cut} 2 {params.dir} | tee {output};"
        "sort -nrk2 {params.log} > {params.slog}"

rule s8_makeMash1:
    input: "beadPool/{sample}.1.log"
    params:
        dir = "beadPool/{sample}/",
    output: "beadPool/{sample}.1.msh"
    benchmark:
        "benchmarks/{sample}.1.mash.benchmark.txt"
    shell:
        "mash sketch {params.dir}*/*.1.fa -o {output}"

rule s8_makeMash2:
    input: "beadPool/{sample}.2.log"
    params:
        dir = "beadPool/{sample}/",
    output: "beadPool/{sample}.2.msh"
    benchmark:
        "benchmarks/{sample}.2.mash.benchmark.txt"
    shell:
        "mash sketch {params.dir}*/*.2.fa -o {output}"

rule s9_filterMash1:
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

rule s9_filterMash2:
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
