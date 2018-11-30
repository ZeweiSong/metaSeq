configfile: "config.yaml"

if config["threads"] > 4:
    thead4fastp = 4
    thread4pigz = 4
    thread4SPAdes = 4
else:
    thead4fastp = config["threads"]
    thread4pigz = config["threads"]
    thread4SPAdes  = config["threads"]

divideNum = config["divideSH"]
parallel = ["{:02d}".format(item) for item in range(divideNum)]

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

rule sum_bead:
    input:
        s1 = "clean/{sample}.fps.1.fq.gz",
        s2 = "clean/{sample}.fps.2.fq.gz",
        cd = "clean/{sample}.code.freq.tsv"
    params:
        dir = "beadPool/{sample}/sumAb10"
    output:
        i1 = "beadPool/{sample}/sumAb10/bead.1.info",
        i2 = "beadPool/{sample}/sumAb10/bead.2.info"
    shell:
        "mkdir -p {params.dir}\n"
        "perl src/beadsWrite3.pl --r1 {input.s1} -L {input.cd} -c 10 -t 8 -p bead.1 -o {params.dir} -f fq -e > {output.i1} &\n"
        "perl src/beadsWrite3.pl --r1 {input.s2} -L {input.cd} -c 10 -t 8 -p bead.2 -o {params.dir} -f fq -e > {output.i2} &\n"
        "wait\n"

rule sum_vsearch:
    input:  "beadPool/{sample}/sumAb10/bead.2.info"
    output: "VSEARCH/{sample}/sumAb10/bead.merge.derep.nonchimeras.otus.pct6.tax.stat"
    params:
        sam= "{sample}",
        threads = thead4fastp,
        REF = config["REF_ITS"]
    shell:
        "sh src/template.vsearch.sum.sh {params.threads} {params.REF} {params.sam} 10 sumAb\n"

rule s3_saveBeadR1:
    input:
        fq = "clean/{sample}.fps.1.fq.gz",
        freq = "clean/{sample}.code.freq.tsv"
    params:
        dir = "beadPool/{sample}",
        cut = config["beadSelect"],
        log = "beadPool/{sample}.s.info"
    threads: thread4pigz
    output: "beadPool/{sample}.1.info"
    benchmark:
        "benchmarks/{sample}.saveBeads1.benchmark.txt"
    shell:
        "perl src/beadsWrite2.pl --r1 {input.fq} -L {input.freq} -c {params.cut} -t {threads} -s 1 -f fafq -o {params.dir} > {output} && "
        "sort -nrk2,2 {output} > {params.log}"

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

rule s4_makeMashDistR1:
    input: "beadPool/{sample}.1.info"
    params:
        dir = "beadPool/{sample}/",
        threshold = config["distCutoff"]
    output:
        msh = "dist/{sample}/1.msh",
        dist = "dist/{sample}/1.dist"
    benchmark:
        "benchmarks/{sample}.1.mash.benchmark.txt"
    shell:
        "mash sketch {params.dir}*/*.1.fa -o {output.msh}\n"
        "mash dist {output.msh} {output.msh} > {output.dist}\n"

rule s4_makeMashDistR2:
    input: "beadPool/{sample}.2.info"
    params:
        dir = "beadPool/{sample}/",
        threshold = config["distCutoff"]
    output:
        msh = "dist/{sample}/2.msh",
        dist = "dist/{sample}/2.dist"
    benchmark:
        "benchmarks/{sample}.2.mash.benchmark.txt"
    shell:
        "mash sketch {params.dir}*/*.2.fa -o {output.msh}\n"
        "mash dist {output.msh} {output.msh} > {output.dist}\n"

rule s4_makeMashDistRB:
    input:
        R1 = "dist/{sample}/1.msh",
        R2 = "dist/{sample}/2.msh"
    output:
        dist = "dist/{sample}/B.dist"
    params:
        threshold = config["distCutoff"]
    shell:
        "mash dist {input.R1} {input.R2} > {output.dist}\n"

rule s4_filterMashDist:
    input:
        d1 = "dist/{sample}/1.dist",
        d2 = "dist/{sample}/2.dist",
        dB = "dist/{sample}/B.dist"
    output:
        f1 = "dist/{sample}/1.filter.dist",
        f2 = "dist/{sample}/2.filter.dist",
        fB = "dist/{sample}/B.filter.dist"
    params:
        threshold = config["distCutoff"]
    shell:
        "perl src/filterMASHoutput.pl {input.d1} {params.threshold} {output.f1}\n"
        "perl src/filterMASHoutput.pl {input.d1} {params.threshold} {output.f1}\n"
        "perl src/filterMASHoutput.pl {output.dB} {params.threshold} {output.fB}\n"


rule topCommunity:
    input:
        d1 = "dist/{sample}/1.dist",
        d2 = "dist/{sample}/2.dist"
    output:
        n1 = "dist/{sample}/1.node.rNum.lst",
        n2 = "dist/{sample}/2.node.rNum.lst"
    params:
        oDIR = "dist/{sample}"
    shell:
        "perl src/filterMASHoutput.pl -i {input.d1} -o {params.oDIR}/1.reNum.dist -M {params.oDIR}/1.map && "
        "convert -i {params.oDIR}/1.reNum.dist -w {params.oDIR}/1.w -o {params.oDIR}/1.bin && "
        "community {params.oDIR}/1.bin -w {params.oDIR}/1.w -l -1 > {params.oDIR}/1.graph.tree && "
        "hierarchy {params.oDIR}/1.graph.tree -l 1 > {params.oDIR}/1.node.lst"
        "perl src/change.id.pl -n {params.oDIR}/1.node.lst -m {params.oDIR}/1.map -v -o {output.n1} &\n"
        "perl src/filterMASHoutput.pl -i {input.d2} -o {params.oDIR}/2.reNum.dist -M {params.oDIR}/2.map && "
        "convert -i {params.oDIR}/2.reNum.dist -w {params.oDIR}/2.w -o {params.oDIR}/2.bin && "
        "community {params.oDIR}/2.bin -w {params.oDIR}/2.w -l -1 > {params.oDIR}/2.graph.tree && "
        "hierarchy {params.oDIR}/2.graph.tree -l 1 > {params.oDIR}/2.node.lst"
        "perl src/change.id.pl -n {params.oDIR}/1.node.lst -m {params.oDIR}/1.map -v -o {output.n2} &\nwait"

rule annoCommunity:
    input:
        n1 = "dist/{sample}/1.node.rNum.lst",
        REF= "REF/ee_its_database.fasta",
        ANN= "REF/ee_its_database.ITSx.anno"
    output: "Assemble/{sample}/Comms.log"
    params:
        sam = "{sample}",
        bDIR = "beadPool/{sample}",
        aDIR = "Assemble/{sample}"
    shell:
        "cNum=`cut -f1 {input.n1}|sort|uniq|wc -l`\n"
        "cNum=`echo $cNum -1 |bc`\n"
        "for i in `seq 0 $cNum`; do\n\t"
        "mkdir -p {params.aDIR}/Comm_$i\n\t"
        "grep \"^$i\" {input.n1} > {params.aDIR}/Comm_$i/beads.lst\n\t"
        "for p in {{1,2}};do\n\t\t"
        "for j in `cut -f2 {params.aDIR}/Comm_$i/beads.lst`;do\n\t\t\t"
        "sed 's/@\(\S*\)\/\(\S*\)\//@\\2\/\\1\//' {params.bDIR}/${{j:0:4}}/$j.$p.fq\n\t\t"
        "done > {params.aDIR}/Comm_$i/beads_$p.fq\n\t"
        "done &&\n\t"
        "spades.py -t 8 -o {params.aDIR}/Comm_$i -1 {params.aDIR}/Comm_$i/beads_1.fq -2 {params.aDIR}/Comm_$i/beads_2.fq &&\n\t"
        "ITSx -i {params.aDIR}/Comm_$i/contigs.fasta -o {params.aDIR}/Comm_$i/ITSx -t F --cpu 2 &&\n\t"
        "cat {params.aDIR}/Comm_$i/ITSx.ITS{{1,2}}.fasta > {params.aDIR}/Comm_$i/ITSx.ITSb.fasta &&\n\t"
        "vsearch --threads 4 --usearch_global {params.aDIR}/Comm_$i/ITSx.ITSb.fasta --db {input.REF} --id 0.97 "
        "--maxaccepts 0 --maxrejects 0 --blast6out {params.aDIR}/Comm_$i/ITSx.ITSb.tax.blast6 "
        "--samout {params.aDIR}/Comm_$i/ITSx.ITSb.tax.sam && \n\t"
        "awk '$4>50' {params.aDIR}/Comm_$i/ITSx.ITSb.tax.blast6|"
        "perl src/refScore.pl -r {input.ANN} -o {params.aDIR}/Comm_$i/ITSx.ITSb.tax.stat &\n"
        "done && wait && echo done > {output}"



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
        "-1 {params.iDir}/\"b1\"/\"$1\".1.fq -2 {params.iDir}/\"b1\"/\"$1\".2.fq\"}}' "
        "{input.inf1} > {output.sh}\n"

rule s5_SPAdes_run:
    input: "Assemble/{sample}.sh"
    output: "Assemble/{sample}.log"
    params:
        iDir = "beadPool/{sample}",
        oDir = "Assemble/{sample}",
        threshold = config["beadSelect"],
        pfx = "Assemble/{sample}_",
        divNum = divideNum
    shell:
        "line=`wc -l {input}|awk -v n={params.divNum} '{{printf \"%0.f\",$1/n}}'`;"
        "split -d -l $line {input} {params.pfx}\n"
        "for i in {params.pfx}??;do sh $i & done && wait && echo done > {output}"

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

rule s7_vsearch_annotation:
    input: "VSEARCH/{sample}.sh"
    output: "VSEARCH/{sample}.anno.sh"
    params:
        pos = config["REF_ITS"]
    shell:
        "awk '{{split($5,tag,\"/\");print \"sh src/template.vsearch.sh {params.threads} "
        "{params.REF} \"tag[2],tag[4]\" EEITS\"}}' {input} > {output} "
