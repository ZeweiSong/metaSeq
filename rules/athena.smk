#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Athena process
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

rule ATHENA_1_barcodeTag:
    input:
        i1 = "{sample}/clean/fastp.sort.1.fq",
        i2 = "{sample}/clean/fastp.sort.2.fq"
    output:  "{sample}/clean/sort_atn_sp.fastq"
    shell:
        """
        perl -e 'open I1,\"< $ARGV[0]\";open I2,\"< $ARGV[1]\";
        while(<I1>){{
          @a = split("/",$_);@b = split("/",<I2>);
          $r1="$a[0]/$a[1] BC:Z:$a[1]-1\\n".<I1>.<I1>.<I1>;
          $r2="$b[0]/$b[1] BC:Z:$b[1]-1\\n".<I2>.<I2>.<I2>;
          unless($a[1] =~ /0000/){{ print \"$r1$r2\"}}
        }}' {input.i1} {input.i2} > {output}
        """

rule ATHENA_2_assemble:
    input:  "{sample}/clean/sort_atn_sp.fastq"
    output: "{sample}/spadesMeta/contigs.fasta"
    params: "{sample}/spadesMeta"
    threads: config['thread']['athena']
    shell:
        "spades.py -t {threads} --meta -o {params} --12 {input}"

rule ATHENA_3_index:
    input:
        idx = "{sample}/spadesMeta/contigs.fasta",
        inf = "{sample}/clean/sort_atn_sp.fastq"
    output: "{sample}/athena/reads.2.metaspades-contigs.bam"
    params: "{sample}/spadesMeta"
    threads: config['thread']['athena']
    shell:
        "bwa index {input}\n"
        "bwa mem -t {threads} -C -p {input.idx} {input.inf}| "
        "samtools sort -@ 8 -o {output} - \n"
        "samtools index {output}"

rule ATHENA_4_config:
    input:
        fqs = "{sample}/clean/sort_atn_sp.fastq",
        ctg = "{sample}/spadesMeta/contigs.fasta",
        bam = "{sample}/athena/reads.2.metaspades-contigs.bam"
    output:   "{sample}/athena/config.json"
    shell:
        """
        echo {{
          "input_fqs": "../clean/sort_atn_sp.fastq",
          "ctgfasta_path": "../spadesMeta/contigs.fasta",
          "reads_ctg_bam_path": "reads.2.metaspades-contigs.bam"
        }} > {output}
        """

rule ATHENA_5_FinalShell:
    input: "{sample}/athena/config.json"
    output:"{sample}/athena/run.athena.sh"
    threads: config['thread']['athena']
    shell:
        "source activate athena_env\n"
        "athena-meta --config {input} --threads {threads}\n" # For now, run it manully

rule ATHENA_6_annotation:
    input:
        r1  = "{sample}/clean/fastp.sort.1.fq.gz",
        r2  = "{sample}/clean/fastp.sort.2.fq.gz",
        scaf = "{sample}/athena/results/olc/flye-asm-1/scaffolds.fasta",
        refFA  = config["REF_FA"],
        refID   = config["REF_ID"]
    params:
        oDir = "{sample}/athena/reAlign",
        b6   = "{sample}/athena/reAlign/scaf2ref.blast6"
    output:"{sample}/athena/reAlign/scaf2ref.blast6.anno"
    threads: config['thread']['athena']
    shell:
        "metabbq reAlign.sh {threads} {input.r1} {input.r2} {input.scaf} {input.refFA} {params.oDir}\n"
        "metabbq anno.pl {input.refID} {params.b6} > {output}"

rule ATHENA_7_quast:
    input:
        scaf  = "{sample}/athena/results/olc/flye-asm-1/scaffolds.fasta",
        refFA = config["REF_FA"],
        read1 = "{sample}/clean/fastp.sort.1.fq",
        read2 = "{sample}/clean/fastp.sort.2.fq"
    params:
        oDir = "{sample}/athena/quast_REF_output",
    output: "{sample}/athena/quast_REF_output/report.pdf"
    threads: config['thread']['athena']
    shell:
        "quast.py -t {threads} {input.scaf} \ "
        "-R {input.refFA} \ "
        "-1 {input.read1} \ "
        "-2 {input.read2} \ "
        "-o {params.oDir} \n"

rule ATHENA_8_circos:
    input:
        anno  = "{sample}/athena/reAlign/scaf2ref.blast6.anno",
        refID = config["REF_ID"]
    params:
        zoom = 10,
        ident = 97,
        iDir = "{sample}/athena/reAlign",
        oDir = "{sample}/athena/circos"
    output: "{sample}/athena/circos/circos.i97.x10.png"
    shell:
        "metabbq circos.sh {params.zoom} {params.ident} {params.iDir} {input.refID} {params.oDir}"
