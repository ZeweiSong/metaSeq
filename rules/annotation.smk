#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beads Isolate Assembly Solution
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

configfile: "config.yaml"

outXc= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".clip.fasta"

if config["sampleType"] == "F":
    kingdom="euk"
    LSU="28S"
    SSU="18S"
else:
    kingdom="bac"
    LSU="23S"
    SSU="16S"

###############################################################################
# cluster clips to create Long OTUs (LOTUs)
###############################################################################
rule LOTU_1_getClipOTU:
    input: outXc
    output:
        uc = "{sample}/CLIP/all.clust.uc"
    threads: 8
    shell:
        "vsearch --threads {threads} --cluster_fast {input} --strand both --fasta_width 0 --sizeout "
        "--relabel LOTU_  --uc {output.uc} --id .99 -iddef 0"

rule LOTU_2_getRobustOTU:
    input:
        fa = outXc,
        uc = "{sample}/CLIP/all.clust.uc"
    output:
        fa = "{sample}/CLIP/all.clust.fa",
        log= "{sample}/CLIP/all.clust.log"
    shell:
        "metabbq beadStat ruc -i {input.fa} -u {input.uc} -o {output.fa} > {output.log} "

rule LOTU_assignmentS:
    input: "{sample}/CLIP/all.LOTU.fa"
    output:
        ssu = "{sample}/CLIP/all.LOTU.map.SSU.m6"
    threads: 8
    shell:
        "blastn -num_threads {threads} -perc_identity 95 "
        "-db $STL/Source/REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta "
        "-query {input} -outfmt '6 std qlen slen' -out {output.ssu}\n"

rule LOTU_assignment_annoS:
    input:
        ssu = "{sample}/CLIP/all.LOTU.map.SSU.m6",
    output:
        ssu = "{sample}/CLIP/all.LOTU.map.SSU.m6.more",
    shell:
        "metabbq anno.pl -l 6 $STL/Source/REF/silva132/SSU/taxmap_embl_ssu_ref_132.tax {input.ssu} > {output.ssu}\n"

rule LOTU_assignmentL:
    input: "{sample}/CLIP/all.LOTU.fa"
    output:
        lsu = "{sample}/CLIP/all.LOTU.map.LSU.m6"
    threads: 8
    shell:
        "blastn -num_threads {threads} -perc_identity 95 "
        "-db $STL/Source/REF/silva132/LSU/SILVA_132_LSURef_tax_RNA.fasta "
        "-query {input} -outfmt '6 std qlen slen' -out {output.lsu}\n"

rule LOTU_assignment_annoL:
    input:
        lsu = "{sample}/CLIP/all.LOTU.map.LSU.m6"
    output:
        lsu = "{sample}/CLIP/all.LOTU.map.LSU.m6.more"
    shell:
        "metabbq anno.pl -l 6 $STL/Source/REF/silva132/LSU/taxmap_embl_lsu_ref_132.tax {input.lsu} > {output.lsu}\n"


###############################################################################
# cluster predicted subunits OTUs (SOTUs)
###############################################################################
outRc= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".rRNA.fasta"

rule SOTU_0_rrnadetect:
    input: outXc
    output: outRc
    log: outRc + ".barrnap"
    params:
        k = kingdom,
    threads: config['thread']['vsearch']
    shell:
        "barrnap --kingdom {params.k} --threads {threads} --reject 0.1 {input} --outseq {output} &> {log}"

rule SOTU_1_getSubunitsOTU:
    input: outRc
    output:
        SSU="{sample}/OTU/pred.SSU.fasta",
        LSU="{sample}/OTU/pred.LSU.fasta"
    params:
        SSU=SSU,
        LSU=LSU
    shell:
        "grep -A1 \"{params.SSU}_rRNA\" {input} > {output.SSU} &\n"
        "grep -A1 \"{params.LSU}_rRNA\" {input} > {output.LSU} \n"


rule SOTU_2_SSU_vsearch_cluster:
    input: "{sample}/OTU/pred.SSU.fasta"
    output:
        ffa = "{sample}/OTU/pred.SSU.1kb.fasta",
        cfa = "{sample}/OTU/pred.SSU.clust.fa",
        cuc = "{sample}/OTU/pred.SSU.clust.uc"
    params:
    threads: 8
    shell:
        "vsearch --threads {threads} --fastx_filter {input} --fastq_minlen 1000 --fastaout {output.ffa}\n"
        "vsearch --threads {threads} --cluster_fast {output.ffa} --strand both --fasta_width 0 --sizeout "
        "--relabel SSU_ --centroids {output.cfa} --uc {output.cuc} --id .99\n"

rule SOTU_2_LSU_vsearch_cluster:
    input: "{sample}/OTU/pred.LSU.fasta"
    output:
        ffa = "{sample}/OTU/pred.LSU.1kb.fasta",
        cfa = "{sample}/OTU/pred.LSU.clust.fa",
        cuc = "{sample}/OTU/pred.LSU.clust.uc"
    params:
    threads: 8
    shell:
        "vsearch --threads {threads} --fastx_filter {input} --fastq_minlen 1000 --fastaout {output.ffa}\n"
        "vsearch --threads {threads} --cluster_fast {output.ffa} --strand both --fasta_width 0 --sizeout "
        "--relabel LSU_ --centroids {output.cfa} --uc {output.cuc} --id .99\n"

rule SOTU_3_getRobustSSU:
    input:
        fa = "{sample}/OTU/pred.SSU.fasta",
        uc = "{sample}/OTU/pred.SSU.clust.uc"
    output:
        fa = "{sample}/OTU/pred.SSU.clust.robust.fa",
        log= "{sample}/OTU/pred.SSU.clust.robust.log"
    shell:
        "metabbq beadStat ruc -l SSU -i {input.fa} -u {input.uc} -o {output.fa} > {output.log} "

rule SOTU_4_assignmentS:
    input: "{sample}/OTU/pred.SSU.clust.robust.fa"
    output:
        ssu = "{sample}/OTU/pred.SSU.clust.robust.m6"
    threads: 8
    shell:
        "blastn -num_threads {threads} -perc_identity 95 "
        "-db $STL/Source/REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta "
        "-query {input} -outfmt '6 std qlen slen' -out {output.ssu}\n"

rule SOTU_5_assignment_annoS:
    input:
        ssu = "{sample}/OTU/pred.SSU.clust.robust.m6",
    output:
        ssu = "{sample}/OTU/pred.SSU.clust.robust.m6.more",
    shell:
        "metabbq anno.pl -l 6 $STL/Source/REF/silva132/SSU/taxmap_embl_ssu_ref_132.tax {input.ssu} > {output.ssu}\n"

rule SOTU_3_getRobustLSU:
    input:
        fa = "{sample}/OTU/pred.LSU.fasta",
        uc = "{sample}/OTU/pred.LSU.clust.uc"
    output:
        fa = "{sample}/OTU/pred.LSU.clust.robust.fa",
        log= "{sample}/OTU/pred.LSU.clust.robust.log"
    shell:
        "metabbq beadStat ruc -l LSU -i {input.fa} -u {input.uc} -o {output.fa} > {output.log} "

rule SOTU_4_assignmentL:
    input: "{sample}/OTU/pred.LSU.clust.robust.fa"
    output:
        LSU = "{sample}/OTU/pred.LSU.clust.robust.m6"
    threads: 8
    shell:
        "blastn -num_threads {threads} -perc_identity 95 "
        "-db $STL/Source/REF/silva132/LSU/SILVA_132_LSURef_tax_RNA.fasta "
        "-query {input} -outfmt '6 std qlen slen' -out {output.LSU}\n"

rule SOTU_5_assignment_annoL:
    input:
        LSU = "{sample}/OTU/pred.LSU.clust.robust.m6",
    output:
        LSU = "{sample}/OTU/pred.LSU.clust.robust.m6.more",
    shell:
        "metabbq anno.pl -l 6 $STL/Source/REF/silva132/LSU/taxmap_embl_lsu_ref_132.tax {input.LSU} > {output.LSU}\n"
