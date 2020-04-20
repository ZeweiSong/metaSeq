#!/usr/bin/env snakemake
# (c) 2016 - 2020 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beads Isolate Assembly Solution
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# Create date:       20 Apr 2020
# ===================================================================

## needs results from ./annotation.smk

###############################################################################
# mapping to Long OTUs (LOTUs)
###############################################################################
rule LOTU_Q_index:
    input:
        fa="{sample}/CLIP/all.clust.fa"
    output:
        fns = "{sample}/CLIP/all.LOTU.fa",
        ann = "{sample}/CLIP/all.LOTU.fa.ann"
    shell:
        "sed 's/ .*$//' {input.fa} > {output.fns} \n bwa index {output.fns} \n"

rule LOTU_Q_map:
    input:
        ann = "{sample}/CLIP/all.LOTU.fa.ann",
        fa1 = "{sample}/clean/fastp.sort.1.fq.gz",
        fa2 = "{sample}/clean/fastp.sort.2.fq.gz"
    output:
        sam = "{sample}/CLIP/all.LOTU.bwa.map.bam"
    params:
        idx = "{sample}/CLIP/all.LOTU.fa"
    threads: 8
    shell:
        "bwa mem -t {threads} {params.idx} {input.fa1} {input.fa2} "
        "| samtools view -b -@ 4 > {output}"

rule LOTU_Q_Quatification:
    input: "{sample}/CLIP/all.LOTU.bwa.map.bam"
    output:
        stat = "{sample}/CLIP/all.LOTU.bwa.stat",
        sbb1 = "{sample}/CLIP/all.LOTU.bwa.sbb1",
        prop = "{sample}/CLIP/all.LOTU.bwa.prop"
    shell:
        "samtools view {input} | metabbq beadStat sam -i - -o {output.stat} -v\n"
        "metabbq beadStat sam2b -i {output.stat} -o {output.sbb1} -v\n"
        "awk 'FNR>1{{print $5}}' {output.sbb1}|sort|uniq -c > {output.prop}\n"

###############################################################################
# mapping to subunits OTUs (SOTUs)
###############################################################################
rule SOTU_Q_index:
    input:
        SSU="{sample}/OTU/pred.SSU.clust.robust.fa",
        LSU="{sample}/OTU/pred.LSU.clust.robust.fa"
    output:
        fa  = "{sample}/OTU/pred.bwa.index.fa",
        ann = "{sample}/OTU/pred.bwa.index.fa.ann"
    shell:
        "cat {input.SSU} {input.LSU}|sed 's/ .*$//' >{output.fa} \n bwa index {output.fa} \n"

rule SOTU_Q_map:
    input:
        ann = "{sample}/OTU/pred.bwa.index.fa.ann",
        fa1 = "{sample}/clean/fastp.sort.1.fq.gz",
        fa2 = "{sample}/clean/fastp.sort.2.fq.gz"
    output:
        sam = "{sample}/OTU/pred.bwa.map.bam"
    params:
        idx = "{sample}/OTU/pred.bwa.index.fa"
    threads: 8
    shell:
        "bwa mem -t {threads} {params.idx} {input.fa1} {input.fa2} "
        "| samtools view -b -@ 4 > {output}"

rule SOTU_Q_Quatification:
    input: "{sample}/OTU/pred.bwa.map.bam"
    output:
        stat = "{sample}/OTU/pred.RSUs.bwa.stat",
        sbb1 = "{sample}/OTU/pred.RSUs.bwa.sbb1",
        prop = "{sample}/OTU/pred.RSUs.bwa.prop"
    shell:
        "samtools view {input} | metabbq beadStat sam -i - -o {output.stat} -v\n"
        "metabbq beadStat sam2b -i {output.stat} -o {output.sbb1} -v\n"
        "awk 'FNR>1{{print $5}}' {output.sbb1}|sort|uniq -c > {output.prop}\n"
