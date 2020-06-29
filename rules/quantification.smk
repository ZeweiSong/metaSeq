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
        fa="{sample}/CLIP/id995def2.clust.fa"
    output:
        fns = "{sample}/PROFILE/LOTU.fa",
        ann = "{sample}/PROFILE/LOTU.fa.ann"
    shell:
        "sed 's/ .*$//' {input.fa} > {output.fns} \n bwa index {output.fns} \n"

rule LOTU_Q_map:
    input:
        ann = "{sample}/PROFILE/LOTU.fa.ann",
        fa1 = "{sample}/clean/fastp.sort.1_00.fq.gz",
        fa2 = "{sample}/clean/fastp.sort.2_00.fq.gz"
    output:
        sam = "{sample}/PROFILE/LOTU.bwa.map.bam"
    params:
        idx = "{sample}/PROFILE/LOTU.fa"
    threads: 16
    shell:
        "bwa mem -t {threads} {params.idx} {input.fa1} {input.fa2} "
        "| samtools view -b -@ 4 > {output}"

rule LOTU_Q_Quatification:
    input: "{sample}/PROFILE/LOTU.bwa.map.bam"
    output:
        stat = "{sample}/PROFILE/LOTU.bwa.stat",
        sbb1 = "{sample}/PROFILE/LOTU.bwa.sbb1",
        prop = "{sample}/PROFILE/LOTU.bwa.prop"
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


rule KRAKEN_classification:
    input:
        fa1 = "{sample}/clean/fastp.sort.1.fq.gz",
        fa2 = "{sample}/clean/fastp.sort.2.fq.gz"
    output: "{sample}/KRAKEN/rDNA.out"
    threads: 8
    shell:
        "kraken2 --threads {threads} --db $LFR/Source/REF/KRAKEN2/rDNA "
        "--paired {input.fa1} {input.fa2} > {output}"

rule KRAKEN_stat:
    input:"{sample}/KRAKEN/rDNA.out"
    output:
        bead = "{sample}/KRAKEN/rDNA.bead",
        prof = "{sample}/KRAKEN/rDNA.prof"
    params:
        db = "$LFR/Source/REF/KRAKEN2/rDNA",
        pfx= "{sample}/KRAKEN/rDNA"
    threads: 1
    shell:
        "metabbq kraken2Quant.pl {params.db} {input} {params.pfx}"

rule KRAKEN_stat_neg:
    input:"{sample}/KRAKEN/rDNA.bead"
    output: "{sample}/KRAKEN/rDNA.bead.neg.stat"
    threads: 1
    shell:
        """
awk '{{if($8 >= 0){{b=0}}else{{b=$7}}printf("%d\\t%d\\n", $3,b)}}' {input}|sort -n |uniq -c > {output}
        """
