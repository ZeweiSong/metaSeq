#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Circos for BC assembles
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

if config["method"]["cluster"]["circos"]:
    rule CC_circos:
        input:
            log = "{sample}/circos_mashBC.log",
            bc = "{sample}/mash/bxxxx.bc.tree.target.cluster.count",
        output: "{sample}/batch.assemble.BC.sh"
        params:
            samDir = "{sample}",
            outDir = "{sample}/Assemble_mashBC",
            threads = config["threads"]
        shell:
            "for i in `sort -k2,2nr {input.bc} | cut -f1`; do "
            "echo metabbq bcPost.template.sh F {params.samDir} mashBC $i BC;"
            "done > {params.samDir}/batch.assemble.BC.sh\n"



# metabbq reAlign.sh {sample}/Assemble_Lv1/BC00000 REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta
# metabbq reAlign.sh {sample}/Assemble_Lv1/BC00000 REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta
# metabbq circos.sh . ../../../REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta.ids 1 97
