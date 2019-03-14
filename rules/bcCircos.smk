#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Circos for BC assembles
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

if config["method"]["cluster"]["circos"]:
    rule circos:
        input:
            log = "{sample}/Assemble_mashBC.log",
            bc = "{sample}/mash/bMin2.bc.tree.target.cluster.count.main",
        output: "{sample}/batch.circos.BC.sh"
        params:
            refDB  = config["REF_GEN"],
            refID  = config["REF_ID"],
            zoom   = config["p_cc_zoom"],
            cut    = config["p_cc_cut"],
            outDir = "{sample}/Assemble_mashBC",
        shell:
            """
            for i in `ls {params.outDir}`; do
              echo metabbq reAlign.sh 8 {params.outDir}/$i/sort.{{1,2}}.fq \
              {params.outDir}/$i/scaffolds.fasta {params.refDB} {params.outDir}/$i/reAlign
              echo metabbq circos.sh {params.zoom} {params.cut} {params.outDir}/$i/reAlign \
              {params.refID} {params.outDir}/$i/circos;
            done > {output}
            """


# metabbq reAlign.sh {sample}/Assemble_Lv1/BC00000 REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta
# metabbq reAlign.sh {sample}/Assemble_Lv1/BC00000 REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta
# metabbq circos.sh . ../../../REF/silva132/SSU/132_SSURef_Nr99_tax_RNA.fasta.ids 1 97
