#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beads Isolate Assembly Solution
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

configfile: "config.yaml"

outXac= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".clip.all.fasta"

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
    input: outXac
    output:
        uc = "{sample}/CLIP/all.clust.uc"
    threads: 8
    shell:
        "vsearch --threads {threads} --cluster_fast {input} --strand both --fasta_width 0 --sizeout "
        "--relabel LOTU_  --uc {output.uc} --id .995 --iddef 2"

rule LOTU_2_getRobustOTU:
    input:
        fa = outXac,
        uc = "{sample}/CLIP/all.clust.uc"
    output:
        fa = "{sample}/CLIP/LOTU.fa",
        log= "{sample}/CLIP/LOTU.log"
    shell:
        "metabbq beadStat ruc -s -i {input.fa} -u {input.uc} -o {output.fa} > {output.log} "

rule LOTU_assignmentS:
    input: "{sample}/CLIP/LOTU.fa"
    output:
        ssu = "{sample}/ANNO/LOTU.map.SSU.m6"
    threads: 8
    params: config["OpenRef"]["silva_SSU_blastdb"]
    shell:
        "blastn -num_threads {threads} -perc_identity 95 -word_size 141 "
        "-db {params} -query {input} -outfmt '6 std qlen slen' -out {output.ssu}\n"

rule LOTU_assignment_annoS:
    input:
        ssu = "{sample}/ANNO/LOTU.map.SSU.m6",
    output:
        ssu = "{sample}/ANNO/LOTU.map.SSU.m6.more"
    params: config["OpenRef"]["silva_SSU_taxomap"]
    shell:
        "metabbq anno.pl -l 6 {params} {input.ssu} > {output.ssu}\n"

rule LOTU_assignmentL:
    input: "{sample}/CLIP/LOTU.fa"
    output:
        lsu = "{sample}/ANNO/LOTU.map.LSU.m6"
    threads: 8
    params: config["OpenRef"]["silva_LSU_blastdb"]
    shell:
        "blastn -num_threads {threads} -perc_identity 95 "
        "-db {params} -query {input} -outfmt '6 std qlen slen' -out {output.lsu}\n"

rule LOTU_assignment_annoL:
    input:
        lsu = "{sample}/ANNO/LOTU.map.LSU.m6"
    output:
        lsu = "{sample}/ANNO/LOTU.map.LSU.m6.more"
    params: config["OpenRef"]["silva_LSU_taxomap"]
    shell:
        "metabbq anno.pl -l 6 {params} {input.lsu} > {output.lsu}\n"

rule LOTU_assignmentUNITE:
    input: "{sample}/CLIP/LOTU.fa"
    output:
        ssu = "{sample}/ANNO/LOTU.map.UNITE.m6"
    threads: 8
    params: config["OpenRef"]["unite_ITS_blastdb"]
    shell:
        "blastn -num_threads {threads} -perc_identity 95 -word_size 77 "
#        "-db $LFR/Source/REF/UNITE/sh_general_release_all_04.02.2020/sh_general_release_dynamic_all "
        "-db {params} -query {input} -outfmt '6 std qlen slen' -out {output.ssu}\n"

rule LOTU_assignment_annoUNITE:
    input:
        its = "{sample}/ANNO/LOTU.map.UNITE.m6",
    output:
        its = "{sample}/ANNO/LOTU.map.UNITE.m6.more"
    params: config["OpenRef"]["unite_ITS_taxomap"]
    shell:
        "metabbq anno.pl -l 6 {params} {input.its} > {output.its}\n"


if config["sampleType"] == "F":
    rule LOTU_assignment_merge_anno:
        input:
            lsu = "{sample}/ANNO/LOTU.map.LSU.m6.more",
            ssu = "{sample}/ANNO/LOTU.map.SSU.m6.more",
            its = "{sample}/ANNO/LOTU.map.UNITE.m6.more"
        output:
            ann = "{sample}/ANNO/LOTU.map.merge.anno",
        shell:
            "metabbq beadsAnno.pl  -m f2b -i {input.lsu},{input.ssu},{input.its} -o {output.ann} -l 6 -v\n"
else:
    rule LOTU_assignment_merge_anno:
        input:
            lsu = "{sample}/ANNO/LOTU.map.LSU.m6.more",
            ssu = "{sample}/ANNO/LOTU.map.SSU.m6.more"
        output:
            ann = "{sample}/ANNO/LOTU.map.merge.anno",
        shell:
            "metabbq beadsAnno.pl  -m b2b -i {input.lsu},{input.ssu} -o {output.ann} -l 6 -v\n"


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

###############################################################################
# Directly use CLIPS for annotation (CLIPS)
###############################################################################
rule CLIP_assignmentS:
    input: outXac
    output:
        ssu = "{sample}/ANNO/CLIP.map.SSU.m6"
    threads: 8
    params: config["OpenRef"]["silva_SSU_blastdb"]
    shell:
        "blastn -num_threads {threads} -perc_identity 95 -word_size 77 "
        "-db {params} -query {input} -outfmt '6 std qlen slen' -out {output.ssu}\n"

rule CLIP_assignment_annoS:
    input:
        ssu = "{sample}/ANNO/CLIP.map.SSU.m6",
    output:
        ssu = "{sample}/ANNO/CLIP.map.SSU.m6.more"
    params: config["OpenRef"]["silva_SSU_taxomap"]
    shell:
        "metabbq anno.pl -l 6 {params} {input.ssu} > {output.ssu}\n"

rule CLIP_assignmentL:
    input: outXac
    output:
        lsu = "{sample}/ANNO/CLIP.map.LSU.m6"
    threads: 8
    params: config["OpenRef"]["silva_LSU_blastdb"]
    shell:
        "blastn -num_threads {threads} -perc_identity 95 -word_size 77 "
        "-db {params} -query {input} -outfmt '6 std qlen slen' -out {output.lsu}\n"

rule CLIP_assignment_annoL:
    input:
        lsu = "{sample}/ANNO/CLIP.map.LSU.m6"
    output:
        lsu = "{sample}/ANNO/CLIP.map.LSU.m6.more"
    params: config["OpenRef"]["silva_LSU_taxomap"]
    shell:
        "metabbq anno.pl -l 6 {params} {input.lsu} > {output.lsu}\n"

rule CLIP_assignmentUNITE:
    input: outXac
    output:
        ssu = "{sample}/ANNO/CLIP.map.UNITE.m6"
    threads: 8
    params: config["OpenRef"]["unite_ITS_blastdb"]
    shell:
        "blastn -num_threads {threads} -perc_identity 95 -word_size 77 "
#        "-db $LFR/Source/REF/UNITE/sh_general_release_all_04.02.2020/sh_general_release_dynamic_all "
        "-db {params} -query {input} -outfmt '6 std qlen slen' -out {output.ssu}\n"

rule CLIP_assignment_annoUNITE:
    input:
        its = "{sample}/ANNO/CLIP.map.UNITE.m6",
    output:
        its = "{sample}/ANNO/CLIP.map.UNITE.m6.more"
    params: config["OpenRef"]["unite_ITS_taxomap"]
    shell:
        "metabbq anno.pl -l 6 {params} {input.its} > {output.its}\n"


if config["sampleType"] == "F":
    rule CLIP_assignment_merge_anno:
        input:
            lsu = "{sample}/ANNO/CLIP.map.LSU.m6.more",
            ssu = "{sample}/ANNO/CLIP.map.SSU.m6.more",
            its = "{sample}/ANNO/CLIP.map.UNITE.m6.more"
        output:
            bead = "{sample}/ANNO/CLIP.map.merge.bead.anno",
            clip = "{sample}/ANNO/CLIP.map.merge.clip.anno"
        shell:
            "metabbq beadsAnno.pl  -m cbr -i {input.lsu},{input.ssu},{input.its} -o {output.clip} -b {output.bead} -l 6 -v\n"
else:
    rule CLIP_assignment_merge_anno:
        input:
            lsu = "{sample}/ANNO/CLIP.map.LSU.m6.more",
            ssu = "{sample}/ANNO/CLIP.map.SSU.m6.more"
        output:
            bead = "{sample}/ANNO/CLIP.map.merge.bead.anno",
            clip = "{sample}/ANNO/CLIP.map.merge.clip.anno"
        shell:
            "metabbq beadsAnno.pl  -m cbr -i {input.lsu},{input.ssu} -o {output.clip} -b {output.bead} -l 6 -v\n"


###############################################################################
# Avoid hybridize beads by clustering clips
###############################################################################
rule KRAKEN_1_ClipClust:
    input: outXac
    output:
        uc = "{sample}/CLIP/id95def4.clust.uc"
    threads: 8
    shell:
        "vsearch --threads {threads} --cluster_fast {input} --strand both "
        " --fasta_width 0 --sizeout  --uc {output.uc} --id .95 --iddef 4 \n"

rule KRAKEN_2_listClustBID:
    input:
        fa = outXac,
        uc = "{sample}/CLIP/id95def4.clust.uc",
        anno="{sample}/ANNO/CLIP.map.merge.bead.anno"
    output: "{sample}/CLIP/id95def4.clust.fa"
    threads: 1
    shell:
        "metabbq IO uch -a {input.anno} -i {input.fa} -u {input.uc} -o {output} -c 999"

rule KRAKEN_3_ClipClust995:
    input: "{sample}/CLIP/id95def4.clust.fa"
    output:
        uc = "{sample}/CLIP/id995def2.clust.uc",
        ofa = "{sample}/CLIP/id995def2.clust.fa",
        cfa = "{sample}/KRAKEN/db/library/added.fna"
    threads: 8
    shell:
        "vsearch --threads {threads} --cluster_fast {input} --strand both --iddef 2 --id .995 "
        " --fasta_width 0 --relabel LOTU_ --uc {output.uc}  --centroids {output.ofa} \n"
        "cp -r {output.ofa} {output.cfa}"



if config["sampleType"] == "F":
    rule KRAKEN_4_makeTaxonomyTree:
        input:
            uc = "{sample}/CLIP/id995def2.clust.uc",
            anno="{sample}/ANNO/CLIP.map.merge.bead.anno"
        output: "{sample}/KRAKEN/db/data/added.txt"
        params: "Eukaryota"
        shell:
            "metabbq IO maketaxon -d {params} -u {input.uc} -a {input.anno} -t $LFR/Source/REF/KRAKEN2/rDNA/data/mergeTaxon.txt -o {output}"
else:
    rule KRAKEN_4_makeTaxonomyTree:
        input:
            uc = "{sample}/CLIP/id995def2.clust.uc",
            anno="{sample}/ANNO/CLIP.map.merge.bead.anno"
        output: "{sample}/KRAKEN/db/data/added.txt"
        params: "Bacteria"
        shell:
            "metabbq IO maketaxon -u {input.uc} -a {input.anno} -t $LFR/Source/REF/KRAKEN2/rDNA/data/mergeTaxon.txt -o {output}"


rule KRAKEN_5_makedb:
    input:
        fa = "{sample}/CLIP/id995def2.clust.fa",
        ad = "{sample}/KRAKEN/db/data/added.txt"
    output: "{sample}/KRAKEN/db/hash.k2d"
    params: "{sample}/KRAKEN/db"
    threads: 4
    shell:
        """
if [[ -e {params}/data/mergeTaxon.txt ]]; then rm {params}/data/mergeTaxon.txt ;fi
cp -r $LFR/Source/REF/KRAKEN2/rDNA/data/mergeTaxon.txt {params}/data/mergeTaxon.txt
metabbq buildKraken2db.sh {params} {threads}
        """

rule KRAKEN_6_mapping:
    input:
        fa1 = "{sample}/clean/fastp.sort.1.fq.gz",
        fa2 = "{sample}/clean/fastp.sort.2.fq.gz",
        db  = "{sample}/KRAKEN/db/hash.k2d"
    output:
        raw = "{sample}/KRAKEN/reads2db.out",
        prf = "{sample}/KRAKEN/reads2db.prof"
    params: "{sample}/KRAKEN/db"
    threads: 4
    shell:
        "kraken2 --db {params} --threads {threads} --paired {input.fa1} {input.fa2} "
        " | tee {output.raw} | metabbq IO kk2prf -d {params} -i - -o {output.prf}"


rule KRAKEN_7_sum2species:
    input:"{sample}/KRAKEN/reads2db.prof"
    output: "{sample}/KRAKEN/reads2db.sp.prf"
    threads: 1
    shell:
        "metabbq IO sumprf -r species -i {input} -o {output}"


kkXpfx="{sample}/KRAKEN/reads2" + str(config["kk2_db"])
rule KRAKEN_X_custom_mapping:
    input:
        fa1 = "{sample}/clean/fastp.sort.1.fq.gz",
        fa2 = "{sample}/clean/fastp.sort.2.fq.gz",
        db  = "{sample}/KRAKEN/" + str(config["kk2_db"])
    output:
        raw = kkXpfx + ".out",
        prf = kkXpfx + ".prof"
    params: "{sample}/KRAKEN/" + str(config["kk2_db"])
    threads: 4
    shell:
        "kraken2 --db {params} --threads {threads} --paired {input.fa1} {input.fa2} "
        " | tee {output.raw} | metabbq IO kk2prf -d {params} -i - -o {output.prf}"

rule KRAKEN_X_sum2species:
    input: kkXpfx + ".prof"
    output: kkXpfx + ".sp.prf"
    threads: 1
    shell:
        "metabbq IO sumprf -r species -i {input} -o {output}"
