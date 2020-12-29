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
        "blastn -num_threads {threads} " #-perc_identity 95"
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
        "blastn -num_threads {threads} " # -perc_identity 95 "
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
        "blastn -num_threads {threads} " #-perc_identity 95 -word_size 77 "
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
    threads: 16
    params: config["OpenRef"]["silva_SSU_blastdb"]
    shell:
        "blastn -num_threads {threads} -word_size 77 " #"-perc_identity 95 "
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
    threads: 16
    params: config["OpenRef"]["silva_LSU_blastdb"]
    shell:
        "blastn -num_threads {threads} -word_size 77 " #"-perc_identity 95 "
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
    threads: 16
    params: config["OpenRef"]["unite_ITS_blastdb"]
    shell:
        "blastn -num_threads {threads} -word_size 77 " #"-perc_identity 95"
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
        uc = "{sample}/CLIP/id90def4.clust.uc"
    threads: 16
    shell:
        "vsearch --threads {threads} --cluster_fast {input} --strand both "
        " --fasta_width 0 --sizeout  --uc {output.uc} --id .90 --iddef 4 \n"

rule KRAKEN_2_listClustBID:
    input:
        fa = outXac,
        uc = "{sample}/CLIP/id90def4.clust.uc",
        anno="{sample}/ANNO/CLIP.map.merge.bead.anno",
    output: "{sample}/CLIP/id90def4.clust.fa"
    threads: 1
    shell:
        "metabbq IO uch -a {input.anno} -i {input.fa} -u {input.uc} -o {output} -c 800 -v"

if config["sampleType"] == "F":
    rule KRAKEN_2_validseq:
        input:
            fa = outXac,
            gff = "{sample}/CLIP/all.barrnap.gff",
            txt = "{sample}/CLIP/all.ITS.positions.txt",
            anno="{sample}/ANNO/CLIP.map.merge.bead.anno",
        output: "{sample}/CLIP/validSeqs.clust.fa"
        params:
            min = config["p_pba_len_min"],
            max = config["p_pba_len_max"]
        threads: 1
        shell:
            "metabbq IO validseq -a {input.anno} -i {input.fa} -g {input.gff} -t {input.txt} -d ALL -m {params.min} -M {params.max} -o {output} -v"

    rule KRAKEN_2_hmmseq:
        input:
            fa = outXac,
            gff = "{sample}/CLIP/all.barrnap.gff",
            txt = "{sample}/CLIP/all.ITS.positions.txt"
        output: "{sample}/CLIP/hmmSeqs.clust.fa"
        params:
            min = config["p_pba_len_min"],
            max = config["p_pba_len_max"]
        threads: 1
        shell:
            "metabbq IO validseq -i {input.fa} -g {input.gff} -t {input.txt} -d ALL -m {params.min} -M {params.max} -o {output} -v"
else:
    rule KRAKEN_2_validseq:
        input:
            fa = outXac,
            gff = "{sample}/CLIP/all.barrnap.gff",
            anno="{sample}/ANNO/CLIP.map.merge.bead.anno",
        output: "{sample}/CLIP/validSeqs.clust.fa"
        params:
            min = config["p_pba_len_min"],
            max = config["p_pba_len_max"]
        threads: 1
        shell:
            "metabbq IO validseq -a {input.anno} -i {input.fa} -g {input.gff} -d SSU -m {params.min} -M {params.max}  -o {output} -v"

    # rule KRAKEN_2_validseq2:
    #     input:
    #         fa = "{sample}/CLIP/all.barrnap.fa",
    #         gff = "{sample}/CLIP/all.barrnap.gff",
    #         anno="{sample}/ANNO/CLIP.map.merge.bead.anno",
    #     output: "{sample}/CLIP/validSeqs.clust.fa"
    #     params:
    #         min = config["p_pba_len_min"],
    #         max = config["p_pba_len_max"]
    #     threads: 1
    #     shell:
    #         "metabbq IO validseq2 -a {input.anno}  -i {input.fa}  -t barrnap -d SSU -m {params.min} -M {params.max}  -o {output} -v"

    rule KRAKEN_2_hmmseq:
        input:
            fa = outXac,
            gff = "{sample}/CLIP/all.barrnap.gff"
        output: "{sample}/CLIP/hmmSeqs.clust.fa"
        params:
            min = config["p_pba_len_min"],
            max = config["p_pba_len_max"]
        threads: 1
        shell:
            "metabbq IO validseq -i {input.fa} -g {input.gff} -d SSU -m {params.min} -M {params.max} -o {output} -v"

# rule KRAKEN_3_ClipClust995:
#     input: "{sample}/CLIP/id95def4.clust.fa"
#     output:
#         uc = "{sample}/CLIP/id995def2.clust.uc",
#         ofa = "{sample}/CLIP/id995def2.clust.fa",
#         cfa = "{sample}/KRAKEN/db/library/added.fna"
#     threads: 8
#     shell:
#         "vsearch --threads {threads} --cluster_fast {input} --strand both --iddef 2 --id .995 "
#         " --fasta_width 0 --relabel LOTU_ --uc {output.uc}  --centroids {output.ofa} \n"
#         "cp -r {output.ofa} {output.cfa}"
rule KRAKEN_3_ClipClust995:
    input:  "{sample}/CLIP/id90def4.clust.fa",
    output: "{sample}/CLIP/id90def4.clade.uc_0.995"
    params: "{sample}/CLIP/id90def4.clade.uc"
    threads: 16
    shell:
        "metabbq IO makeclade -i {input} -o {params} -t {threads} -v"


if config["sampleType"] == "F":
    rule KRAKEN_4_makeTaxonomyTree:
        input:
            fa = "{sample}/CLIP/id90def4.clust.fa",
            uc = "{sample}/CLIP/id90def4.clade.uc_0.995",
            anno="{sample}/ANNO/CLIP.map.merge.bead.anno"
        output:
            path = "{sample}/KRAKEN/db/data/added.txt",
            fna  = "{sample}/KRAKEN/db/library/added.fna"
        params:
            db = "{sample}/KRAKEN/db",
            mode = "euk",
            kk2db  = "$LFR/Source/REF/KRAKEN2/rDNA",
            ucpfx  = "{sample}/CLIP/id90def4.clade.uc"
        shell:
            "metabbq IO clade2tree -i {input.fa} -a {input.anno} -d {params.kk2db} -s {params.ucpfx} -m {params.mode} -o {params.db} -v"
else:
    rule KRAKEN_4_makeTaxonomyTree:
        input:
            fa = "{sample}/CLIP/id90def4.clust.fa",
            uc = "{sample}/CLIP/id90def4.clade.uc_0.995",
            anno="{sample}/ANNO/CLIP.map.merge.clip.anno"
        output:
            path = "{sample}/KRAKEN/db/data/added.txt",
            fna  = "{sample}/KRAKEN/db/library/added.fna"
        params:
            db = "{sample}/KRAKEN/db",
            mode = "bac",
            kk2db  = "$LFR/Source/REF/KRAKEN2/rDNA",
            ucpfx  = "{sample}/CLIP/id90def4.clade.uc"
        shell:
            "metabbq IO clade2tree -i {input.fa} -a {input.anno} -d {params.kk2db} -m {params.mode} -s {params.ucpfx} -o {params.db} -v"

rule KRAKEN_5_makedb:
    input:
        path = "{sample}/KRAKEN/db/data/added.txt",
        fna = "{sample}/KRAKEN/db/library/added.fna"
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

kkAdpfx="{sample}/KRAKEN/reads2" + str(config["kk2_db"])
rule KRAKEN_Ad_custom_mapping:
    input:
        fa1 = "{sample}/clean/fastp.sort.1.fq.gz",
        fa2 = "{sample}/clean/fastp.sort.2.fq.gz",
        db  = "{sample}/KRAKEN/" + str(config["kk2_db"])
    output:
        rpt = kkAdpfx + ".kreport2",
        out = kkAdpfx + ".kraken2",
        prf = kkAdpfx + ".k2.prof",
        bed = kkAdpfx + ".k2.bead"
    threads: 4
    shell:
        "kraken2 --db {input.db} --threads {threads} --paired {input.fa1} {input.fa2} "
        "--report {output.rpt} | tee {output.out} | metabbq IO kk2prf -d {input.db} -i - -o {output.prf}"
#        " | tee {output.raw} | metabbq IO kk2prf -d {params} -i - -o {output.prf}"
#"kraken2 --db=SAM/SBR/KRAKEN/db --threads 4 --report ${SAMPLE}/KRAKEN/reads2SBR.kreport2 \
#--paired ${SAMPLE}/clean/fastp.sort.1.fq.gz ${SAMPLE}/clean/fastp.sort.2.fq.gz > ${SAMPLE}/KRAKEN/reads2SBR.kraken2"
rule KRAKEN_Ad_custom_kreport:
    input:
        bed = kkAdpfx + ".k2.bead",
        db  = "{sample}/KRAKEN/" + str(config["kk2_db"])
    output:
        rpt = kkAdpfx + ".bead.kreport2"
    threads: 4
    shell:
        "awk '$5>0.9&&$8<10' {input.bed} | metabbq kraken-report --d={input.db} - -v > {output.rpt}"

rule KRAKEN_Ad_sum2species:
    input:
        report = kkAdpfx + ".bead.kreport2",
        db  = "{sample}/KRAKEN/" + str(config["kk2_db"])
    output:
        profile = kkAdpfx + ".bead.bracken",
        report = kkAdpfx + ".bead.breport2"
    threads: 1
    shell:
        "bracken -d {input.db} -i {input.report} -o {output.profile} -l S > {output.report}"


kkXpfx="{sample}/KRAKEN/rrad2" + str(config["kk2_db"])
rule KRAKEN_X_custom_mapping:
    input:
        fa1 = "{sample}/clean/fastp.sort.-ad.1.fq.gz",
        fa2 = "{sample}/clean/fastp.sort.-ad.2.fq.gz",
        db  = "{sample}/KRAKEN/" + str(config["kk2_db"])
    output:
        rpt = kkXpfx + ".kreport2",
        out = kkXpfx + ".kraken2",
        prf = kkXpfx + ".k2.prof",
        bed = kkXpfx + ".k2.bead"
    threads: 8
    shell:
        "kraken2 --db {input.db} --threads {threads} --paired {input.fa1} {input.fa2} "
        "--report {output.rpt} | tee {output.out} | metabbq IO kk2prf -d {input.db} -i - -o {output.prf}"
#        " | tee {output.raw} | metabbq IO kk2prf -d {params} -i - -o {output.prf}"
#"kraken2 --db=SAM/SBR/KRAKEN/db --threads 4 --report ${SAMPLE}/KRAKEN/reads2SBR.kreport2 \
#--paired ${SAMPLE}/clean/fastp.sort.1.fq.gz ${SAMPLE}/clean/fastp.sort.2.fq.gz > ${SAMPLE}/KRAKEN/reads2SBR.kraken2"
rule KRAKEN_X_custom_kreport:
    input:
        bed = kkXpfx + ".k2.bead",
        db  = "{sample}/KRAKEN/" + str(config["kk2_db"])
    output:
        rpt = kkXpfx + ".bead.kreport2"
    threads: 8
    shell:
        "awk '$5>0.9&&$8<10' {input.bed} | metabbq kraken-report --d={input.db} - -v > {output.rpt}"

rule KRAKEN_X_sum2species:
    input:
        report = kkXpfx + ".bead.kreport2",
        db  = "{sample}/KRAKEN/" + str(config["kk2_db"])
    output:
        profile = kkXpfx + ".bead.bracken",
        report = kkXpfx + ".bead.breport2"
    threads: 1
    shell:
        "bracken -d {input.db} -i {input.report} -o {output.profile} -l S > {output.report}"
#        "metabbq IO sumprf -r species -n {input.db}/taxonomy/names.dmp -i {input.prf} -o {output}"
