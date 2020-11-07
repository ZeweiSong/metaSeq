#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beads Isolate Assembly Solution
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

configfile: "config.yaml"

rule BIA_0_cutoffOLD:
    input:
        idx = "{sample}/clean/fastp.sort.1.fq.idx",
        bbs = "{sample}/clean/BB.stat"
    output: "{sample}/Assemble_BI/ID.lst.old"
    params:
        p_rpb  = config["p_rpb_min"]
    shell:
        "rpb3k=`awk '($3>3000){{print $2;exit}}' {input.bbs}`\n"
        "awk -v p3=$rpb3k 'BEGIN{{if(p3<{params.p_rpb}){{rpb=p3}}else{{rpb={params.p_rpb}}}}}(!/0000/&&(($3-$2+1)/4)>rpb){{print FNR\"\\t\"$1\"\\t\"($3-$2+1)/4}}' {input.idx} > {output}"

rule BIA_0_sketchOLD:
    input:
        id = "{sample}/clean/fastp.sort.1.fq.idx",
        x1 = "{sample}/clean/fastp.sort.1.fq",
        x2 = "{sample}/clean/fastp.sort.2.fq",
        bb = "{sample}/clean/BB.stat"
    output:  "{sample}/mash/allBeads.msh"
    params:
        pfx = "{sample}/mash/allBeads",
        k   = config["p_dist_k"],
        s   = config["p_dist_s"]
    resources:
        mem=120
    threads: 6
    shell:
        "mash sketch -p {threads} -k {params.k} -s {params.s} -r -B {input.x1} {input.x2} -o {params.pfx}"

rule BIA_0_cutoff:
	input:
		id = "{sample}/clean/fastp.sort.1.fq.idx",
		x1 = "{sample}/clean/fastp.sort.1.fq",
		x2 = "{sample}/clean/fastp.sort.2.fq",
		bb = "{sample}/clean/BB.stat"
	output:
		x1 = "{sample}/mash/BI.1.fq",
		x2 = "{sample}/mash/BI.2.fq"
	params:
		pfx = "{sample}/mash/BI",
		minR= config['p_cluster_minR'],
		maxR= config['p_cluster_maxR'],
		topB= config['p_cluster_topB'],
		ranP= config['p_cluster_ranP'],
		k   = config["p_dist_k"],
		s   = config["p_dist_s"]
	threads: 2
	shell:
		"export maxC={params.maxR}\nexport minC={params.minR}\n"
		"echo choose minc = $minC , maxc = $maxC \n"
		"metabbq binWrite fqpick -x {input.id} -c {params.minR} -m {params.maxR} -b {params.topB} -r {params.ranP} -i {input.x1} -o {output.x1} & \n"
		"metabbq binWrite fqpick -x {input.id} -c {params.minR} -m {params.maxR} -b {params.topB} -r {params.ranP} -i {input.x2} -o {output.x2} && wait\n"

rule BIA_0_sketch:
	input:
		x1 = "{sample}/mash/BI.1.fq",
		x2 = "{sample}/mash/BI.2.fq",
	output:  "{sample}/mash/BI.msh"
	params:
		pfx = "{sample}/mash/BI",
		k   = config["p_dist_k"],
		s   = config["p_dist_s"]
	resources:
		mem=120
	threads: 36
	shell:
		"mash sketch -p 8 -k {params.k} -s {params.s} -r -B {input.x1} {input.x2} -o {params.pfx}"

rule BIA_1_mashstat:
    input: "{sample}/mash/BI.msh"
    output: "{sample}/mash/BI.msh.tsv"
    resources:
        mem=80
    threads: 36
    shell:
        "mash info -t {input}| perl -ne 'if($_=~/\d+\\t(\d+)\\t(\S+)\\t\[(\d+) /){{print \"$2\\t$3\\t$1\\n\"}}' > {output}\n"

rule BIA_1_getBIs:
    input: "{sample}/mash/BI.msh.tsv"
    output: "{sample}/Assemble_BI/ID.lst"
    params:
        k   = config["p_asm_minKmers"],
        n   = (config["i_rdLen"] - config["p_dist_k"] + 1) * config["i_seqPair"],
        c   = config["p_asm_minKmerCoverage"]
    shell:
        "awk '$1!~/0000/&&$3>{params.k}&&$2*{params.n}/$3>{params.c}{{print FNR\"\\t\"$0}}' {input} > {output}"

rule BIA_1_write:
    input:
        s1 = "{sample}/mash/BI.1.fq",
        s2 = "{sample}/mash/BI.2.fq",
        bi = "{sample}/Assemble_BI/ID.lst"
    output: "{sample}/Assemble_BI/log"
    params:
        outDir = "{sample}/Assemble_BI"
    threads: 2
    shell:
        "metabbq beadsWrite.pl -d {input.bi} -f fq -t {threads} -o {params.outDir} -v "
        "--r1 {input.s1}  --r2 {input.s2} &> {output}\n"

rule BIA_1_prepCFG:
    input: "{sample}/Assemble_BI/ID.lst"
    output: "{sample}/primers.cfg"
    params:
        Ad = config["AdRCA"],
        FWD= config["AdFw"],
        REV= config["AdRv"]
    shell:
        """
echo -e "Ad={params.Ad}\nFWD={params.FWD}\nREV={params.REV}" > {output}
        """

rule BIA_2_initASMsh:
    input:
        log = "{sample}/Assemble_BI/log",
        bi = "{sample}/Assemble_BI/ID.lst",
        cfg= "{sample}/primers.cfg"
    output: "{sample}/batch.assemble.BI.sh"
    params:
        samDir = "{sample}",
        outDir = "{sample}/Assemble_BI",
        mode   = config["method"]["assemble"]["mode"],
        mem    = config["p_asm_mem"],
        cpu    = config["p_asm_cpu"]
    shell:
        "for i in `cut -f1 {input.bi}`; do "
        "echo metabbq RCAasm.sh {params.mode} {params.samDir} BI $i BI {params.mem} {params.cpu} {input.cfg}; done > {output}\n"

outPfx="{sample}/summary.BI." + str(config["method"]["assemble"]["mode"])
outXa= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.fasta"
outXn= "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.anno"

outXc= outPfx + ".clip.fasta"

rule BIA_3_assemble:
    input: "{sample}/batch.assemble.BI.sh"
    output: "{sample}/log/batch.assemble.BI.done"
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        dir  = "{sample}"
    threads: config['divideSH'] * config["p_asm_cpu"]
    shell:
        """
# 0. Run assembly one by one
dl=$[$(wc -l {input}|tr ' ' '\\n' |head -1) / {params.divide}]
split -d -l $dl {input} {params.dpfx}
for i in {params.dpfx}??;do sh $i & done > {params.dir}/running.BI.asm.log && wait
mv {params.dir}/running.BI.asm.log {params.dir}/log/batch.assemble.BI.done
        """

rule BIA_4_summary1:
    input: "{sample}/log/batch.assemble.BI.done"
    output:
        fa   = outXa
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        outDir = "{sample}/Assemble_BI",
        mAsm   = config["p_asm_min"],
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
# 1. pick robust assemble fasta
for i in `ls {params.outDir}|grep BI`;do
  if [ -s {params.outDir}/$i/{params.mode}/final.contigs.fa ]; then
    metabbq clusterHelper asmPick -l {params.mAsm} -b $i -i {params.outDir}/$i/{params.mode}/final.contigs.fa
  fi
done > {output.fa}
        """
outXs = "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".contig.tsv"
rule BIA_4_summary2:
    input: "{sample}/log/batch.assemble.BI.done"
    output:
        stat = outXs
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        outDir = "{sample}/Assemble_BI",
        mAsm   = config["p_asm_min"],
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
# 2. stat assemble results
for i in `ls {params.outDir}|grep BI`;do
  if [ -s {params.outDir}/$i/{params.mode}/final.contigs.fa ]; then
    awk -v bc=$i '/^>/{{print bc"\\t"$0}}' {params.outDir}/$i/{params.mode}/final.contigs.fa | sed 's/>//;s/ flag=/\\t/;s/ multi=/\\t/;s/ len=/\\t/'
  fi
done > {output.stat}
        """


rule BIA_4_summary3:
    input: "{sample}/log/batch.assemble.BI.done"
    output:
        fa   = outXc
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        outDir = "{sample}/Assemble_BI",
        mAsm   = config["p_asm_min"],
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
# 3. pick RCA clip fasta
for i in `ls {params.outDir}|grep BI`;do
  if [ -s {params.outDir}/$i/{params.mode}/RCAclip.fa ]; then
    awk -v b=$i 'FNR%2==1{{sub(">",">"b"_",$0);id=$0}}FNR%2==0{{if(length($0)>={params.mAsm}){{print id"\\n"$0}}}}' {params.outDir}/$i/{params.mode}/RCAclip.fa
  fi
done > {output.fa}
        """
outXctsv = "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".clip.metadata.tsv"
rule BIA_4_summary4:
    input: "{sample}/log/batch.assemble.BI.done"
    output:
        tsv   = outXctsv
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        outDir = "{sample}/Assemble_BI",
        mAsm   = config["p_asm_min"],
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
# 4. pick all RCA clip info
for i in `awk '{{printf("BI%08d\\n",$1)}}' {params.outDir}/ID.lst`;do
  if [ -s {params.outDir}/$i/{params.mode}/RCAclip.log ]; then
    awk -v b=$i -F ' +|\\t|_|=' '/multi=/{{print b"\\t"$2"\\t"$3"\\t"$5"\\t"$7"\\t"$9"\\t"$10"\\t"$11"\\t"$12"\\t"$13}}' {params.outDir}/$i/{params.mode}/RCAclip.log
  fi
done > {output.tsv}
        """
outXcafa = "{sample}/summary.BI." + str(config["method"]["assemble"]["mode"]) + ".clip.all.fasta"
rule BIA_4_summary5:
    input: "{sample}/log/batch.assemble.BI.done"
    output:
        fa   = outXcafa
    params:
        divide = config['divideSH'],
        dpfx = "{sample}/sp.BI.",
        outDir = "{sample}/Assemble_BI",
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
# 5. pick RCA clip fasta
for i in `awk '{{printf("BI%08d\\n",$1)}}' {params.outDir}/ID.lst`;do
  if [ -s {params.outDir}/$i/megahit/RCAclip.fa ]; then
    awk -v b=$i 'FNR%2==1{{sub(">",">"b"_",$0);id=$0}}FNR%2==0{{print id"\\n"$0}}' {params.outDir}/$i/{params.mode}/RCAclip.fa
  fi
done > {output.fa}
        """

# 9. clean intermediate files
rule BIA_4_cleanfiles:
    input: "{sample}/log/batch.assemble.BI.done"
    output:
        log   = "{sample}/log/delete.assemble.BI.done"
    params:
        dpfx = "{sample}/sp.BI.",
        outDir = "{sample}/Assemble_BI",
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
# 5. pick RCA clip fasta
for i in `awk '{{printf("BI%08d\\n",$1)}}' {params.outDir}/ID.lst`;do
  if [ -s {params.outDir}/$i ]; then
    rm -rf {params.outDir}/$i
  fi
done > {output.log}
        """

#Following were rules for Bead isolated scaffolding method:
## 1 discard abnormal contigs and rename sequence id:
rule BIS_0_filter:
    input: outXc
    output: "{sample}/CLIP/RCA.clip.rename.fasta"
    shell:
        "metabbq scafHelper ffa < {input} > {output}"

rule BIS_1_vsearch_cluster:
    input: "{sample}/CLIP/RCA.clip.rename.fasta"
    output:
        fa="{sample}/CLIP/RCA.clip.rename.clust.fa",
        uc="{sample}/CLIP/RCA.clip.rename.clust.uc",
        log="{sample}/CLIP/RCA.clip.rename.clust.log"
    params:
        pct = config['p_BIS_clust_id']
    threads: 16
    shell:
        "vsearch --threads {threads} --cluster_fast  {input} --strand both --fasta_width 0 "
        "--centroids {output.fa} --uc  {output.uc} --id {params.pct} --iddef 0 &> {output.log}"

rule BIS_2_blastn_makedb:
    input: "{sample}/CLIP/RCA.clip.rename.clust.fa"
    output:
        nsq="{sample}/CLIP/RCA.clip.rename.clust.fa.nsq",
        nin="{sample}/CLIP/RCA.clip.rename.clust.fa.nin",
        nhr="{sample}/CLIP/RCA.clip.rename.clust.fa.nhr"
    shell:
        "makeblastdb -dbtype nucl -in {input}"

rule BIS_2_blastn_mappingEachOther:
    input:
        fa = "{sample}/CLIP/RCA.clip.rename.clust.fa",
        db = "{sample}/CLIP/RCA.clip.rename.clust.fa.nsq",
    output: "{sample}/CLIP/RCA.clip.rename.clust.blastn.m6"
    params:
        pct = 99
    threads: 8
    shell:
        "blastn -num_threads {threads} -db {input.fa} -query {input.fa} -out {output} "
        "-perc_identity {params.pct}  -outfmt '6 std qlen slen' "

rule BIS_2_AdjGet_scaffolding:
    input:
        fa = "{sample}/CLIP/RCA.clip.rename.clust.fa",
        uc = "{sample}/CLIP/RCA.clip.rename.clust.uc",
        m6 = "{sample}/CLIP/RCA.clip.rename.clust.blastn.m6"
    output: "{sample}/CLIP/RCA.clip.rename.clust.ADJ.log"
    params:
        pfx = "{sample}/CLIP/RCA.clip.rename.clust.adj"
    shell:
        "metabbq AdjGet -i {input.fa} -u {input.uc} -m {input.m6} -o {params.pfx} -v &> {output}"

rule BIS_3_Scaffold_clust:
    input: "{sample}/CLIP/RCA.clip.rename.clust.adj.merge.fa"
    output:
        fa="{sample}/CLIP/RCA.clip.rename.clust.adj.merge.clust.fa",
        uc="{sample}/CLIP/RCA.clip.rename.clust.adj.merge.clust.uc",
        log="{sample}/CLIP/RCA.clip.rename.clust.adj.merge.clust.log"
    params:
        pct = config['p_BIS_clust_id']
    threads: 16
    shell:
        "sed 's/ /_/g' {input} > {input}.tmp.fa \n"
        "vsearch --threads {threads} --cluster_fast {input}.tmp.fa --strand both --fasta_width 0 --relabel_keep "
        "--centroids {output.fa} --uc  {output.uc} --id {params.pct} --iddef 0 &> {output.log}"

rule BIS_4_Scaffold_quast:
    input: "{sample}/CLIP/RCA.clip.rename.clust.adj.merge.clust.fa"
    output: "{sample}/CLIP/RCA.clip.rename.clust.adj.merge.clust_quast.log"
    params:
        dir = "{sample}/CLIP/RCA.clip.rename.clust.adj.merge.clust_quast"
    threads: 16
    shell:
        "quast.py -m 3000 --min-identity 98 -t {threads} {input} -R $LFR/Source/REF/zymo/D6305.fix.rRNA.fa -o {params.dir} &> {output}"
