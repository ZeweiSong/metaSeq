#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beads Cluster Generation
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

if config["method"]["cluster"]["mash"]["do"]:
    rule BC_0_sketch:
        input:
            id = "{sample}/clean/fastp.sort.1.fq.idx",
            x1 = "{sample}/clean/fastp.sort.1.fq",
            x2 = "{sample}/clean/fastp.sort.2.fq"
        output:  "{sample}/mash/bMin2.msh"
        params:  "{sample}/mash/bMin2"
        shell:
            "metabbq binWrite fqpick -x {input.id} -c 2 -i {input.x1} -o {params}.sort.1.fq & \n"
            "metabbq binWrite fqpick -x {input.id} -c 2 -i {input.x2} -o {params}.sort.2.fq & \n"
            "wait && mash sketch -r -B {params}.sort.1.fq {params}.sort.2.fq -o {params}"

    rule BC_1_distRaw:
        input:  "{sample}/mash/bMin2.msh"
        output: "{sample}/mash/bMin2.raw.dist"
        shell:  "mash dist -p 24 -d 0.5 {input} {input}  > {output}"

    rule BC_2_distReNum:
        input:  "{sample}/mash/bMin2.raw.dist"
        output:
            d = "{sample}/mash/bMin2.bc.dist",
            m = "{sample}/mash/bMin2.bc.map"
        shell:  "metabbq filterMASHoutput.pl -i  {input} -o {output.d} -M  {output.m}"

    rule BC_3_convert:
        input:   "{sample}/mash/bMin2.bc.dist"
        output:
            b  = "{sample}/mash/bMin2.bc.bin",
            w  = "{sample}/mash/bMin2.bc.w"
        shell:   "convert -i {input} -w {output.w} -o {output.b}"

    rule BC_4_community:
        input:
            b  = "{sample}/mash/bMin2.bc.bin",
            w  = "{sample}/mash/bMin2.bc.w"
        output:   "{sample}/mash/bMin2.bc.tree"
        shell:   "community {input.b} -w {input.w} -l -1 -v > {output}"

    rule BC_5_clusterMAP:
        input:
            t = "{sample}/mash/bMin2.bc.tree",
            m = "{sample}/mash/bMin2.bc.map"
        output: "{sample}/mash/bMin2.bc.tree.target.cluster"
        shell:
            "export maxLv=$[ `hierarchy {input.t} | wc -l` - 2 ]\n"
            "hierarchy {input.t} -l $maxLv > {input.t}.Lv$maxLv.lst\n"
            "metabbq change.id.pl -n {input}.Lv$maxLv.lst "
            "-m {input.m} -v -o {input.t}.Lv$maxLv.cluster\n"
            "ln -s bMin2.bc.tree.Lv$maxLv.cluster {output}"

    rule BC_6_stat:
        input:
            i = "{sample}/clean/fastp.sort.1.fq.idx",
            c = "{sample}/mash/bMin2.bc.tree.target.cluster"
        output: "{sample}/mash/bMin2.bc.tree.target.cluster.count"
        shell:
            "perl -e 'open IN,\"<'{input.i}'\";while(<IN>){{chomp;@a=split;"
            "$num=($a[2]-$a[1]+1)/4;$HASH{{$a[0]}}=$num}};while(<>){{chomp;@a=split;"
            "$HB{{$a[0]}}{{R}}+=$HASH{{$a[1]}};$HB{{$a[0]}}{{C}}++}};"
            "foreach my $c (sort {{$a<=>$b}} keys %HB){{"
            "print \"$c\t$HB{{$c}}{{C}}\t$HB{{$c}}{{R}}\n\"}}' < {input.c} > {output}"
    rule BC_7_write:
        input:
            s1 = "{sample}/clean/fastp.sort.1.fq",
            s2 = "{sample}/clean/fastp.sort.2.fq",
            bc = "{sample}/mash/bMin2.bc.tree.target.cluster"
        output: "{sample}/Assemble_mashBC.log"
        params:
            outDir = "{sample}/Assemble_mashBC",
            threads = config["threads"]
        shell:
            "metabbq beadsWrite3.pl -b {input.bc} -f fq -t 4 -o {params.outDir} -v "
            "--r1 {input.s1}  --r2 {input.s2} > {output}\n"

    rule BC_8_initASMsh:
        input:
            log = "{sample}/Assemble_mashBC.log",
            bc = "{sample}/mash/bMin2.bc.tree.target.cluster.count",
        output: "{sample}/batch.assemble.BC.sh"
        params:
            samDir = "{sample}",
            outDir = "{sample}/Assemble_mashBC",
        shell:
            "for i in `sort -k2,2nr {input.bc} | cut -f1`; do "
            "echo metabbq bcPost.template.sh F {params.samDir} mashBC $i BC;"
            "done > {params.samDir}/batch.assemble.BC.sh\n"

# Previous method
if config["method"]["cluster"]["dups"]["do"]:
    rule BC_1_method:
        input:
            s1 = "{sample}/clean/fastp.sort.1.fq",
            s2 = "{sample}/clean/fastp.sort.2.fq",
            idx= "{sample}/clean/fastp.sort.1.fq.idx"
        output:
            bm = "{sample}/VSEARCH/read.merge.derep.2T.bc.cluster_Lv1.main",
            bc = "{sample}/VSEARCH/read.merge.derep.2T.bc.cluster_Lv1.count.main",
            bi = "{sample}/VSEARCH/read.individual.beads.list"
        params:
            pfx= "{sample}/VSEARCH/read",
            threads =config["thread"]["vsearch"]
        shell:
            "metabbq template.vsearch.preCluster.sh {params.threads} {input.s1} {input.s2} {params.pfx}\n"

    rule BC_2_write:
        input:
            s1 = "{sample}/clean/fastp.sort.1.fq",
            s2 = "{sample}/clean/fastp.sort.2.fq",
            bc = "{sample}/VSEARCH/read.merge.derep.2T.bc.cluster_Lv1.main",
            bi = "{sample}/VSEARCH/read.individual.beads.list"
        output: "{sample}/Assemble_Lv1.log"
        params:
            outDir = "{sample}/Assemble_Lv1"
        shell:
            "metabbq beadsWrite3.pl -b {input.bc} -d {input.bi} -f fq -t 4 -o {params.outDir} -v "
            "--r1 {input.s1}  --r2 {input.s2} > {output}\n"

    rule BC_3_Assemble_init_sh:
        input:
            log = "{sample}/Assemble_Lv1.log",
            bc = "{sample}/VSEARCH/read.merge.derep.2T.bc.cluster_Lv1.count.main",
            bi = "{sample}/VSEARCH/read.individual.beads.list"
        output: "{sample}/batch.assemble.BC.sh"
        params:
            samDir = "{sample}",
            outDir = "{sample}/Assemble_Lv1"
        shell:
            "for i in `sort -k2,2nr {input.bc} | cut -f1`; do "
            "echo metabbq bcPost.template.sh B {params.samDir} 1 $i BC;"
            "done > {params.samDir}/batch.assemble.BC.sh\n"
            "for i in `sort -k2,2nr {input.bi} | cut -f1`; do "
            "echo metabbq bcPost.template.sh B {params.samDir} 1 $i BI;"
            "done > {params.samDir}/batch.assemble.BI.sh\n"
