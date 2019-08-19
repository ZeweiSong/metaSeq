#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beads Cluster Generation
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

rule BC2M_0_sketch:
    input:
        id = "{sample}/clean/fastp.sort.1.fq.idx",
        x1 = "{sample}/clean/fastp.sort.1.fq",
        x2 = "{sample}/clean/fastp.sort.2.fq",
        bb = "{sample}/clean/BB.stat"
    output:  "{sample}/BC2M/mash/bMin2.msh"
    params:  "{sample}/BC2M/mash/bMin2"
    shell:
        "export maxC=$[ `awk '$3>=10000{{print $2;exit}}' {input.bb}` ]\n"
        "export minC=$[ `awk '$3>=20000{{print $2;exit}}' {input.bb}` ]\n"
        "echo choose minc = $minC , maxc = $maxC \n"
        "metabbq binWrite fqpick -x {input.id} -c $minC -m $maxC -i {input.x1} -o {params}.sort.1.fq & \n"
        "metabbq binWrite fqpick -x {input.id} -c $minC -m $maxC -i {input.x2} -o {params}.sort.2.fq & \n"
        "wait && mash sketch -r -B {params}.sort.1.fq {params}.sort.2.fq -o {params}"

rule BC2M_1_distRaw:
    input:  "{sample}/BC2M/mash/bMin2.msh"
    output: "{sample}/BC2M/mash/bMin2.raw.dist"
    shell:  "mash dist -p 24 -d 0.5 {input} {input}  > {output}"

rule BC2M_2_distReNum:
    input:  "{sample}/BC2M/mash/bMin2.raw.dist"
    output:
        d = "{sample}/BC2M/mash/bMin2.bc.dist",
        m = "{sample}/BC2M/mash/bMin2.bc.map"
    shell:  "metabbq filterMASHoutput.pl -i  {input} -o {output.d} -M  {output.m}"

rule BC2M_3_convert:
    input:   "{sample}/BC2M/mash/bMin2.bc.dist"
    output:
        b  = "{sample}/BC2M/mash/bMin2.bc.bin",
        w  = "{sample}/BC2M/mash/bMin2.bc.w"
    shell:   "convert -i {input} -w {output.w} -o {output.b}"

rule BC2M_4_community:
    input:
        b  = "{sample}/BC2M/mash/bMin2.bc.bin",
        w  = "{sample}/BC2M/mash/bMin2.bc.w"
    output:   "{sample}/BC2M/mash/bMin2.bc.tree"
    shell:   "community {input.b} -w {input.w} -l -1 -v > {output}"

rule BC2M_5_clusterMAP:
    input:
        t = "{sample}/BC2M/mash/bMin2.bc.tree",
        m = "{sample}/BC2M/mash/bMin2.bc.map"
    output: "{sample}/BC2M/mash/bMin2.bc.tree.target.cluster"
    shell:
        "export maxLv=$[ `hierarchy {input.t} | wc -l` - 2 ]\n"
        "hierarchy {input.t} -l $maxLv > {input.t}.Lv$maxLv.lst\n"
        "metabbq change.id.pl -n {input.t}.Lv$maxLv.lst "
        "-m {input.m} -v -o {input.t}.Lv$maxLv.cluster\n"
        "ln -s bMin2.bc.tree.Lv$maxLv.cluster {output}"

rule BC2M_6_stat:
    input:
        i = "{sample}/clean/fastp.sort.1.fq.idx",
        c = "{sample}/BC2M/mash/bMin2.bc.tree.target.cluster"
    output: "{sample}/BC2M/mash/bMin2.bc.tree.target.cluster.count"
    shell:
        "perl -e 'open IN,\"<'{input.i}'\";while(<IN>){{chomp;@a=split;"
        "$num=($a[2]-$a[1]+1)/4;$HASH{{$a[0]}}=$num}};while(<>){{chomp;@a=split;"
        "$HB{{$a[0]}}{{R}}+=$HASH{{$a[1]}};$HB{{$a[0]}}{{C}}++}};"
        "foreach my $c (sort {{$a<=>$b}} keys %HB){{"
        "print \"$c\t$HB{{$c}}{{C}}\t$HB{{$c}}{{R}}\n\"}}' < {input.c} > {output}"
rule BC2M_7_write:
    input:
        s1 = "{sample}/clean/fastp.sort.1.fq",
        s2 = "{sample}/clean/fastp.sort.2.fq",
        bc = "{sample}/BC2M/mash/bMin2.bc.tree.target.cluster"
    output: "{sample}/BC2M/Assemble_mashBC.log"
    params:
        outDir = "{sample}/BC2M/Assemble_mashBC",
        threads = config["threads"]
    shell:
        "metabbq beadsWrite3.pl -b {input.bc} -f fq -t 4 -o {params.outDir} -v "
        "--r1 {input.s1}  --r2 {input.s2} > {output}\n"

rule BC2M_8_initASMsh:
    input:
        log = "{sample}/BC2M/Assemble_mashBC.log",
        bc = "{sample}/BC2M/mash/bMin2.bc.tree.target.cluster.count",
    output: "{sample}/BC2M/batch.assemble.BC.sh"
    params:
        samDir = "{sample}/BC2M",
        outDir = "{sample}/BC2M/Assemble_mashBC",
    shell:
        "for i in `sort -k2,2nr {input.bc} | cut -f1`; do "
        "echo metabbq bcPost.template.sh F {params.samDir} mashBC $i BC;"
        "done > {params.samDir}/batch.assemble.BC.sh\n"
