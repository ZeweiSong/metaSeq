# MetaSeq

This is  a sequencing data treatment pipeline mainly implemented with stLFR technology.

## Requirment

**Environment**: `python >= 3.6` `perl >= 5` `R3.4`

**Developing projects**：[metaSeq](https://github.com/ZeweiSong/metaSeq) ([biogit](https://biogit.cn/Fangchao/metaSeq))、[cOMG](https://biogit.cn/Fangchao/Omics_pipeline)、[fastp](https://github.com/OpenGene/fastp) ([biogit](https://biogit.cn/PUB/fastp))

> make sure `metaSeqinit`, `cOMG` and `fastp` canbe found in the PATH

**Third party program:**

- **Snakemake** - a pythonic workflow system ([bitbucket](https://bitbucket.org/snakemake/snakemake))

## Prepare
**init pipeline**
**configure**
**Show pipeline directed acyclic graph(dag)**
```
snakemake --dag | dot -Tsvg > dag.svg
```

**test stLFR process**
```
snakemake -j -rp benchmarks/stLFR_summary.txt
# -j make the jobs execuated paralled under suitable cores/threads
```
