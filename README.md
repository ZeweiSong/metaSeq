# MetaSeq

This is  a sequencing data treatment pipeline mainly implemented with stLFR technology.

## Requirment

**Environment**: `python >= 3.6` `perl >= 5` `R3.4`

**Developing projects**：

> Note: the biogit is an internal site and only accessable from intranet at present.

**metaSeq** ( [github](https://github.com/ZeweiSong/metaSeq) | [biogit](https://biogit.cn/Fangchao/metaSeq) )

**cOMG** ( [biogit](https://biogit.cn/Fangchao/Omics_pipeline) )

**fastp** ( [github](https://github.com/OpenGene/fastp) | [biogit](https://biogit.cn/PUB/fastp) )

**Mash**( [github](https://github.com/marbl/Mash) | [biogit](https://biogit.cn/PUB/Mash) )

**Community ** ( [source](https://sites.google.com/site/findcommunities/) | [biogit](https://biogit.cn/PUB/community) )

> make sure  above commands can be found in the PATH

**Third party program:**

- **Snakemake** - a pythonic workflow system ([bitbucket](https://bitbucket.org/snakemake/snakemake))

- **SPAdes** - SPAdes Genome Assembler ( [about ](http://cab.spbu.ru/software/spades/)| [github ](https://github.com/ablab/spades) )

  - **configure**

    ```bash
    #SPAdes 3.13.0
    export PATH="$MOPT/SPAdes-3.13.0-Linux/bin":$PATH
    ```

## Prepare
**Deploy pipeline**

```
cd /path/to/work/dir
git clone https://github.com/ZeweiSong/metaSeq.git
```

**Show pipeline directed acyclic graph(dag)**

```bash
cd metaSeq/instance
snakemake --dag | dot -Tsvg > dag.svg
```

**test stLFR process**

```
snakemake -j -rp benchmarks/stLFR_summary.txt
# -j make the jobs execuated paralled under suitable cores/threads
```