# MetaSeq (Early access stage)

This is a sequencing data treatment pipeline mainly implemented with stLFR technology.

## Requirment

**Environment**: `python >= 3.6` `perl >= 5` `R3.4`

**Requirements**：

> Note: the biogit is an internal site and only accessable from intranet at present.

**metabbq** ( [github](https://github.com/ZeweiSong/metaSeq) | [biogit](https://biogit.cn/Fangchao/metaSeq) )  
metabbq means "METAgenome Bead Barcode Quantification", which is a launcher to initiate workdir and calling sub functions.

~~**cOMG** ( [biogit](https://biogit.cn/Fangchao/Omics_pipeline) )~~

**fastp** ( [github](https://github.com/OpenGene/fastp) | [biogit](https://biogit.cn/PUB/fastp) )  
I've modified `fastp` to speed up the split barcodes process

~~**Mash**( [github](https://github.com/marbl/Mash) | [biogit](https://biogit.cn/PUB/Mash) )~~

**Community** ([source](https://sites.google.com/site/findcommunities/) | [biogit](https://biogit.cn/PUB/community))

> make sure  above commands can be found in the PATH

**Third party program:**

- **Snakemake** - a pythonic workflow system ([bitbucket](https://bitbucket.org/snakemake/snakemake))

- **SPAdes** - SPAdes Genome Assembler ( [about ](http://cab.spbu.ru/software/spades/)| [github ](https://github.com/ablab/spades) )

  - **configure**

    ```bash
    #SPAdes 3.13.0
    export PATH="$MOPT/SPAdes-3.13.0-Linux/bin":$PATH
    ```

- **MEGAHIT** -  An ultra-fast and memory-efficient NGS assembler ([github](https://github.com/voutcn/megahit))


## Usage

**Deploy pipeline**
```
cd /path/to/your/dir
git clone https://github.com/ZeweiSong/metaSeq.git

export PATH="/path/to/your/dir/metaSeq":$PATH
```
**Prepare configs**
```bash
cd instance
metabbq cfg  
```
This command will create a `default.cfg` in your current dir. 
You should modifed it to let the launcher know the required files and parameters

**Initiating a project**
```
metabbq -i input.list -c default.cfg -V
```
By default, the launcher will initate the workshop in current directory.

**Show pipeline directed acyclic graph(dag)**
```bash
snakemake --dag | dot -Tsvg > dag.svg
```

**test pipeline**

```
snakemake -j -np T01/VSEARCH/read.merge.derep.2T.bc.graph.tree
# -j make the jobs execuated paralled under suitable cores/threads
# -n mean dry-run with a preview of "what needs to be run". Remove it to really run the pipeline.
```
