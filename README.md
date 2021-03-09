# METAgenomic Beads Barcoding Quantification (METABBQ) pipeline.

This is a data processing pipeline to achieve bacterial/fungal long amplicons from complex environmental samples. Two experiments were mainly implemented: *sing-tube Long Fragment Reads (stLFR)*  and *Rolling Circle Replication (RCR)*.

### Installation 

**Prerequisites**  
- python >= 3.6
- perl >= 5
- [fastp](https://github.com/Scelta/fastp/tree/stlfr) - A modified version which implemented a module to split the stlfr barcodes.
- [Mash](https://github.com/marbl/Mash) ([dev repo](https://biogit.cn/PUB/Mash) ) - dev version mandatory since I've modified it to fit stLFR data 
- [Snakemake](https://bitbucket.org/snakemake/snakemake) - a pythonic workflow system.  
- [blast](https://blast.ncbi.nlm.nih.gov) - The classic alignment tool finding regions of similarity between biological sequences.
- **Assemble methods**  
  - [SPAdes ](https://github.com/ablab/spades) - SPAdes Genome Assembler
  - [MEGAHIT](https://github.com/voutcn/megahit) - An ultra-fast and memory-efficient NGS assembler

I recommend to install above tools in an virtual env via [conda](https://conda.anaconda.org/):
1. create and install part of them:
```bash
conda create -n metaseq -c bioconda -c conda-forge snakemake pigz megahit blast
source activate metaseq
```
2. According to the corresponding documents, install `fastp`, `SPAdes` and `community`, etc. under env `metaseq`

Make sure above commands (executables) can be found in the `PATH`. 

**install**  
Clone the launcher to initiate the work dir as well as to call sub-functions.
```
cd /path/to/your/dir
git clone https://github.com/ZeweiSong/metaSeq.git
export PATH="/path/to/your/dir/metaSeq":$PATH
```

I haven't yet write any testing module to check abve prerequesites. At present you may need to test it yourself.

### Usage
**Prepare configs**
```bash
cd instance
metabbq cfg  
```
This command will create a `default.cfg` in your current dir.
You should modifed it to let the launcher know the required files and parameters

**Initiating a project**
Prepare an `input.list` file to describe the sample name and input sequence file path.  
```
metabbq -i input.list -c default.cfg -V
```
By default, the `metabbq` will create a directory with the name of {sample} and a sub-directory named `input` under it.  

#### Run Quality-Contorl module
```
metabbq smk -j -np {sample}/clean/BB.stat
# -j make the jobs execuated paralled under suitable cores/threads
# -n mean dry-run with a preview of "what needs to be run". Remove it to really run the pipeline.
```

#### Run pre-binning assembly module
You need to select a assemble tool in the configure file and the corresponding output file name in following:
```
metabbq smk -j -np {sample}/summary.BI.megahit.contig.fasta
metabbq smk -j -np {sample}/summary.BI.idba.contig.fasta
metabbq smk -j -np {sample}/summary.BI.spades.contig.fasta
```

### Troubleshooting  
Feedback are welcome to submit in the issue page.
