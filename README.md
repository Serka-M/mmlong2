<p align="center">
<img align="center" width="250" height="250" src="msc/mmlong2-logo.png" alt="logo" style="zoom:100%;" />
</p>

Automated long-read metagenomics workflow, using either PacBio HiFi or Nanopore sequencing reads as input to generate characterized MAGs.
The mmlong2 workflow is a continuation of [mmlong](https://github.com/SorenKarst/mmlong).
<br/>

## Workflow description
### Core features
* [Snakemake](https://snakemake.readthedocs.io) workflow running dependencies from a [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) container for enhanced reproducibility
* Bioinformatics tool and parameter optimizations for high complexity metagenomics samples
* Circular prokaryotic genome extraction as separate genome bins
* Eukaryotic contig removal for reduced prokaryotic genome contamination
* Differential coverage support for improved prokaryotic genome recovery
* Iterative ensemble binning strategy for improved prokaryotic genome recovery
* Genome quality classification according to [MIMAG guidelines](https://www.nature.com/articles/nbt.3893)
* Expanded prokaryotic genome quality assessment, including microdiversity approximation and chimerism checks
* Taxonomic classification at prokaryotic genome, contig and 16S rRNA levels
* Generation of analysis-ready dataframes

### Schematic overview
<img src="msc/mmlong2-np-wf.png" alt="mmlong2-np" style="zoom:100%;" />
<br/>

## Installation (Bioconda)
### Bioconda
The recommended way of installing mmlong2 is by setting up a [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) environment through [Bioconda](https://bioconda.github.io/):
```
mamba install -c bioconda mmlong2
```

### From source (Conda) 
A Conda environment with the latest workflow code can also be created by using the following code:
```
mamba create --prefix mmlong2 -c conda-forge -c bioconda snakemake=8.2.3 singularity=3.8.6 zenodo_get pv pigz tar yq ncbi-amrfinderplus -y
mamba activate ./mmlong2 || source activate ./mmlong2
git clone https://github.com/Serka-M/mmlong2 mmlong2/repo
cp -r mmlong2/repo/src/* mmlong2/bin
chmod +x mmlong2/bin/mmlong2
mmlong2 -h 
```

### Databases and bioinformatics software
Bioinformatics tools and other software dependencies will be automatically installed when running the workflow for the first time.
By default, a pre-built [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) container will be downloaded and set up, although pre-defined Conda environments can also be used by running the workflow with the '--conda_envs_only' setting.
<br/>
To acquire prokaryotic genome taxonomy and annotation results, databases are necessary and can be automatically installed by running the following command:
```
mmlong2 --install_databases
```
If some of the databases are already installed, they can also be used by the workflow without a new download (e.g. '--database_gtdb' option). Alternatively, a guide for [manual](msc/mmlong2-db.md) database installation is also provided.
<br/>

**Quick-start (AAU bioserver users):**
```
conda activate /projects/microflora_danica/mmlong2/conda/mmlong2-v0.9.2
mmlong2 -h
```
<br/>

**Usage example for Nanopore-only mode:**
```
mmlong2 -np [Nanopore_reads.fastq] -p [Processes/Threads] -o [Output_dir]
```
<br/>

**Full usage:**
```
MAIN INPUTS:
-np     --nanopore_reads        Path to Nanopore reads (default: none)
-pb     --pacbio_reads          Path to PacBio HiFi reads (default: none)
-o      --output_dir            Output directory name (default: mmlong2)
-p      --processes             Number of processes/multi-threading (default: 3)
-cov    --coverage              CSV dataframe for differential coverage binning (e.g. NP/PB/IL,/path/to/reads.fastq)
-run    --run_until             Run pipeline until a specified stage completes 
				(e.g. assembly polishing binning taxonomy annotation variants)

ADDITIONAL INPUTS:
-tmp    --temporary_dir         Directory for temporary files (default: none)
-med1   --medaka_model_polish   Medaka polishing model (default: r1041_e82_400bps_sup_v4.2.0)
-med2	--medaka_model_variant	Medaka variant calling model (default: r1041_e82_400bps_sup_variant_v4.2.0)
-sem    --semibin_model         Binning model for SemiBin (default: global)
-fmo    --flye_min_ovlp         Minimum overlap between reads used by Flye assembler (default: auto)
-fmc    --flye_min_cov          Minimum initial contig coverage used by Flye assembler (default: 3)
-mlc    --min_len_contig        Minimum assembly contig length (default: 3000)
-mlb    --min_len_bin           Minimum genomic bin size (default: 250000)
-slv    --silva                 Silva database to use (default: none)
-mds    --midas                 Midas database to use (default: none)
-gnc    --gunc                  Gunc database to use (default: none)
-bkt    --bakta                 Bakta database to use (default: none)
-kj     --kaiju                 Kaiju database to use (default: none)
-gdb    --gtdb                  GTDB-tk database to use (default: none)
-x1     --extra_inputs1         Extra inputs for the MAG production part of the Snakemake workflow (default: none)
-x2     --extra_inputs2         Extra inputs for the MAG processing part of the Snakemake workflow (default: none)

MISCELLANEOUS INPUTS:
-h      --help                  Print help information
-v      --version               Print workflow version number
```
<br/>

**Overview of result files:**
* `assembly.fasta` - assembled and polished metagenome
* `rRNA.fa` - rRNA sequences, recovered from the polished metagenome
* `rRNA_16S.fa` - 16S rRNA sequences, recovered from the polished metagenome
* `<name>_contigs.tsv` - dataframe for metagenome contig metrics
* `<name>_bins.tsv` - dataframe for automated binning results
* `<name>_general.tsv` - workflow results, summarized into a single row
* `dependencies.csv`- list of dependencies used and their versions
* `bins` - directory for metagenome assembled genomes
* `bakta` - directory, containing bin annotation results from [bakta](https://github.com/oschwengers/bakta)
<br/>

**Additional documentation:**
* [Dataframe description](msc/mmlong2-dfs.md)
* [Dependency list](msc/mmlong2-dep.md)
* [Database setup](msc/mmlong2-db.md)
<br/>

**Comments:**
* The workflow assumes that the input reads have been quality-filtered and adapter/barcode sequences have been trimmed off.
* The workflow is long-read-based and requires either Nanopore or PacBio HiFi reads. It doesn't feature an Illumina-only mode.
* If the workflow crashes, it can be resumed by re-running the same command. Some of the intermediary files might have to be removed for compatibility.
* It is recommended to run the workflow from a screen session. This can be achieved with e.g. `screen -R mmlong2` and then running the workflow.
* Multiple large-scale databases are utilized by mmlong2 for genome bin analysis. If you are only interested in getting the MAGs, check out [mmlong2-lite](https://github.com/Serka-M/mmlong2-lite).
<br/>

**Future improvements**<br/>
Suggestions on improving the workflow or fixing bugs are always welcome.<br/>
Please use the GitHub `Issues` section or e-mail to mase@bio.aau.dk for providing feedback.
