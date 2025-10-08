<p align="center">
<img align="center" width="250" height="250" src="msc/mmlong2-logo.png" alt="logo" style="zoom:100%;" />
</p>

Genome-centric long-read metagenomics workflow for automated recovery and analysis of prokaryotic genomes with Nanopore or PacBio HiFi sequencing data.
The mmlong2 workflow is a continuation of [mmlong](https://github.com/SorenKarst/mmlong).
<br/>

## Workflow description
### Core features
* [Snakemake](https://snakemake.readthedocs.io) workflow running dependencies from a [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) container for enhanced reproducibility
* Bioinformatics tool and parameter optimizations for processing high complexity metagenomic samples
* Circular prokaryotic genome extraction as separate genome bins
* Eukaryotic contig removal for reduced prokaryotic genome contamination
* Differential coverage support for improved prokaryotic genome recovery
* Iterative ensemble binning strategy for improved prokaryotic genome recovery
* Recovered genome quality classification according to [MIMAG guidelines](https://www.nature.com/articles/nbt.3893)
* Supplemental genome quality assessment, including microdiversity approximation and chimerism checks
* Automated taxonomic classification at genome, contig and 16S rRNA levels
* Generation of analysis-ready [dataframes](msc/mmlong2-dfs.md) at genome and contig levels

### Schematic overview
<img src="msc/mmlong2-np-wf.png" alt="mmlong2-np" style="zoom:100%;" />

## Installation
The mmlong2 pipeline is currently only available on Linux and was developed, tested on Ubuntu 22.04.

### Bioconda
The recommended way of installing mmlong2 is by setting up a [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) environment through [Bioconda](https://bioconda.github.io/) (ETA 20 minutes):
```
conda create -c conda-forge -c bioconda mmlong2 -n mmlong2
```

### From source (Conda) 
Alternatively, a local Conda environment with the latest workflow code can also be created by using the following code (ETA 20 minutes):
```
conda create --prefix mmlong2 -c conda-forge -c bioconda snakemake=9.12.0 singularity=3.8.6 zenodo_get=1.6.1 pv=1.6.6 pigz=2.8 tar=1.35 rsync=3.4.1 yq=3.4.3 ncbi-amrfinderplus=4.0.23 -y
conda activate ./mmlong2 || source activate ./mmlong2
git clone https://github.com/Serka-M/mmlong2 mmlong2/repo
cp -r mmlong2/repo/src/* mmlong2/bin
chmod +x mmlong2/bin/mmlong2
mmlong2 -h 
```

### Databases and bioinformatics software
Bioinformatics tools and other software dependencies will be automatically installed when running the workflow for the first time (ETA 40 minutes).
By default, a pre-built [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) container will be downloaded and set up, although pre-defined Conda environments can also be used by running the workflow with the `--conda_envs_only` setting.
<br/><br/>
To acquire prokaryotic genome taxonomy and annotation results, databases are necessary and can be automatically installed by running the following command:
```
mmlong2 --install_databases
```
If some of the databases are already installed, they can also be re-used by the workflow without downloading (e.g. `--database_gtdb` option). Alternatively, a guide for [manual](msc/mmlong2-db.md) database installation is also provided.
<br/><br/>

## Running mmlong2
### Usage examples
When running mmlong2 for the first time, it is highly recommended to first perform a testrun with a small Nanopore or PacBio read dataset. Such example datasets are available at [Zenodo](https://zenodo.org/records/12168493). When the mmlong2 environment is active, these can be downloaded by running:
```
zenodo_get -r 12168493
```
For the initial run of the pipeline, please launch only one instance of mmlong2, as to not interfere with the automated dependency installation (ETA 40 minutes).
<br/><br/><br/>
To test the workflow in Nanopore mode up until the genome binning completes (ETA 2 hours, 110 Gb peak RAM):
```
mmlong2 -np mmlong2_np.fastq.gz -o mmlong2_testrun_np -p 60 -run binning
```
Once the run successfully finishes, a directory `mmlong2_testrun_np` can be expected in the current working directory. Inside, a `results` directory with the main pipeline output (described [here](https://github.com/Serka-M/mmlong2/tree/main?tab=readme-ov-file#overview-of-workflow-results)) can be found.
Users can expect 19 MAGs to be recovered (4 circular, 1 single-contig, 14 multi-contig).
<br/><br/><br/>
To test the workflow in PacBio HiFi mode using [metaMDBG](https://www.nature.com/articles/s41587-023-01983-6) as the assembler and perform genome recovery and analysis (ETA 4.5 hours, 170 Gb peak RAM):
```
mmlong2 -pb mmlong2_pb.fastq.gz -o mmlong2_testrun_pb -p 60 -dbg
```
A total of 31 MAGs (3 circular, 6 single-contig, 22 multi-contig) should be recovered.
<br/><br/>
### Full usage
```
MAIN INPUTS:
-np     --nanopore_reads        Path to Nanopore reads
-pb     --pacbio_reads          Path to PacBio HiFi reads
-o      --output_dir            Output directory name (default: mmlong2)
-p      --processes             Number of processes/multi-threading (default: 3)

OPTIONAL SETTINGS:
-db     --install_databases     Install missing databases used by the workflow
-dbd    --database_dir          Output directory for database installation (default: current working directory)
-cov    --coverage              CSV dataframe for differential coverage binning (e.g. NP/PB/IL,/path/to/reads.fastq)
-run    --run_until             Run pipeline until a specified stage completes (e.g.  assembly polishing filtering singletons coverage binning taxonomy annotation extraqc stats)
-tmp    --temporary_dir         Directory for temporary files (default: current working directory)
-dbg    --use_metamdbg          Use metaMDBG for assembly of PacBio reads (default: use metaFlye)
-med    --medaka_model          Medaka polishing model (default: r1041_e82_400bps_sup_v5.0.0)
-mo     --medaka_off            Do not run Medaka polishing with Nanopore assemblies (default: use Medaka)
-vmb    --use_vamb              Use VAMB for binning (default: use GraphMB)
-sem    --semibin_model         Binning model for SemiBin (default: global)
-mlc    --min_len_contig        Minimum assembly contig length (default: 3000)
-mlb    --min_len_bin           Minimum genomic bin size (default: 250000)
-rna    --database_rrna         16S rRNA database to use
-gunc   --database_gunc         Gunc database to use
-bkt    --database_bakta        Bakta database to use
-kj     --database_kaiju        Kaiju database to use
-gtdb   --database_gtdb         GTDB-tk database to use
-h      --help                  Print help information
-v      --version               Print workflow version number

ADVANCED SETTINGS:
-fmo    --flye_min_ovlp         Minimum overlap between reads used by Flye assembler (default: auto)
-fmc    --flye_min_cov          Minimum initial contig coverage used by Flye assembler (default: 3)
-env    --conda_envs_only       Use conda environments instead of container (default: use container)
-n      --dryrun                Print summary of jobs for the Snakemake workflow
-t      --touch                 Touch Snakemake output files
-r1     --rule1                 Run specified Snakemake rule for the MAG production part of the workflow
-r2     --rule2                 Run specified Snakemake rule for the MAG processing part of the workflow
-x1     --extra_inputs1         Extra inputs for the MAG production part of the Snakemake workflow
-x2     --extra_inputs2         Extra inputs for the MAG processing part of the Snakemake workflow
-xb     --extra_inputs_bakta    Extra inputs (comma-separated) for MAG annotation using Bakta
```

### Using differential coverage binning
To perform genome recovery with differential coverage, prepare a 2-column comma-separated dataframe, indicating the additional read datatype (`NP` for Nanopore, `PB` for PacBio, `IL` for short reads) and read file location.<br/>
Dataframe example:
```
PB,/path/to/your/reads/file1.fastq
NP,/path/to/your/reads/file2.fastq
IL,/path/to/your/reads/file3.fastq.gz
```
The prepared dataframe can be provided to the workflow through the `-cov` option.

### Overview of main output (results directory)
* `<output_name>_assembly.fasta` - metagenome assembly file
* `<output_name>_16S.fa` - 16S rRNA gene sequences, recovered from the assembly
* `<output_name>_bins.tsv` - per-bin results [dataframe](msc/mmlong2-dfs.md#column-names-for-output_name_binstsv)
* `<output_name>_contigs.tsv` - per-contig results [dataframe](msc/mmlong2-dfs.md#column-names-for-output_name_contigstsv)
* `<output_name>_general.tsv` - workflow result summary as a single row [dataframe](msc/mmlong2-dfs.md#column-names-for-output_name_generaltsv)
* `dependencies.csv`- list of dependencies used and their versions
* `bins` - directory for metagenome assembled genomes
* `bakta` - directory, containing genome annotation results from [bakta](https://github.com/oschwengers/bakta)

## Additional documentation
* [Frequently asked questions](msc/mmlong2-faq.md)
* [Result dataframe documentation](msc/mmlong2-dfs.md)
* [Dependency list](msc/mmlong2-dep.md) and the specific [versions](arc)
* [Manual database setup](msc/mmlong2-db.md)

## Future improvements
Suggestions on improving the workflow or fixing bugs are always welcome.<br/>
Please use the GitHub `Issues` section or e-mail to mase@bio.aau.dk for providing feedback.

## Citation
If you use mmlong2 in a publication, please cite:
> Sereika, M., Mussig, A.J., Jiang, C. et al. Genome-resolved long-read sequencing expands known microbial diversity across terrestrial habitats. Nat Microbiol (2025). [https://doi.org/10.1038/s41564-025-02062-z](https://doi.org/10.1038/s41564-025-02062-z)

[//]: # (Written by Mantas Sereika)
