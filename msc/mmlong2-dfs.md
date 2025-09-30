## Description of mmlong2 result dataframes

### Column names for <output_name>_general.tsv

| Category | Description |
| --- | --- |
| reads_n | Number of reads, reported by Nanoq |
| reads_yield_bp | Read yield, reported by Nanoq |
| reads_n50_bp | Read length N50 value, reported by Nanoq |
| reads_len_max_bp | Length of the longest read, reported by Nanoq |
| reads_len_min_bp | Length of the shortest read, reported by Nanoq |
| reads_mean_len_bp | Mean read length, reported by Nanoq |
| reads_median_len_bp | Median read length, reported by Nanoq |
| reads_mean_q | Mean read Phred quality score, reported by Nanoq |
| reads_median_q | Median read Phred quality score, reported by Nanoq |
| contigs_n | Number of assembled contigs, reported by Nanoq |
| contigs_yield_bp |  Metagenome assembly size, reported by Nanoq |
| contigs_n50_bp | Contig length N50 value, reported by Nanoq |
| contigs_len_max_bp | Length of the longest contig, reported by Nanoq |
| contigs_len_min_bp | Length of the shortest contig, reported by Nanoq |
| contigs_mean_len_bp | Mean contig length, reported by Nanoq |
| contigs_median_len_bp | Median contig length, reported by Nanoq |
| map_ident_median | Mean read identity (%) to the assembled metagenome, reported by Cramino |
| map_ident_mean | Median read identity (%) to the assembled metagenome, reported by Cramino |
| contigs_circ | Number of circular contigs, reported by the assembler |
| bins_all | Total number of genome bins |
| bins_circ | Number of circular genome bins |
| bins_hq | Number of high-quality genome bins |
| bins_mq | Number of medium-quality genome bins |
| bins_lq | Number of low-quality genome bins |
| bins_cont | Number of contaminated genome bins |
| bin_cov_median | Median genome bin coverage |
| yield_assembled | Percentage of sequenced data assembled into contigs |
| yield_binned | Percentage of sequenced data that was binned |
| assembly_binned | Percent of assembly size that was binned |
| wf_name | Workflow output name |
| wf_read_mode | Workflow long-read mode |
| wf_binning_mode | Workflow long-read mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

<br/>

### Column names for <output_name>_contigs.tsv

| Category | Description |
| --- | --- |
| contig | Contig identifier |
| len_bp | Contig length |
| cov | Contig coverage, reported by the assembler |
| status_circular | Contig circularity status, reported by the assembler |
| metabuli_tax | Contig taxonomy, reported by Metabuli |
| rrna_tax | Contig taxonomy, based on 16S rRNA classification to the rRNA reference database |
| rrna_tophit | Percent identity of the top hit in the rRNA reference database, reported by Usearch |
| rrna_alnlen | Alignment length (bp) of the top hit in the rRNA reference database, reported by Usearch |
| bin | Identifier for genome bin that contains the contig |
| gc | Contig guanine-cytosine content (%) |
| var_n | Number of detected nucleotide variants, reported by Longshot |
| var_perc | Percentage of contig with nucleotide variants |
| wf_name | Workflow output name |
| wf_read_mode | Workflow long-read mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

<br/>

### Column names for <output_name>_bins.tsv

| Category | Description |
| --- | --- |
| bin | Genome bin (MAG) identifier |
| completeness_checkm1 | Genome bin completeness estimate, reported by CheckM |
| contamination_checkm1 | Genome bin contamination estimate, reported by CheckM |
| strain_heterogeneity_checkm1 | Genome bin strain heterogeneity index, reported by CheckM |
| completeness_checkm2 | Genome bin completeness estimate, reported by CheckM2 |
| contamination_checkm2 | Genome bin contamination estimate, reported by CheckM2 |
| contigs | Number of contigs, reported by Quast |
| longest_contig | Longest contig length (bp), reported by Quast |
| genome_size | Genome bin size (bp), reported by Quast |
| gc | Genome bin guanine-cytosine content (%), reported by Quast |
| contig_n50 | Genome bin N50 (bp), reported by Quast |
| contig_n90 | Genome bin N90 (bp), reported by Quast |
| [aun](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) | Nx area under the curve, reported by Quast |
| n_per_100kb | Rate of Ns in a genome bin per 100 kbp, reported by Quast |
| cov | Genome bin coverage, reported by CoverM |
| r_abund | Genome bin relative abundance (%), reported by CoverM |
| bakta_cds_all | Number of protein coding genes, reported by Bakta |
| bakta_cds_hyp | Number of hypothetical protein coding genes, reported by Bakta |
| bakta_cds_dens | Genome coding density (%), reported by Bakta |
| bakta_trna_all | Number of all tRNA genes genes, reported by Bakta |
| bakta_trna_uniq | Number of unique tRNA genes genes, reported by Bakta  |
| bakta_16S | Number of 16S rRNA genes genes, reported by Bakta |
| bakta_23S | Number of 23S rRNA genes genes, reported by Bakta |
| bakta_5S | Number of 5S rRNA genes genes, reported by Bakta |
| barrnap_16S | Number of 16S rRNA genes genes, reported by Barrnap |
| barrnap_23S | Number of 23S rRNA genes genes, reported by Barrnap |
| barrnap_5S | Number of 5S rRNA genes genes, reported by Barrnap |
| custom_trna_uniq | Number of unique tRNA genes genes, reported by tRNAscan-SE |
| bin_status | Genome bin quality ranking according to [MIMAG standards](https://www.nature.com/articles/nbt.3893) |
| cbin_status | Genome bin circularity status |
| gunc_css | [Clade separation score](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02393-0), reported by GUNC |
| gunc_rrs | [Reference representation score](https://grp-bork.embl-community.io/gunc/output.html#output-columns), reported by GUNC |
| gunc_pass | Chimerism test status, reported by GUNC |
| var_n | Number of detected nucleotide variants, reported by Longshot |
| var_perc | Percentage of genome bin with nucleotide variants |
| gtdb_tax | Genome bin taxonomy, reported by GTDB-tk |
| gtdb_ref | Identifier for the closest genome reference match, reported by GTDB-tk |
| gtdb_ani | Average nucleotide identity for the closest match, reported by GTDB-tk |
| gtdb_af | Alignment fraction for the closest match, reported by GTDB-tk |
| gtdb_msa | Percentage of the multi-sequence alignment spanned by the genome bin, reported by GTDB-tk |
| gtdb_red | Genome bin relative evolutionary divergence, reported by GTDB-tk |
| gtdb_warning | Warning message, reported by GTDB-tk |
| rrna_tax | Genome bin taxonomy, based on 16S rRNA classification to the rRNA reference database |
| rrna_tophit | Percent identity of the top hit in the rRNA reference database, reported by Usearch |
| rrna_alnlen | Alignment length (bp) of the top hit in the rRNA reference database, reported by Usearch |
| wf_name | Workflow output name |
| wf_read_mode | Workflow long-read mode |
| wf_binning_mode | Workflow long-read mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

<br/>

### Column names for <output_name>_usage.tsv

| Category | Description |
| --- | --- |
| stage | Workflow stage (sequentially sorted) |
| step | Workflow step (alphabetically sorted) |
| s | Running time in seconds |
| h:m:s | Running time in hours:minutes:seconds format |
| max_rss | Maximum Resident Set Size (MB) - peak amount of non-swapped physical memory used by the process |
| max_vms | Maximum Virtual Memory Size (MB) - total amount of virtual memory used by the process |
| max_uss | Maximum Unique Set Size (MB) - peak memory unique to the process |
| max_pss | Maximum Proportional Set Size (MB) - peak shared memory, divided equally among processes sharing it |
| io_in | The number of MB read (cumulative) |
| io_out | The number of MB written (cumulative) |
| mean_load | CPU usage over time, divided by the total running time |
| cpu_time | Total CPU time (user + system), in seconds |

<br/>

### Row names for column `stage` in <output_name>_usage.tsv
| stage | Description |
| --- | --- |
| assembly | Metagenomic assembly |
| polishing | Optional polishing of the assembled metagenome |
| curation | Assembly curation |
| filtering | Metagenome filtering based on contig length and domain-level classification |
| singletons | Recovery of single-contig genomes |
| coverage | Read mapping and contig coverage calculation |
| binning | Recovery of multi-contig metagenomic bins |
| summary | Summary statistics for the recovered genomes |
| annotation | Annotation of the recovered genomes |
| taxonomy | Taxonomic classification (bin, contig, rRNA) |
| extraqc | Additional bin QC (chimerism and microdiversity) |
| stats | Miscellaneous workflow statistics |

<br/>

### Row names for column `step` in <output_name>_usage.tsv
| stage | step | Description |
| --- | --- | --- |
| assembly | myloasm | Metagenomic assembly with myloasm |
| assembly | metaflye | Metagenomic assembly with metaFlye |
| assembly | metamdbg | Metagenomic assembly with metaMDBG |
| assembly | custom | Import of custom assembly |
| polishing | prep | Preparation for polishing metagenome with Medaka |
| polishing | consensus | Generation of polishing consensus (per-batch) |
| polishing | stitch | Generation of the full consensus sequence |
| curation | map | Read multi-mapping to the assembly |
| curation | screening | Misassembly screening |
| curation | aggregate | Selection of contigs to curate |
| curation | selection | Generation of curated assembly |
| filtering | length| Metagenome filtering based on contig length |
| filtering | tiara | Detection of eukaryotic contigs with Tiara |
| filtering | whokaryote | Detection of eukaryotic contigs with Whokaryote |
| filtering | eukaryotes | Separation of eukaryotic contigs from the metagenome |
| singletons | circ | Extraction of circular contigs as separate genomes |
| singletons | lin | Extraction of linear contigs as separate genomes |
| singletons | qc | Quality control of single-contig genomes |
| coverage | prep | Preparation for read mapping to metagenome |
| coverage | map | Read mapping and contig coverage calculation (per-dataset) |
| coverage | aggregate | Aggregation of contig coverage profiles |
| binning | prep | Preparation for contig binning (per-round) |
| binning | metabat2 | Metagenomic binning with MetaBAT2 (per-round) |
| binning | semibin2 | Metagenomic binning with SemiBin2 (per-round) |
| binning | vamb | Metagenomic binning with VAMB (per-round) |
| binning | binette | Ensemble bin selection with Binette (per-round) |
| binning | prep-comebin-innit | Preparation for contig binning with COMEBin (per-round) |
| binning | prep-comebin | Read mapping preparation for contig binning with COMEBin (per-dataset, per-round) |
| binning | comebin | Metagenomic binning with COMEBin (per-round) |
| binning | qc | Quality control of the recovered bins (per-round) |
| binning | aggregate | Aggregation of the recovered bins |
| binning | qc2 | Secondary quality control of the recovered bins |
| summary | cov | Genome coverage and relative abundance calculation with CoverM |
| summary | stats | Genome summary statistic calculation with QUAST |

[//]: # (Written by Mantas Sereika)
