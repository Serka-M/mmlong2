## Description of mmlong2 custom dataframes

### [output_name]_general.tsv

| Category | Description |
| --- | --- |
| assembly_size_bp | Metagenome assembly size, reported by Flye |
| assembly_contigs_n | No. of contigs, reported by Flye |
| assembly_contigs_N50_bp | Contig N50, reported by Flye |
| assembly_contigs_meanCOV | Mean contig coverage, reported by Flye |
| assembly_reads_mapped | Fraction of read data that maps to the assembly, reported by Flye |
| num_reads_ilm | No. of short reads |
| bp_total_ilm | Short read yield |
| N50_ilm | Short read N50 |
| longest_reads_ilm | Length for the longest short read |
| shortest_reads_ilm | Length for the shortest short read |
| mean_length_ilm | Mean short read length |
| median_length_ilm | Median short read length |
| mean_q_ilm | Mean short read Phred quality score |
| median_q_ilm | Median short read Phred quality score |
| num_reads_long | No. of long reads |
| bp_total_long | Long read yield |
| N50_long | Long read N50 |
| longest_reads_long | Length for the longest long read |
| shortest_read_long | Length for the shortest long read |
| mean_length_long | Mean long read length |
| median_length_long | Median long read length |
| mean_q_long | Mean long read Phred quality score |
| median_q_long | Median long read Phred quality score |
| mags_workflow_name | Experiment name (same as output directory) |
| mags_workflow_date | Timestamp for dataframe generation |
| mags_workflow_mode | Workflow mode (PacBio/Nanopore/Nanopore-Illumina) |
| contigs_circ | No. of circular contigs |
| contigs_circ_above_HalfMb | No. of circular contigs that are above 0.5 Mb |
| contigs_abund_bac | Relative bacterial abundance, based on Kaiju classification |
| contigs_abund_arc | Relative archaeal abundance, based on Kaiju classification |
| all_mags | Total no. of genome bins |
| circ_mags | No. of circular genome bins |
| hq_mags | No. of high-quality genome bins |
| mq_mags | No. of medium-quality genome bins |
| lq_mags | No. of low-quality genome bins |
| contaminated_mags | No. of contaminated genome bins |
| asm_binned | Percent of assembly size that was binned |
| Rabund_binned_long | Percent long read abundance explained by bins |
| Rabund_binned_ilm | Percent short read abundance explained by bins |
| mag_cov_med_hq_mq | Median coverage of high-quality and medium-quality bins |
| mag_cov_med_hq | Median coverage of high-quality bins |
| mag_cov_med_mq | Median coverage of medium-quality binss |
| mag_cov_mad_hq_mq | Median absolute deviation in coverage of high-quality and medium-quality bins |
| mag_cov_mad_hq | Median absolute deviation in coverage of high-quality bins |
| mag_cov_mad_mq | Median absolute deviation in coverage of medium-quality binss |
| contigs_med_hq_mq | Median no. of contigs in high-quality and medium-quality bins |
| contigs_med_hq | Median no. of contigs in high-quality bins |
| contigs_med_mq | Median no. of contigs in medium-quality binss |

<br/>

### [output_name]_contigs.tsv

| Category | Description |
| --- | --- |
| contig | Contig ID |
| contig_len_bp | Contig length |
| cov_long | Contig coverage by long reads |
| cov_long_var | Variance in contig coverage by long reads |
| cov_ilm | Contig coverage by short reads|
| cov_ilm_var | Variance in contig coverage by short reads |
| GC | Contig guanine-cytosine percent content |
| status_circular | Contig circularity status |
| status_repeat | Contig repeat structure status |
| graph_path | Links between contig and edges of the assembly graph |
| medak_var_n | No. of variants determined by Medaka |
| medak_var_perc | Percent rate of variants determined by Medaka |
| domain | Contig domain taxonomy determined by Kaiju |
| phylum | Contig phylum taxonomy determined by Kaiju |
| class | Contig class taxonomy determined by Kaiju |
| order | Contig order taxonomy determined by Kaiju |
| family | Contig family taxonomy determined by Kaiju |
| species | Contig species taxonomy determined by Kaiju (not very reliable) |
| mags_workflow_name | Experiment name (same as output directory) |
| mags_workflow_date | Timestamp for dataframe generation |
| mags_workflow_mode | Workflow mode (PacBio/Nanopore/Nanopore-Illumina) |
| bin | ID for genome bin that contains the contig |

<br/>

### [output_name]_bins.tsv

| Category | Description |
| --- | --- |
| bin | Genome bin ID |
| CDS_all | No. of protein coding genes reported by Bakta |
| CDS_hyp | No. of hypothetical protein coding genes reported by Bakta |
| tRNA_all | No. of all tRNA genes genes reported by Bakta |
| tRNA_all | No. of unique tRNA genes genes reported by Bakta |
| rRNA_16S | No. of 16S rRNA genes genes reported by Bakta |
| rRNA_23S | No. of 23S rRNA genes genes reported by Bakta |
| rRNA_5S | No. of 5S rRNA genes genes reported by Bakta |
| contigs | No. of contigs |
| Longest_contig_bp | Longest contig length |
| Genome_size_bp | Genome bin size |
| GC | Genome bin guanine-cytosine content |
| N50_bp | Contig N50 |
| [auN](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) | Nx area under the curve |
| N_per_100kb | Rate of Ns in a genome bin |
| Coding_density | Genome bin coding density estimate by CheckM2 |
| Completeness | Genome bin completeness estimate by CheckM2 |
| Contamination | Genome bin contamination estimate by CheckM2 |
| gunc_contamination | Fraction of single-copy-genes from other clades by GUNC |
| gunc_status | Chimeric MAG status by GUNC |
| gtdb_classification | Bin taxonomy by GTDB-tk (specific) |
| fastani_ani | Average nucleotide identity for a close match by GTDB-tk |
| fastani_af | Alignment fraction for close match by GTDB-tk |
| closest_placement_reference | ID for reference genome matched by GTDB-tk |
| closest_placement_ani | Top ANI match thats is < 95 % by GTDB-tk  |
| closest_placement_af | Alignment fraction for matches with < 95 % ANI by GTDB-tk  |
| msa_percent | Percentage of amino acids in the multi-sequence alignment |
| red_value | Relative Evolutionary Divergence for a distant match by GTDB-tk |
| gtdb_warning | Warning message by GTB-tk |
| silva_classification | 16S rRNA classification to [Silva](https://www.arb-silva.de/) database |
| silva_identity | Percent identity to the reference sequence |
| midas_classification | 16S rRNA classification to [MiDAS](https://www.midasfieldguide.org/guide) database |
| midas_identity | Percent identity to the reference sequence |
| mags_workflow_name | Experiment name (same as output directory) |
| mags_workflow_date | Timestamp for dataframe generation |
| mags_workflow_mode | Workflow mode (PacBio/Nanopore/Nanopore-Illumina) |
| MAG_status | Genome bin quality ranking according to [MIMAG standards](https://www.nature.com/articles/nbt.3893) |
| cov_long | Bin coverage by long reads |
| cov_ilm | Bin coverage by short reads |
| Rabund_long | Bin relative abundance in long read data |
| Rabund_ilm | Bin relative abundance in short read data  |
| medak_var_n | No. of variants determined by Medaka |
| medak_var_perc | Percent rate of variants determined by Medaka |
| cMAG_status | No. of circular Genome bins |
| kaiju_domain | Genome bin domain taxonomy by Kaiju |
| kaiju_phylum | Genome bin phylum taxonomy by Kaiju |

[//]: # (Written by Mantas Sereika)
