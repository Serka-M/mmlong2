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
| Completeness | Genome bin completeness estimate, reported by CheckM2 |
| Contamination | Genome bin contamination estimate, reported by CheckM2 |
| Completeness_Model_Used | Model used by CheckM2 |
| Coding_Density | Genome bin gene coding density, reported by CheckM2 |
| Contig_N50 | Genome bin N50 (bp), reported by CheckM2 |
| Average_Gene_Length | Average gene length, reported by CheckM2 |
| Genome_Size | Genome bin size (bp), reported by CheckM2 |
| GC_Content | Genome bin guanine-cytosine content (fraction), reported by CheckM2 |
| Total_Coding_Sequences | Number of coding sequences, reported by CheckM2 |
| contigs | Number of contigs, reported by Quast |
| longest_contig | Length of the longest contig (bp), reported by Quast |
| N90 | Genome bin N90 (bp), reported by Quast |
| [auN](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) | Nx area under the curve, reported by Quast |
| N_per_100kb | Rate of Ns in a genome bin per 100 kb, reported by Quast |
| cov | Average genome bin coverage, reported by CoverM |
| r_abund | Average genome bin relative abundance (%), reported by CoverM  |
| bakta_CDS_all | Number of protein coding genes, reported by Bakta |
| bakta_CDS_hyp | Number of hypothetical protein coding genes, reported by Bakta |
| bakta_tRNA_all | Number of all tRNA genes genes, reported by Bakta |
| bakta_tRNA_uniq | Number of unique tRNA genes genes, reported by Bakta  |
| bakta_16S | Number of 16S rRNA genes genes, reported by Bakta |
| bakta_23S | Number of 23S rRNA genes genes, reported by Bakta |
| bakta_5S | Number of 5S rRNA genes genes, reported by Bakta |
| barrnap_16S | Number of 16S rRNA genes genes, reported by Barrnap |
| barrnap_23S | Number of 23S rRNA genes genes, reported by Barrnap |
| barrnap_5S | Number of 5S rRNA genes genes, reported by Barrnap |
| MAG_status | Genome bin quality ranking according to [MIMAG standards](https://www.nature.com/articles/nbt.3893) |
| cMAG_status | MAG circularity status |
| gunc_contamination | Fraction of single-copy-genes from other clades, reported by GUNC |
| gunc_status | Chimeric MAG status, reported by GUNC |
| tax_gtdb| MAG taxonomy, reported by GTDB-tk |
| fastani_ani | Average nucleotide identity for a close match, reported by GTDB-tk |
| fastani_af | Alignment fraction for close match, reported by GTDB-tk |
| closest_placement_reference | ID for reference genome matched, reported by GTDB-tk |
| closest_placement_ani | Top ANI match thats is < 95 %, reported by GTDB-tk  |
| closest_placement_af | Alignment fraction for matches with < 95 % ANI, reported by GTDB-tk  |
| msa_percent | Percentage of amino acids in the multi-sequence alignment, reported by GTDB-tk |
| red_value | Relative Evolutionary Divergence for a distant match, reported by GTDB-tk |
| gtdb_warning | Warning message, reported by GTB-tk |
| tax_silva | MAG taxonomy based on 16S rRNA classification to [Silva](https://www.arb-silva.de/) database |
| tophit_silva | Percent identity for best match to the reference sequence, reported by vsearch |
| tax_midas | MAG taxonomy based on  16S rRNA classification to [MiDAS](https://www.midasfieldguide.org/guide) database |
| tophit_midas | Percent identity for best match to the reference sequence, reported by vsearch  |
| var_n | Number of detected variants in the MAG |
| var_perc | Fraction of MAG length as variants |
| wf_name | Workflow output name |
| wf_mode | Workflow mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

[//]: # (Written by Mantas Sereika)
