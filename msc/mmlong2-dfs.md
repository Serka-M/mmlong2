## Description of mmlong2 custom dataframes

### [output_name]_general.tsv

| Category | Description |
| --- | --- |
| reads_n | Number of reads, reported by Nanoq |
| reads_size_bp | Read yield, reported by Nanoq |
| reads_N50_bp | Read N50, reported by Nanoq |
| reads_len_max_bp | Length of the longest read, reported by Nanoq |
| reads_len_min_bp | Length of the shortest read, reported by Nanoq |
| reads_mean_len_bp | Mean read length, reported by Nanoq |
| reads_median_len_bp | Median read length, reported by Nanoq |
| reads_mean_q | Mean read Phred quality score, reported by Nanoq |
| reads_median_q | Median read Phred quality score, reported by Nanoq |
| assembly_size_bp | Metagenome assembly size, reported by Flye |
| assembly_contigs_n | Number of contigs, reported by Flye |
| assembly_contigs_N50_bp | Contig N50, reported by Flye |
| assembly_contigs_meanCOV | Mean contig coverage, reported by Flye |
| assembly_reads_mapped | Fraction of read data that maps to the assembly, reported by Flye |
| contigs_circ | Number of circular contigs |
| all_mags | Total number of genome bins |
| circ_mags | Number of circular genome bins |
| hq_mags | Number of high-quality genome bins |
| mq_mags | Number of medium-quality genome bins |
| lq_mags | Number of low-quality genome bins |
| contaminated_mags | Number of contaminated genome bins |
| asm_binned | Percent of assembly size that was binned |
| r_abund_bins | Percent of read relative abundance explained by the bins |
| mag_cov_med_hq_mq | Median coverage of high-quality and medium-quality bins |
| mag_cov_med_hq | Median coverage of high-quality bins |
| mag_cov_med_mq | Median coverage of medium-quality binss |
| mag_cov_mad_hq_mq | Median absolute deviation in coverage of high-quality and medium-quality bins |
| mag_cov_mad_hq | Median absolute deviation in coverage of high-quality bins |
| mag_cov_mad_mq | Median absolute deviation in coverage of medium-quality binss |
| contigs_med_hq_mq | Median no. of contigs in high-quality and medium-quality bins |
| contigs_med_hq | Median no. of contigs in high-quality bins |
| contigs_med_mq | Median no. of contigs in medium-quality binss |
| wf_name | Workflow output name |
| wf_mode | Workflow mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

<br/>

### [output_name]_contigs.tsv

| Category | Description |
| --- | --- |
| contig | Contig ID |
| len_bp | Contig length |
| cov | Contig coverage by long reads |
| status_circular | Contig circularity status, reported by Flye |
| status_repeat | Contig repeat structure status, reported by Flye |
| graph_path | Links between contigs and edges of the assembly graph, reported by Flye |
| tax_kaiju | Contig taxonomy, reported by Kaiju |
| tax_silva | Contig taxonomy, based on 16S rRNA classification to [Silva](https://www.arb-silva.de/) database |
| tophit_silva | Percent identity for best match to the reference sequence, reported by vsearch |
| tax_midas | Contig taxonomy, based on  16S rRNA classification to [MiDAS](https://www.midasfieldguide.org/guide) database |
| tophit_midas | Percent identity for best match to the reference sequence, reported by vsearch |
| bin | ID for genome bin that contains the contig |
| GC | Contig guanine-cytosine percent content |
| var_n | Number of detected variants in the contig |
| var_perc | Fraction of contig length as variants |
| wf_name | Workflow output name |
| wf_mode | Workflow mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

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
| Total_Contigs | Number of contigs, reported by CheckM2 |
| Max_Contig_Length | Length of the longest contig (bp), reported by CheckM2 |
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
| tax_silva | MAG taxonomy, based on 16S rRNA classification to [Silva](https://www.arb-silva.de/) database |
| tophit_silva | Percent identity for best match to the reference sequence, reported by vsearch |
| tax_midas | MAG taxonomy, based on  16S rRNA classification to [MiDAS](https://www.midasfieldguide.org/guide) database |
| tophit_midas | Percent identity for best match to the reference sequence, reported by vsearch |
| var_n | Number of detected variants in the MAG |
| var_perc | Fraction of MAG length as variants |
| wf_name | Workflow output name |
| wf_mode | Workflow mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

[//]: # (Written by Mantas Sereika)
