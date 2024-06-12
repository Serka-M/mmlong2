## Description of mmlong2 custom dataframes

### <output_name>_general.tsv

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

### <output_name>_contigs.tsv

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

### <output_name>_bins.tsv

| Category | Description |
| --- | --- |
| bin | Genome bin ID |
| completeness_checkm1 | Genome bin completeness estimate, reported by CheckM |
| contamination_checkm1 | Genome bin contamination estimate, reported by CheckM |
| strain_heterogeneity_checkm1 | Genome bin strain heterogeneity index, reported by CheckM |
| completeness_checkm2 | Genome bin completeness estimate, reported by CheckM2 |
| contamination_checkm2 | Genome bin contamination estimate, reported by CheckM2 |
| coding_density | Genome bin gene coding density, reported by CheckM2 |
| genome_size | Genome bin size (bp), reported by CheckM2 |
| contigs | Number of contigs, reported by CheckM2 |
| gc |  Genome bin guanine-cytosine content (%), reported by Quast |
| contig_n50 | Genome bin N50 (bp), reported by Quast |
| contig_n90 | Genome bin N90 (bp), reported by Quast |
| [aun](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) | Nx area under the curve, reported by Quast |
| n_per_100kb | Rate of Ns in a genome bin per 100 kb, reported by Quast |
| cov | Genome bin coverage, reported by CoverM |
| r_abund | Genome bin relative abundance (%), reported by CoverM  |
| bakta_cds_all | Number of protein coding genes, reported by Bakta |
| bakta_cds_hyp | Number of hypothetical protein coding genes, reported by Bakta |
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
| gunc_css | Clade separation [score](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02393-0), reported by GUNC |
| gunc_rrs | Reference representation [score](https://grp-bork.embl-community.io/gunc/output.html#output-columns), reported by GUNC |
| gunc_pass | Chimerism test status, reported by GUNC |
| var_n | Number of detected nucleotide variants, reported by Longshot |
| var_perc | Percentage of genome bin with nucleotide variants |
| tax_gtdb| Genome bin taxonomy, reported by GTDB-tk |
| ref_gtdb | Reference for the closest genome bin match, reported by GTDB-tk |
| ani_gtdb | Average nucleotide identity for the closest match, reported by GTDB-tk |
| af_gtdb | Alignment fraction for the closest match, reported by GTDB-tk |
| msa_gtdb | Percentage of the multi-sequence alignment spanned by the genome bin, reported by GTDB-tk |
| red_gtdb | Genome bin relative evolutionary divergence, reported by GTDB-tk |
| warning_gtdb | Warning message, reported by GTDB-tk |
| tax_rrna | Genome bin taxonomy, based on 16S rRNA classification to a rRNA database |
| tophit_rrna | Percent identity for the best match to a reference sequence in a rRNA database, reported by Vsearch |
| wf_name | Workflow output name |
| wf_mode | Workflow mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

[//]: # (Written by Mantas Sereika)
