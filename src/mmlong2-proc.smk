# DESCRIPTION: Snakemake workflow for analysing metagenome assembled genomes
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

import os
import re
import sys
import glob
import pandas as pd
import numpy as np

shell.executable("/bin/bash")
singularity: config["sing"]

wf_v=config["version"]
mode=config["mode"]
modes=['Nanopore-simplex', 'PacBio-HiFi']

dbs_barrnap=['bac', 'arc', 'euk', 'mito']
dbs_trnascan=['G', 'B', 'A']

proc=config["proc"]
proc_sub=config["proc_sub"]

loc=config["loc"]
sample=config["sample"]
fastq=config["fastq"]

def get_bins():	
    bins_full=glob.glob(os.path.join(loc, sample, "results/bins/*.fa"))
    bins = [os.path.splitext(os.path.basename(file))[0] for file in bins_full]
    return bins

def get_splits(proc,proc_sub,split_max):
    if int(proc/proc_sub) < split_max: return int(proc/proc_sub)
    else: return split_max

def get_mapping(mode):
    if mode == "PacBio-HiFi": return "PB"
    else: return "NP"

def get_bakta_params(bakta_extra):
    if bakta_extra == "FALSE": return " "
    else: return re.sub(","," ",bakta_extra)

onstart:
    from snakemake.utils import min_version
    min_version("9.9.0")
    if not os.path.exists(config["db_metabuli"]): sys.exit(print("Metabuli database input (((",config["db_metabuli"],"))) not found. Aborting..."))
    if not os.path.exists(config["db_bakta"]): sys.exit(print("Bakta database input (((",config["db_bakta"],"))) not found. Aborting..."))
    if not os.path.exists(config["db_gtdb"]): sys.exit(print("GTDB-tk database input (((",config["db_gtdb"],"))) not found. Aborting..."))
    if not os.path.exists(config["db_rrna"]): sys.exit(print("16S rRNA database input (((",config["db_rrna"],"))) not found. Aborting..."))
    if not os.path.exists(config["db_gunc"]): sys.exit(print("GUNC database input (((",config["db_gunc"],"))) not found. Aborting..."))
    if not os.path.exists(fastq): sys.exit(print("Read input (((",fastq,"))) not found. Aborting..."))
    if not config["db_barrnap"] in dbs_barrnap: sys.exit(print("Provided Barrnap database (((",db_barrnap,"))) not recognised. Aborting..."))
    if not config["db_trnascan"] in dbs_trnascan: sys.exit(print("Provided tRNAscan-SE database (((",db_trnascan,"))) not recognised. Aborting..."))
    if not mode in modes: sys.exit(print("Provided workflow mode (((",mode,"))) not recognised. Aborting..."))
    if len(os.path.join(loc, sample)) > 85: sys.exit(print("Path for provided output too long: (((",os.path.join(loc, sample),")))\nPlease re-run with different output location."))
    if not os.path.exists(os.path.join(loc, sample)): os.makedirs(os.path.join(loc, sample))
    if not os.path.exists(os.path.join(loc, sample, "results")): os.makedirs(os.path.join(loc, sample, "results"))
    if not os.path.exists(os.path.join(loc, sample, "tmp")): os.makedirs(os.path.join(loc, sample, "tmp"))
    if not os.path.exists(os.path.join(loc, sample, "tmp", "logs")): os.makedirs(os.path.join(loc, sample, "tmp", "logs"))
    if not os.path.exists(os.path.join(loc, sample, "tmp", "dep_mmlong2-proc.csv")): 
	    with open(os.path.join(loc, sample, "tmp", "dep_mmlong2-proc.csv"), 'w') as f:
		    f.write("dependency,version\n")
    if not os.path.exists(os.path.join(loc, sample, "tmp", "db_mmlong2-proc.csv")): 
	    with open(os.path.join(loc, sample, "tmp", "db_mmlong2-proc.csv"), 'w') as f:
		    f.write("database,path\n")

onsuccess:
    from datetime import datetime
    now = datetime.now()
    print("MAG processing with version",wf_v,"completed at",now.strftime("%Y/%m/%d %H:%M:%S"))
    print("Thank you for using mmlong2")

onerror:
    print("An error has occurred. Inspect Snakemake log files for troubleshooting.")

rule Finalise:
    conda: config["env_10"]
    input:
        gen_stats=os.path.join(loc, sample, "tmp/stats/gen_stats.tsv"),
        contigs_stats=os.path.join(loc, sample, "tmp/stats/contigs_stats.tsv"),
        contigs_extraqc=os.path.join(loc, sample, "tmp/extra_qc/contigs_extraqc.tsv"),
        contigs_taxonomy=os.path.join(loc, sample, "tmp/taxa/contigs_taxonomy.tsv"),
        bins=os.path.join(loc, sample, "tmp/binning/bins_mmlong2-lite.tsv"),
        bins_extraqc=os.path.join(loc, sample, "tmp/extra_qc/bins_extraqc.tsv"),
        bins_taxonomy=os.path.join(loc, sample, "tmp/taxa/bins_taxonomy.tsv"),
        bins_annotation=os.path.join(loc, sample, "tmp/annotation/bins_annotation.tsv"),
        usage=os.path.join(loc, sample, "tmp/logs/summary_mmlong2-proc.tsv"),
        clean=os.path.join(loc, sample, "tmp/logs/cleanup-proc.txt"),
    output:
        dep=os.path.join(loc, sample, "results/dependencies.csv"),
        db=os.path.join(loc, sample, "results/databases.csv"),
        gen=os.path.join(loc, sample, "results", f"{sample}_general.tsv"),
        contigs=os.path.join(loc, sample, "results", f"{sample}_contigs.tsv"),
        bins=os.path.join(loc, sample, "results", f"{sample}_bins.tsv"),
    shell:
        """
        cat {loc}/{sample}/tmp/dep_mmlong2-lite.csv > {output.dep}
        tail -n+2 {loc}/{sample}/tmp/dep_mmlong2-proc.csv >> {output.dep}
        cat {loc}/{sample}/tmp/db_mmlong2-proc.csv > {output.db}

        R --no-echo --silent --args << 'df' > /dev/null 2>&1
        # Load data
        bins <- read.delim("{input.bins}", sep="\t", header=T)
        bins_annot <- read.delim("{input.bins_annotation}", sep="\t", header=T)
        bins_tax <- read.delim("{input.bins_taxonomy}", sep="\t", header=T)
        bins_qc <- read.delim("{input.bins_extraqc}", sep="\t", header=T)
        
        contigs <- read.delim("{input.contigs_taxonomy}", sep="\t", header=T)
        contigs_stats <- read.delim("{input.contigs_stats}", sep="\t", header=T)
        contigs_qc <- read.delim("{input.contigs_extraqc}", sep="\t", header=T)
        
        gen <- read.delim("{input.gen_stats}", sep="\t", header=T)
        
        # Combine data
        contigs <- merge(contigs,contigs_stats,by="contig")
        contigs <- merge(contigs,contigs_qc,by="contig", all=TRUE)
        contigs[is.na(contigs$var_n),]$var_n <- 0
        contigs$var_perc <- round(contigs$var_n/contigs$len_bp*100,3)
        contigs$wf_name <- "{sample}"
        contigs$wf_read_mode <- "{mode}"
        contigs$wf_v <- "{wf_v}"
        contigs$wf_date <- Sys.Date()
        
        bins <- merge(bins,bins_annot, by="bin")
        bins$bin_status <- ifelse((bins$completeness_checkm2 >= 90 & bins$contamination_checkm2 <= 5 &
			(bins$bakta_trna_uniq >= 18 | bins$custom_trna_uniq) &
			(bins$bakta_5s >= 1 | bins$barrnap_5s >= 1) &
			(bins$bakta_16s >= 1 | bins$barrnap_16s >= 1) &
			(bins$bakta_23s >= 1 | bins$barrnap_23s >= 1)),"HQ",
			ifelse(bins$completeness_checkm2 >= 50 & bins$contamination_checkm2 <= 10, "MQ",
			ifelse(bins$completeness_checkm2 <= 50 & bins$contamination_checkm2 <= 10, "LQ", "Contaminated")))
        
        bins$cbin_status <- ifelse(grepl("bin.c",bins$bin),"Y","N")
        bins <- merge(bins,bins_qc , by="bin")
        bins$var_perc <- round(bins$var_n/bins$genome_size*100,3)
        bins <- merge(bins,bins_tax, by="bin")
        
        wf_read_mode <- unique(bins$wf_read_mode)
        wf_binning_mode <- unique(bins$wf_binning_mode)
        bins$wf_name <- NULL
        bins$wf_v <- NULL
        bins$wf_read_mode <- NULL
        bins$wf_binning_mode <- NULL
        bins$wf_date <- NULL
        
        bins$wf_name <- "{sample}"
        bins$wf_read_mode <- wf_read_mode
        bins$wf_binning_mode <- wf_binning_mode
        bins$wf_v <- "{wf_v}"
        bins$wf_date <- Sys.Date()
        
        gen$contigs_circ <- nrow(contigs[(contigs$status_circular == "Y"), ])
        gen$bins_all <- nrow(bins)
        gen$bins_circ <- nrow(bins[(bins$cbin_status == "Y"), ])
        gen$bins_hq <- nrow(bins[(bins$bin_status == "HQ"), ])
        gen$bins_mq <- nrow(bins[(bins$bin_status == "MQ"), ])
        gen$bins_lq <- nrow(bins[(bins$bin_status == "LQ"), ])
        gen$bins_cont <- nrow(bins[(bins$bin_status == "Contaminated"), ])
        
        gen$bin_cov_median <- median(bins$cov)
        gen$yield_assembled <- round(gen$map_yield_gbp*10**9/gen$reads_yield_bp*100,3)
        gen$yield_binned <- round(gen$yield_assembled*sum(bins$r_abund)/100,3)
        gen$assembly_binned <- round(sum(bins$genome_size)/gen$contigs_yield_bp*100,3)
        gen[,"map_yield_gbp"] <- NULL
        
        gen$wf_name <- "{sample}"
        gen$wf_read_mode <- wf_read_mode
        gen$wf_binning_mode <- wf_binning_mode
        gen$wf_v <- "{wf_v}"
        gen$wf_date <- Sys.Date()
        
        # Save data
        write.table(gen,"{output.gen}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        write.table(contigs,"{output.contigs}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        write.table(bins,"{output.bins}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Stats_aggregate:
    conda: config["env_10"]
    input:
        reads=os.path.join(loc, sample, "tmp/stats/reads.txt"),
        contigs=os.path.join(loc, sample, "tmp/stats/contigs.txt"),
        gc=os.path.join(loc, sample, "tmp/stats/gc.tsv"),
        map=os.path.join(loc, sample, "tmp/stats/map_filt.tsv"),
    output:
        contigs=os.path.join(loc, sample, "tmp/stats/contigs_stats.tsv"),
        gen=os.path.join(loc, sample, "tmp/stats/gen_stats.tsv")
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_stats_aggregate.tsv")
    shell:
        """
        R --no-echo --silent --args << 'df' > /dev/null 2>&1
        # Sample-level stats
        reads=read.delim("{input.reads}", sep=" ", header=F)
        colnames(reads) <- c("reads_n", "reads_yield_bp", "reads_n50_bp","reads_len_max_bp","reads_len_min_bp", "reads_mean_len_bp","reads_median_len_bp","reads_mean_q","reads_median_q")
        
        contigs_gen=read.delim("{input.contigs}", sep=" ", header=F)
        colnames(contigs_gen) <- c("contigs_n", "contigs_yield_bp", "contigs_n50_bp","contigs_len_max_bp","contigs_len_min_bp", "contigs_mean_len_bp","contigs_median_len_bp","del1","del2")   
        contigs_gen$del1 <- NULL
        contigs_gen$del2 <- NULL
        
        map=read.delim("{input.map}", sep="\t", header=T)
        colnames(map) <- c("mean_cov","map_yield_gbp", "map_ident_median", "map_ident_mean")
        
        gen <- cbind(reads,contigs_gen,map)
        write.table(gen,"{output.gen}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        
        # Contig-level stats
        contigs=read.delim("{input.gc}", sep="\t", header=T)
        colnames(contigs) <- c("contig", "gc")
        write.table(contigs,"{output.contigs}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Stats_seq:
    conda: config["env_13"]
    input:
        reads={fastq},
        contigs=os.path.join(loc, sample, "results", f"{sample}_assembly.fasta"),
    output:
        reads=os.path.join(loc, sample, "tmp/stats/reads.txt"),
        contigs=os.path.join(loc, sample, "tmp/stats/contigs.txt"),
        gc=os.path.join(loc, sample, "tmp/stats/gc.tsv")
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_stats_seq.tsv")
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/stats" ]; then mkdir {loc}/{sample}/tmp/stats; fi
        seqkit fx2tab -n --gc {loc}/{sample}/results/{sample}_assembly.fasta > {output.gc}
        nanoq -i {input.contigs} -s -t {threads} --report {output.contigs}
        nanoq -i {input.reads} -s -t {threads} --report {output.reads}
        if ! grep -q "nanoq" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "nanoq " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule Stats_map:
    conda: config["env_13"]
    input:
        expand("{loc}/{sample}/tmp/binning/mapping/1-{map}.bam",sample=sample,loc=loc,map=get_mapping(mode))
    output:
        innit=os.path.join(loc, sample, "tmp/stats/map.tsv"),
        filt=os.path.join(loc, sample, "tmp/stats/map_filt.tsv"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_stats_map.tsv")
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/stats" ]; then mkdir {loc}/{sample}/tmp/stats; fi
        cramino --threads {threads} {input} > {output.innit}
        awk 'BEGIN {{FS = "\t"}} {{col1[NR] = $1; col2[NR] = $2}} END {{print col1[6] "\t" col1[5] "\t" col1[13] "\t" col1[14]; print col2[6] "\t" col2[5] "\t" col2[13] "\t" col2[14];}}' {output.innit} > {output.filt}
        if ! grep -q "cramino" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "cramino " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule ExtraQC_aggregate:
    conda: config["env_10"]
    input:
        variants=os.path.join(loc, sample, "tmp/extra_qc/variants/variants.tsv"),
        links=os.path.join(loc, sample, "tmp/binning/contig_bin.tsv"),
        gunc=os.path.join(loc, sample, "tmp/extra_qc/gunc.tsv"),
    output:
        contigs=os.path.join(loc, sample, "tmp/extra_qc/contigs_extraqc.tsv"),
        bins=os.path.join(loc, sample, "tmp/extra_qc/bins_extraqc.tsv")
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_extraqc_aggregate.tsv")
    shell:
        """
        R --no-echo --silent --args << 'df' > /dev/null 2>&1
        # Contig-level stats
        variants <- read.delim("{input.variants}", sep="\t", header=F)
        colnames(variants) <- c("contig","var_n")
        write.table(variants,"{output.contigs}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        
        # Bin-level stats
        links <- read.delim("{input.links}", sep="\t", header=F)
        colnames(links) <- c("contig","bin")
        
        gunc <- read.delim("{input.gunc}", sep="\t", header=T)
        colnames(gunc) <- c("bin", "gunc_css", "gunc_rrs", "gunc_pass")
        
        variants <- merge(variants,links,by="contig",all=TRUE)
        variants[is.na(variants$var_n),]$var_n <- 0
        
        bins <- aggregate(variants$var_n, by=list(variants$bin), FUN=sum)
        colnames(bins) <- c("bin", "var_n")
        
        bins <- merge(gunc,bins,by="bin")
        write.table(bins,"{output.bins}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule ExtraQC_variants_prep:
    conda: config["env_14"]
    input:
        expand("{loc}/{sample}/results/bins/{bin}.fa",sample=sample,loc=loc,bin=get_bins())
    output:
        contigs=os.path.join(loc, sample, "tmp/extra_qc/variants/contigs.fasta"),
        list=os.path.join(loc, sample, "tmp/extra_qc/variants/contigs.txt"),
        ids=expand("{loc}/{sample}/tmp/extra_qc/variants/contigs_{split}.txt",sample=sample,loc=loc,split=[f"{i:02}" for i in range(1, get_splits(proc,proc_sub,config["split_max"])+1)]),
    params:
        splits=get_splits(proc,proc_sub,config["split_max"]),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_extraqc_variants-prep.tsv")
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/extra_qc" ]; then mkdir {loc}/{sample}/tmp/extra_qc; fi
        if [ -d "{loc}/{sample}/tmp/extra_qc/variants" ]; then rm -r {loc}/{sample}/tmp/extra_qc/variants; fi
        mkdir {loc}/{sample}/tmp/extra_qc/variants
        
        cut -f1 {loc}/{sample}/tmp/binning/contig_bin.tsv > {output.list}
        split -n l/{params.splits} --numeric-suffixes=1 --additional-suffix=.txt -d {output.list} {loc}/{sample}/tmp/extra_qc/variants/contigs_
        
        cat  {loc}/{sample}/results/bins/*.fa > {output.contigs}
        samtools faidx {output.contigs}
        """

rule ExtraQC_variants_call:
    conda: config["env_14"]
    input:
        "{loc}/{sample}/tmp/extra_qc/variants/contigs_{split}.txt",
    output:
        map="{loc}/{sample}/tmp/extra_qc/variants/contigs_{split}.bam",
        var="{loc}/{sample}/tmp/extra_qc/variants/contigs_{split}.vcf",
    params:
        map=get_mapping(mode),
        min_cov=config["min_cov"],
        min_mapq=config["min_mapq"],
        min_frac=config["min_frac"],
        min_count=config["min_count"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_extraqc_variants-{split}.tsv"
    resources: usage=config["longshot_usage"]
    threads: proc_sub
    shell:
        """
        samtools view --threads $(({threads} / 2)) -bh {loc}/{sample}/tmp/binning/mapping/1-{params.map}.bam $(cat {input}) | samtools sort --threads $(({threads} / 2)) --write-index -o {output.map}##idx##{output.map}.bai
        longshot --bam {output.map} --ref {loc}/{sample}/tmp/extra_qc/variants/contigs.fasta --out {output.var} -n --min_cov {params.min_cov} --min_mapq {params.min_mapq} --min_alt_frac {params.min_frac} --min_alt_count {params.min_count}
        """

rule ExtraQC_variants_sum:
    conda: config["env_14"]
    input:
        expand("{loc}/{sample}/tmp/extra_qc/variants/contigs_{split}.vcf",sample=sample,loc=loc,split=[f"{i:02}" for i in range(1, get_splits(proc,proc_sub,config["split_max"])+1)]),
    output:
        var_sort=os.path.join(loc, sample, "tmp/extra_qc/variants/variants_sorted.vcf"),
        var_sum=os.path.join(loc, sample, "tmp/extra_qc/variants/variants.tsv"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_extraqc_variants-sum.tsv")
    threads: proc
    shell:
        """
        grep -hv '^#' {input} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' > {output.var_sort}
        awk -F "\t" '{{ if ($7 == "PASS") {{print $1}} }}' {output.var_sort} | uniq -c - | awk '{{$1=$1;print}}' | awk -F" " '{{print $2,$1}}' OFS='\t' > {output.var_sum}
        if ! grep -q "longshot" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "longshot " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule ExtraQC_contamination:
    conda: config["env_13"]
    input:
        expand("{loc}/{sample}/results/bins/{bin}.fa",sample=sample,loc=loc,bin=get_bins())
    output:
        raw=os.path.join(loc, sample, "tmp/extra_qc/gunc/gunc_raw.tsv"),
        trunc=os.path.join(loc, sample, "tmp/extra_qc/gunc.tsv")
    params:
        db=config["db_gunc"]
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_extraqc_gunc.tsv")
    resources: usage=100
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/extra_qc" ]; then mkdir {loc}/{sample}/tmp/extra_qc; fi
        if [ ! -d "{loc}/{sample}/tmp/extra_qc/tmp" ]; then mkdir {loc}/{sample}/tmp/extra_qc/tmp; fi
        if [ ! -d "{loc}/{sample}/tmp/extra_qc/gunc" ]; then rm -r {loc}/{sample}/tmp/extra_qc/gunc; fi

        gunc run --db_file {params.db} --input_dir {loc}/{sample}/results/bins --out_dir {loc}/{sample}/tmp/extra_qc/gunc --threads {threads} --temp_dir {loc}/{sample}/tmp/extra_qc/tmp
        mv {loc}/{sample}/tmp/extra_qc/gunc/*.tsv {output.raw}
        cut -f1,8,12,13 {output.raw} > {output.trunc}

        if grep -q "GUNC," {loc}/{sample}/tmp/db_mmlong2-proc.csv; then sed -i '/GUNC,/d' {loc}/{sample}/tmp/db_mmlong2-proc.csv; fi
        echo "GUNC,{params.db}" >> {loc}/{sample}/tmp/db_mmlong2-proc.csv
        if ! grep -q "gunc" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "gunc " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule Annotation_aggregate:
    conda: config["env_10"]
    input:
        bakta=os.path.join(loc, sample, "tmp/annotation/bakta_stats.csv"),
        barrnap=os.path.join(loc, sample, "tmp/annotation/barrnap_stats.csv"),
        trnascan=os.path.join(loc, sample, "tmp/annotation/trna_stats.csv")
    output:
        df=os.path.join(loc, sample, "tmp/annotation/bins_annotation.tsv")
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_annotation_aggregate.tsv")
    shell:
        """
        R --no-echo --silent --args << 'df' > /dev/null 2>&1
        bakta=read.delim("{input.bakta}", sep=",", header=T)
        barrnap=read.delim("{input.barrnap}", sep=",", header=T)
        trnascan=read.delim("{input.trnascan}", sep=",", header=T)
        annot=merge(merge(bakta,barrnap,by="bin"),trnascan,by="bin")
        write.table(annot,"{output.df}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Annotation_rrna:
    conda: config["env_10"]
    input:
        expand("{loc}/{sample}/results/bins/{bin}.fa",sample=sample,loc=loc,bin=get_bins())
    output:
        df=os.path.join(loc, sample, "tmp/annotation/barrnap_stats.csv")
    params:
        db=config["db_barrnap"]
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_annotation_barrnap.tsv")
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/annotation" ]; then mkdir {loc}/{sample}/tmp/annotation; fi
        if [ ! -d "{loc}/{sample}/tmp/annotation/tmp" ]; then mkdir {loc}/{sample}/tmp/annotation/tmp; fi
        if [ -d "{loc}/{sample}/tmp/annotation/trna" ]; then rm -r {loc}/{sample}/tmp/annotation/trna; fi
        if [ -d "{loc}/{sample}/tmp/annotation/bins" ]; then rm -r {loc}/{sample}/tmp/annotation/bins; fi
        mkdir {loc}/{sample}/tmp/annotation/trna
        mkdir {loc}/{sample}/tmp/annotation/bins
        
        function rrna_stats {{
        FILE=$1
        FILE_NAME=${{FILE##*/}}
        FILE_NAME=${{FILE_NAME%%.*}}
        BIN=$(basename $FILE | sed "s/.fa//")
        S=`barrnap --threads {threads} --kingdom {params.db} --quiet $FILE |\
        awk -F "=" -v bin=$BIN '
        NR == FNR {{a[$1]=0; next}}
        /^[^#]/{{gsub(/ .*/, "", $3); a[$3]++}}
        END{{OFS = ","; print bin, a["16S"], a["23S"], a["5S"]}}
        ' <(printf "%s\n" 16S 23S 5S) -` 
        echo "$S"; }}
        
        echo "bin, barrnap_16s, barrnap_23s, barrnap_5s" > {output.df}
        find {loc}/{sample}/results/bins -name "*.fa" | while read file; do rrna_stats "$file" >> {output.df}; done
        if ! grep -q "barrnap" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "barrnap " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule Annotation_trna:
    conda: config["env_12"]
    input:
        bin="{loc}/{sample}/results/bins/{bin}.fa",
        rrna="{loc}/{sample}/tmp/annotation/barrnap_stats.csv",
    output:
        "{loc}/{sample}/tmp/annotation/trna/trna_{bin}.txt"
    params:
        db=config["db_trnascan"],
        bin="{bin}"
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_annotation_trnascan-{bin}.tsv"
    resources: usage=config["trnascan_usage"]
    threads: proc_sub
    shell:
        """
        if [ -f "{loc}/{sample}/tmp/annotation/trna/trna_{params.bin}.txt" ]; then rm {loc}/{sample}/tmp/annotation/trna/trna_{params.bin}.txt; fi
        if [ -f "{loc}/{sample}/tmp/annotation/trna/stats_{params.bin}.txt" ]; then rm {loc}/{sample}/tmp/annotation/trna/stats_{params.bin}.txt; fi
        tRNAscan-SE -{params.db} -o {loc}/{sample}/tmp/annotation/trna/trna_{params.bin}.txt -m {loc}/{sample}/tmp/annotation/trna/stats_{params.bin}.txt -d {input.bin} --thread {threads}
        sed -i -e '1,3'd -e "s/$/\t{params.bin}/g" {loc}/{sample}/tmp/annotation/trna/trna_{params.bin}.txt
        """

rule Annotation_trna_sum:
    conda: config["env_12"]
    input:
        expand("{loc}/{sample}/tmp/annotation/trna/trna_{bin}.txt",sample=sample,loc=loc,bin=get_bins())
    output:
        df=os.path.join(loc, sample, "tmp/annotation/trna_stats.csv")
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_annotation_trnascan-sum.tsv")
    threads: proc
    shell:
        """
        cat {loc}/{sample}/tmp/annotation/trna/trna_*.txt | cut -f11,5 | sed 's/\t/,/g' | sed '/Undet/d' | sed '/Sup/d' > {loc}/{sample}/tmp/annotation/trna/stats.csv
        awk -F $',' ' {{ t = $1; $1 = $2; $2 = t; print; }} ' OFS=$',' {loc}/{sample}/tmp/annotation/trna/stats.csv > {loc}/{sample}/tmp/annotation/trna/stats_inv.csv
        echo "bin,custom_trna_uniq" > {output.df}
        sort -u < {loc}/{sample}/tmp/annotation/trna/stats_inv.csv | cut -d, -f1 | uniq -c | awk -vOFS=, '{{ print $2,$1 }}' >> {output.df}
        if ! grep -q "trnascan-se" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "trnascan-se " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi  
        """

rule Annotation_bins:
    conda: config["env_12"]
    input:
        bin="{loc}/{sample}/results/bins/{bin}.fa",
        rrna="{loc}/{sample}/tmp/annotation/barrnap_stats.csv",
    output:
        "{loc}/{sample}/tmp/annotation/bins/{bin}/{bin}.tsv"
    params:
        db=config["db_bakta"],
        extra=get_bakta_params(config["bakta_extra"]),
        bin="{bin}"
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_annotation_bakta-{bin}.tsv"
    resources: usage=config["bakta_usage"]
    threads: proc_sub
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/annotation/bins/{params.bin}" ]; then rm -r {loc}/{sample}/tmp/annotation/bins/{params.bin}; fi
        export TMPDIR={loc}/{sample}/tmp/annotation/tmp
        bakta --db {params.db} --prefix {params.bin} --output {loc}/{sample}/tmp/annotation/bins/{params.bin} --keep-contig-headers --tmp-dir {loc}/{sample}/tmp/annotation/tmp --threads {threads} {input.bin} --meta --force {params.extra}
        bakta_plot --prefix {params.bin} --output {loc}/{sample}/tmp/annotation/bins/{params.bin} --force {loc}/{sample}/tmp/annotation/bins/{params.bin}/{params.bin}.json
        """

rule Annotation_bins_sum:
    conda: config["env_12"]
    input:
        expand("{loc}/{sample}/tmp/annotation/bins/{bin}/{bin}.tsv",sample=sample,loc=loc,bin=get_bins())
    output:
        df=os.path.join(loc, sample, "tmp/annotation/bakta_stats.csv")
    params:
        db=config["db_bakta"]
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_annotation_bakta-sum.tsv")
    threads: proc
    shell:
        """
        echo "bin,bakta_cds_all,bakta_cds_hyp,bakta_cds_dens,bakta_trna_all,bakta_trna_uniq,bakta_16s,bakta_23s,bakta_5s" > {output.df}
        for file in {loc}/{sample}/tmp/annotation/bins/*; do
        name=$(basename $file )
        CDS_all=$(awk -F "\t" '{{ if ($2 == "cds") {{print $2}} }}' $file/${{name}}.tsv | grep -c "cds" -) || true
        CDS_hyp=$(awk -F "\t" '{{ if ($2 == "cds") {{print $8}} }}' $file/${{name}}.tsv | grep -c "hypothetical protein" -) || true
        CDS_dens=$(grep "coding density:" $file/${{name}}.txt | cut -f3 -d" ") || true
        tRNA_all=$(awk -F "\t" '{{ if ($2 == "tRNA") {{print $7}} }}' $file/${{name}}.tsv | grep -c "trn" -) || true
        tRNA_uniq=$(awk -F "\t" '{{ if ($2 == "tRNA") {{print $7}} }}' $file/${{name}}.tsv | sort -u - | grep -c "trn" -) || true
        rRNA_16S=$(awk -F "\t" '{{ if ($7 == "rrs") {{print $2}} }}' $file/${{name}}.tsv | grep -c "rRNA" -) || true
        rRNA_23S=$(awk -F "\t" '{{ if ($7 == "rrl") {{print $2}} }}' $file/${{name}}.tsv | grep -c "rRNA" -) || true
        rRNA_5S=$(awk -F "\t" '{{ if ($7 == "rrf") {{print $2}} }}' $file/${{name}}.tsv | grep -c "rRNA" -) || true
        echo "$name,$CDS_all,$CDS_hyp,$CDS_dens,$tRNA_all,$tRNA_uniq,$rRNA_16S,$rRNA_23S,$rRNA_5S" >> {output.df}; done

        if [ -d "{loc}/{sample}/results/bakta" ]; then rm -r "{loc}/{sample}/results/bakta"; fi
        mkdir {loc}/{sample}/results/bakta
        rsync -r {loc}/{sample}/tmp/annotation/bins/* {loc}/{sample}/results/bakta/.

        if grep -q "Bakta," {loc}/{sample}/tmp/db_mmlong2-proc.csv; then sed -i '/Bakta,/d' {loc}/{sample}/tmp/db_mmlong2-proc.csv; fi
        echo "Bakta,{params.db}" >> {loc}/{sample}/tmp/db_mmlong2-proc.csv
        if ! grep -q "bakta" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "bakta " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi  
        """

rule Taxonomy_aggregate:
    conda: config["env_10"]
    input:
        rrna=os.path.join(loc, sample, "tmp/taxa/rrna/rrna.tsv"),
        metabuli=os.path.join(loc, sample, "tmp/taxa/metabuli/metabuli_c.tsv"),
        gtdb=os.path.join(loc, sample, "tmp/taxa/gtdb/gtdb.tsv"),
    output:
        contigs=expand("{loc}/{sample}/tmp/taxa/contigs_taxonomy.tsv",sample=sample,loc=loc),
        bins=expand("{loc}/{sample}/tmp/taxa/bins_taxonomy.tsv",sample=sample,loc=loc)
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_taxonomy_aggregate.tsv")
    shell:
        """
        R --no-echo --silent --args << 'df' > /dev/null 2>&1
        # Load and wrangle
        metabat <- read.delim("{loc}/{sample}/tmp/binning/mapping/cov_all.tsv", sep="\t", header=T)
        metabat <- metabat[,c(1,2,3)]
        colnames(metabat) <- c("contig","len_bp","cov")

        assembly_info <- read.delim("{loc}/{sample}/tmp/filtering/assembly_info.tsv", sep="\t", header=F)
        assembly_info <- assembly_info[,c(1,4)]
        colnames(assembly_info) <- c("contig","status_circular")
        assembly_info <- merge(metabat,assembly_info,by="contig")
        
        metabuli <- read.delim("{input.metabuli}", sep="\t", header=F)
        colnames(metabuli) <- c("contig","metabuli_tax")
	
        rrna <- read.delim("{input.rrna}", sep="\t", header=F)
        rrna <- cbind(rrna, data.frame(do.call('rbind', strsplit(as.character(rrna$V1),':',fixed=TRUE))))
        rrna <- rrna[c("X3","V2","V3","V4")]
        colnames(rrna) <- c("contig","rrna_tax","rrna_tophit","rrna_alnlen")
        
        links <- read.delim("{loc}/{sample}/tmp/binning/contig_bin.tsv", sep="\t", header=F)
        colnames(links) <- c("contig","bin")
        
        gtdb <- read.delim("{input.gtdb}", sep="\t", header=T)
        gtdb <- gtdb[,c(1,2,8,11,12,17,19,20)]
        colnames(gtdb)[1] <- "bin"
        colnames(gtdb)[2] <- "gtdb_tax"
        colnames(gtdb)[3] <- "gtdb_ref"
        colnames(gtdb)[4] <- "gtdb_ani"
        colnames(gtdb)[5] <- "gtdb_af"
        colnames(gtdb)[6] <- "gtdb_msa"
        colnames(gtdb)[7] <- "gtdb_red"
        colnames(gtdb)[8] <- "gtdb_warning"
        
        # SSU rRNA taxonomy for contigs
        rrna_tophit <- rrna
        rrna_tophit$id <- paste(rrna_tophit$contig,rrna_tophit$rrna_tax,sep="_")
        rrna_tophit <- rrna_tophit[c("id","rrna_tophit","rrna_alnlen")]
        rrna_tophit <- rrna_tophit[order(rrna_tophit$rrna_alnlen,decreasing=TRUE),]
        rrna_tophit <- rrna_tophit[order(rrna_tophit$rrna_tophit,decreasing=TRUE),]
        rrna_tophit <- rrna_tophit[!duplicated(rrna_tophit$id),]
        
        rrna_contig <- aggregate(rrna$contig, by=list(rrna$contig, rrna$rrna_tax), FUN=length)
        rrna_contig$id <- paste(rrna_contig$Group.1,rrna_contig$Group.2,sep="_")
        rrna_contig <- merge(rrna_contig,rrna_tophit, by="id")
 
        rrna_contig <- rrna_contig[order(rrna_contig$rrna_alnlen,decreasing=TRUE),]
        rrna_contig <- rrna_contig[order(rrna_contig$rrna_tophit,decreasing=TRUE),]
        rrna_contig <- rrna_contig[!duplicated(rrna_contig$Group.1),]
        
        rrna_contig <- rrna_contig[,c(2,3,5,6)]
        colnames(rrna_contig) <- c("contig","rrna_tax","rrna_tophit","rrna_alnlen")
        
        # Make contigs dataframe
        contigs <- merge(assembly_info,metabuli, by="contig", all=TRUE)
        contigs <- merge(contigs,rrna_contig, by="contig", all=TRUE)
        contigs <- merge(contigs,links, by="contig", all=TRUE)
        write.table(contigs,"{output.contigs}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        
        # SSU rRNA taxonomy for bins
        rrna_bin <- merge(rrna,links, by="contig")
        rrna_bin$bin <- as.factor(rrna_bin$bin)
        
        rrna_bin_tophit <- rrna_bin
        rrna_bin_tophit$id <- paste(rrna_bin_tophit$bin,rrna_bin_tophit$rrna_tax,sep="_")
        rrna_bin_tophit <- rrna_bin_tophit[c("id","rrna_tophit","rrna_alnlen")]
        rrna_bin_tophit <- rrna_bin_tophit[order(rrna_bin_tophit$rrna_alnlen,decreasing=TRUE),]
        rrna_bin_tophit <- rrna_bin_tophit[order(rrna_bin_tophit$rrna_tophit,decreasing=TRUE),]
        rrna_bin_tophit <- rrna_bin_tophit[!duplicated(rrna_bin_tophit$id),]
        
        rrna_bin <- aggregate(rrna_bin$bin, by=list(rrna_bin$bin, rrna_bin$rrna_tax), FUN=length)
        rrna_bin$id <- paste(rrna_bin$Group.1,rrna_bin$Group.2,sep="_")
        rrna_bin <- merge(rrna_bin,rrna_bin_tophit, by="id")	
        
        rrna_bin <- rrna_bin[order(rrna_bin$rrna_alnlen,decreasing=TRUE),]
        rrna_bin <- rrna_bin[order(rrna_bin$rrna_tophit,decreasing=TRUE),]
        rrna_bin <- rrna_bin[!duplicated(rrna_bin$Group.1),]
        
        rrna_bin <- rrna_bin[,c(2,3,5,6)]	
        colnames(rrna_bin) <- c("bin","rrna_tax","rrna_tophit","rrna_alnlen")

        # Make bins dataframe
        bins <- merge(gtdb,rrna_bin, by="bin", all=TRUE)
        write.table(bins,"{output.bins}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Taxonomy_bins:
    conda: config["env_11"]
    input:
        expand("{loc}/{sample}/results/bins/{bin}.fa",sample=sample,loc=loc,bin=get_bins())
    output:
        os.path.join(loc, sample, "tmp/taxa/gtdb/gtdb.tsv")
    params:
        db=config["db_gtdb"]
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_taxonomy_gtdb.tsv")
    resources: usage=100
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/taxa" ]; then mkdir {loc}/{sample}/tmp/taxa; fi
        if [ -d "{loc}/{sample}/tmp/taxa/gtdb" ]; then rm -r {loc}/{sample}/tmp/taxa/gtdb; fi
        export GTDBTK_DATA_PATH="{params.db}"
		
        org_dir=$(pwd) && cd {loc}/{sample}/tmp/taxa
        gtdbtk classify_wf --cpus {threads} --genome_dir {loc}/{sample}/results/bins --extension .fa --out_dir gtdb --tmpdir $(pwd) --skip_ani_screen
        cat {loc}/{sample}/tmp/taxa/gtdb/classify/gtdbtk.bac120.summary.tsv > {output}
        if [ -f "{loc}/{sample}/tmp/taxa/gtdb/classify/gtdbtk.ar53.summary.tsv" ]; then tail -n+2 {loc}/{sample}/tmp/taxa/gtdb/classify/gtdbtk.ar53.summary.tsv >> {output}; fi
        cd $org_dir

        if grep -q "GTDB," {loc}/{sample}/tmp/db_mmlong2-proc.csv; then sed -i '/GTDB,/d' {loc}/{sample}/tmp/db_mmlong2-proc.csv; fi
        echo "GTDB,{params.db}" >> {loc}/{sample}/tmp/db_mmlong2-proc.csv
        if ! grep -q "gtdbtk" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "gtdbtk " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi  
        """

rule Taxonomy_contigs:
    conda: config["env_10"]
    input:
        os.path.join(loc, sample, "results", f"{sample}_assembly.fasta")
    output:
        res=os.path.join(loc, sample, "tmp/taxa/metabuli/contigs_classifications.tsv"),
        res_filt=os.path.join(loc, sample, "tmp/taxa/metabuli/metabuli_c.tsv")
    params:
        db=config["db_metabuli"],
        link=config["link_metabuli"],
        apptainer=config["apptainer_status"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_taxonomy_metabuli.tsv")
    resources: usage=100
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/taxa" ]; then mkdir {loc}/{sample}/tmp/taxa; fi
        if [ -d "{loc}/{sample}/tmp/taxa/metabuli" ]; then rm -r {loc}/{sample}/tmp/taxa/metabuli; fi
        if [ "{params.apptainer}" == "FALSE" ]; then if [ ! -d $CONDA_PREFIX/metabuli ]; then wget -qO- {params.link} | tar xvz && mv metabuli $CONDA_PREFIX/.; fi; fi

        $CONDA_PREFIX/metabuli/bin/metabuli classify {input} {params.db} {loc}/{sample}/tmp/taxa/metabuli contigs --threads {threads} --seq-mode 3 --lineage 1
        cut -f2,7 {output.res} > {output.res_filt}

        if grep -q "metabuli," {loc}/{sample}/tmp/db_mmlong2-proc.csv; then sed -i '/metabuli,/d' {loc}/{sample}/tmp/db_mmlong2-proc.csv; fi
        echo "metabuli,{params.db}" >> {loc}/{sample}/tmp/db_mmlong2-proc.csv
        version=$(grep "### Update in" $CONDA_PREFIX/metabuli/README.md | cut -f4 -d" " | head -1)
        if ! grep -q "metabuli" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then echo "metabuli,${{version}}" >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi  
        """

rule Taxonomy_rrna:
    conda: config["env_10"]
    input:
        contigs=os.path.join(loc, sample, "results", f"{sample}_assembly.fasta"),
        links=os.path.join(loc, sample, "tmp/binning/contig_bin.tsv"),
    output:
        rrna=os.path.join(loc, sample, "tmp/taxa/rrna/rRNA.fa"),
        ssu=os.path.join(loc, sample, "tmp/taxa/rrna/rRNA_ssu.fa"),
        id=os.path.join(loc, sample, "tmp/taxa/rrna/ssu_id.txt"),
        ssu2=os.path.join(loc, sample, "tmp/taxa/rrna/rRNA_ssu_renamed.fa"),
        tax=os.path.join(loc, sample, "tmp/taxa/rrna/rrna.tsv"),
    params:
        min_len_ssu=config["min_len_ssu"],
        min_id_ssu=config["min_id_ssu"],
        ssu=config["ssu"],
        db=config["db_rrna"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_taxonomy_usearch.tsv")
    resources: usage=100
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/taxa" ]; then mkdir {loc}/{sample}/tmp/taxa; fi
        if [ -d "{loc}/{sample}/tmp/taxa/rrna" ]; then rm -r {loc}/{sample}/tmp/taxa/rrna; fi
        mkdir {loc}/{sample}/tmp/taxa/rrna
        if [ -f "{loc}/{sample}/results/{sample}_{params.ssu}.fa" ]; then rm {loc}/{sample}/results/{sample}_{params.ssu}.fa; fi
        
        barrnap {input.contigs} --threads {threads} --outseq {output.rrna}
        if [ -f "{input.contigs}.fai" ]; then rm {input.contigs}.fai; fi

        grep -A1 ">{params.ssu}" {output.rrna} > {output.ssu}
        grep ">" {output.ssu} | cut -c2- > {output.id}
        cat {output.ssu} > {output.ssu2}
        
        cat {output.id} | while read id; do
        contig=$(echo "$id" | cut -f3 -d":" -)
        if grep -wq $contig {input.links}; then bin=$(grep -w $contig {input.links} | cut -f2); else bin="{sample}.unbinned"; fi
        sed -i "s/${{id}}/${{bin}}:${{id}}/g" {output.ssu2}; done
        
        sed -i "s/::/:/g" {output.ssu2}
        rsync {output.ssu2} {loc}/{sample}/results/{sample}_{params.ssu}.fa
        usearch --usearch_global {output.ssu} --db {params.db} --threads {threads} --strand both --id {params.min_id_ssu} -mincols {params.min_len_ssu} --top_hit_only -maxaccepts 0 -maxrejects 0 --blast6out {output.tax}

        if grep -q "rRNA," {loc}/{sample}/tmp/db_mmlong2-proc.csv; then sed -i '/rRNA,/d' {loc}/{sample}/tmp/db_mmlong2-proc.csv; fi
        echo "rRNA,{params.db}" >> {loc}/{sample}/tmp/db_mmlong2-proc.csv
        if ! grep -q "usearch" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "usearch " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule Summary_usage:
    conda: config["env_10"]
    input:
        gen_stats=os.path.join(loc, sample, "tmp/stats/gen_stats.tsv"),
        bins_extraqc=os.path.join(loc, sample, "tmp/extra_qc/bins_extraqc.tsv"),
        bins_taxonomy=os.path.join(loc, sample, "tmp/taxa/bins_taxonomy.tsv"),
        bins_annotation=os.path.join(loc, sample, "tmp/annotation/bins_annotation.tsv"),
    output:
        os.path.join(loc, sample, "tmp/logs/summary_mmlong2-proc.tsv")
    threads: proc
    shell:
        """
        head -n1 {loc}/{sample}/tmp/logs/summary_mmlong2-lite.tsv > {output}.tmp
        for log in {loc}/{sample}/tmp/logs/usage_*.tsv; do
            name=$(basename $log)
            IFS=_ read _ stage step <<< "${{name%.tsv}}"
            sed 1d $log | sed "s/^/${{stage}}\t${{step}}\t/" >> {output}.tmp
        done

        awk 'BEGIN {{
        order["stage"]=1; order["assembly"]=2; order["polishing"]=3; order["curation"]=4; order["filtering"]=5; order["singletons"]=6; order["coverage"]=7; order["binning"]=8; order["summary"]=9; order["stats"]=10; order["extraqc"]=11; order["taxonomy"]=12; order["annotation"]=13;
        }}{{
            key = ($1 in order) ? order[$1] : 9999
            print key "\t" $0
        }}' {output}.tmp | sort -k1,1n | cut -f2- > {output}

        rm {output}.tmp

        if [ -f "{loc}/{sample}/tmp/logs/usage_summary_stats.tsv" ]; then
        rsync {output} {loc}/{sample}/results/{sample}_usage.tsv
        else cat {loc}/{sample}/tmp/logs/summary_mmlong2-lite.tsv <(tail -n +2 {output}) > {loc}/{sample}/results/{sample}_usage.tsv; fi
        """

rule Cleanup:
    conda: config["env_10"]
    input:
        os.path.join(loc, sample, "tmp/logs/summary_mmlong2-proc.tsv")
    output:
        os.path.join(loc, sample, "tmp/logs/cleanup-proc.txt")
    params:
        cleanup=config["cleanup_status"],
    threads: proc
    shell:
        """
        if [ "{params.cleanup}" == "TRUE" ]; then

        size_pre=$(du -sb {loc}/{sample}/tmp | awk '{{printf "%.1f", $1/1024/1024/1024}}')
        files_pre=$(find {loc}/{sample}/tmp -type f -printf '.' | wc -c)
        echo "Status before cleanup (MAG analysis): $files_pre files and $size_pre GB of storage" > {output}

        if [ -d "{loc}/{sample}/tmp/extra_qc/tmp" ]; then rm -r {loc}/{sample}/tmp/extra_qc/tmp; fi
        if [ -d "{loc}/{sample}/tmp/extra_qc/gunc/gene_calls" ]; then rm -r {loc}/{sample}/tmp/extra_qc/gunc/gene_calls; fi
        if [ -f "{loc}/{sample}/tmp/extra_qc/variants/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/extra_qc/variants/contigs.fasta; fi
        if ls {loc}/{sample}/tmp/extra_qc/variants/*.bam >/dev/null 2>&1; then rm {loc}/{sample}/tmp/extra_qc/variants/*.bam*; fi

        if [ -d "{loc}/{sample}/tmp/annotation/tmp" ]; then rm -r {loc}/{sample}/tmp/annotation/tmp; fi
        if [ -d "{loc}/{sample}/tmp/annotation/bins" ]; then rm -r {loc}/{sample}/tmp/annotation/bins; fi
        if [ -f "{loc}/{sample}/tmp/taxa/rrna/rRNA.fa" ]; then pigz {loc}/{sample}/tmp/taxa/rrna/*.fa; fi

        if [ -f "{loc}/{sample}/results/{sample}_bins_annotation.tsv" ]; then rm {loc}/{sample}/results/{sample}_bins_annotation.tsv; fi
        if [ -f "{loc}/{sample}/results/{sample}_bins_taxonomy.tsv" ]; then rm {loc}/{sample}/results/{sample}_*_taxonomy.tsv; fi
        if [ -f "{loc}/{sample}/results/{sample}_bins_extraqc.tsv" ]; then rm {loc}/{sample}/results/{sample}_*_extraqc.tsv; fi
        if [ -f "{loc}/{sample}/results/{sample}_gen_stats.tsv" ]; then rm {loc}/{sample}/results/{sample}_*_stats.tsv; fi

        if ls {loc}/{sample}/tmp/logs/usage_*.tsv >/dev/null 2>&1; then
        rsync {loc}/{sample}/tmp/logs/summary_mmlong2-proc.tsv {loc}/{sample}/tmp/logs/summary_mmlong2-proc_$(date +"%Y%m%d-%Hh%Mm%Ss").tsv 
        rm {loc}/{sample}/tmp/logs/usage_*.tsv; fi

        size_post=$(du -sb {loc}/{sample}/tmp | awk '{{printf "%.1f", $1/1024/1024/1024}}')
        files_post=$(find {loc}/{sample}/tmp -type f -printf '.' | wc -c)
        echo "Status after cleanup (MAG analysis): $files_post files and $size_post GB of storage" >> {output}

        else echo "Cleanup skipped for MAG analysis" > {output}; fi
        """
