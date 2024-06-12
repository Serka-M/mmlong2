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
    min_version("8.0.0")
    if not os.path.exists(config["db_kaiju"]): sys.exit(print("Kaiju database input (((",config["db_kaiju"],"))) not found. Aborting..."))
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
    if not os.path.exists(os.path.join(loc, sample, "tmp", "dep_mmlong2-proc.csv")): 
	    with open(os.path.join(loc, sample, "tmp", "dep_mmlong2-proc.csv"), 'w') as f:
		    f.write("dependency,version\n")

onsuccess:
    from datetime import datetime
    now = datetime.now()
    print("MAG processing with version",wf_v,"completed at",now.strftime("%Y/%m/%d %H:%M:%S"))
    print("Thank you for using mmlong2")

onerror:
    print("An error has occurred. Inspect Snakemake log files for troubleshooting.")

rule Finalise:
    conda: config["env_9"]
    input:
        gen_stats=expand("{loc}/{sample}/tmp/stats/gen_stats.tsv",sample=sample,loc=loc),
        contigs_stats=expand("{loc}/{sample}/tmp/stats/contigs_stats.tsv",sample=sample,loc=loc),
        contigs_extraqc=expand("{loc}/{sample}/tmp/extra_qc/contigs_extraqc.tsv",sample=sample,loc=loc),
        contigs_taxonomy=expand("{loc}/{sample}/tmp/taxa/contigs_taxonomy.tsv",sample=sample,loc=loc),
        bins=expand("{loc}/{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample,loc=loc),
        bins_extraqc=expand("{loc}/{sample}/tmp/extra_qc/bins_extraqc.tsv",sample=sample,loc=loc),
        bins_taxonomy=expand("{loc}/{sample}/tmp/taxa/bins_taxonomy.tsv",sample=sample,loc=loc),
        bins_annotation=expand("{loc}/{sample}/tmp/annotation/bins_annotation.tsv",sample=sample,loc=loc),
    output:
        dep=expand("{loc}/{sample}/results/dependencies.csv",sample=sample,loc=loc),
        gen=expand("{loc}/{sample}/results/{sample}_general.tsv",sample=sample,loc=loc),
        contigs=expand("{loc}/{sample}/results/{sample}_contigs.tsv",sample=sample,loc=loc),
        bins=expand("{loc}/{sample}/results/{sample}_bins.tsv",sample=sample,loc=loc),
    shell:
        """
        cat {loc}/{sample}/tmp/dep_mmlong2-lite.csv > {output.dep}
        tail -n+2 {loc}/{sample}/tmp/dep_mmlong2-proc.csv >> {output.dep}
        
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
        contigs$wf_mode <- "{mode}"
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
        
        bins$wf_name <- NULL
        bins$wf_mode <- NULL
        bins$wf_v <- NULL
        bins$wf_date <- NULL
        
        bins$wf_name <- "{sample}"
        bins$wf_mode <- "{mode}"
        bins$wf_v <- "{wf_v}"
        bins$wf_date <- Sys.Date()
        
        gen$contigs_circ <- nrow(contigs[(contigs$status_circular == "Y"), ])
        gen$all_bins <- nrow(bins)
        gen$circ_bins <- nrow(bins[(bins$cbin_status == "Y"), ])
        gen$hq_bins <- nrow(bins[(bins$bin_status == "HQ"), ])
        gen$mq_bins <- nrow(bins[(bins$bin_status == "MQ"), ])
        gen$lq_bins <- nrow(bins[(bins$bin_status == "LQ"), ])
        gen$contaminated_bins <- nrow(bins[(bins$bin_status == "Contaminated"), ])
        
        gen$bin_cov_median <- median(bins$cov)
        gen$yield_assembled <- round(gen$mapped_yield_gp*10**9/gen$reads_yield_bp*100,3)
        gen$yield_binned <- round(gen$yield_assembled*sum(bins$r_abund)/100,3)
        gen$assembly_binned <- round(sum(bins$genome_size)/gen$contigs_yield_bp*100,3)
        gen[,"mapped_yield_gp"] <- NULL
        
        gen$wf_name <- "{sample}"
        gen$wf_mode <- "{mode}"
        gen$wf_v <- "{wf_v}"
        gen$wf_date <- Sys.Date()
        
        # Save data
        write.table(gen,"{output.gen}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        write.table(contigs,"{output.contigs}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        write.table(bins,"{output.bins}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Stats_aggregate:
    conda: config["env_9"]
    input:
        reads=expand("{loc}/{sample}/tmp/stats/reads.txt",sample=sample,loc=loc),
        contigs=expand("{loc}/{sample}/tmp/stats/contigs.txt",sample=sample,loc=loc),
        gc=expand("{loc}/{sample}/tmp/stats/gc.tsv",sample=sample,loc=loc),
        map=expand("{loc}/{sample}/tmp/stats/map_filt.tsv",sample=sample,loc=loc),
    output:
        contigs=expand("{loc}/{sample}/tmp/stats/contigs_stats.tsv",sample=sample,loc=loc),
        gen=expand("{loc}/{sample}/tmp/stats/gen_stats.tsv",sample=sample,loc=loc)
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
        colnames(map) <- c("mapped_yield_gp", "mapped_identity_median", "mapped_identity_mean")
        
        gen <- cbind(reads,contigs_gen,map)
        write.table(gen,"{output.gen}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        
        # Contig-level stats
        contigs=read.delim("{input.gc}", sep="\t", header=T)
        colnames(contigs) <- c("contig", "gc")
        write.table(contigs,"{output.contigs}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Stats_seq:
    conda: config["env_12"]
    input:
        expand("{loc}/{sample}/results/{sample}_assembly.fasta",sample=sample,loc=loc),
        {fastq}
    output:
        reads=expand("{loc}/{sample}/tmp/stats/reads.txt",sample=sample,loc=loc),
        contigs=expand("{loc}/{sample}/tmp/stats/contigs.txt",sample=sample,loc=loc),
        gc=expand("{loc}/{sample}/tmp/stats/gc.tsv",sample=sample,loc=loc)
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/stats" ]; then mkdir {loc}/{sample}/tmp/stats; fi
        seqkit fx2tab -n --gc {loc}/{sample}/results/{sample}_assembly.fasta > {output.gc}
        nanoq -i {loc}/{sample}/results/{sample}_assembly.fasta -s -t {threads} --report {output.contigs}
        nanoq -i {fastq} -s -t {threads} --report {output.reads}
        if ! grep -q "nanoq" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "nanoq " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule Stats_map:
    conda: config["env_12"]
    input:
        expand("{loc}/{sample}/tmp/binning/mapping/1_{map}.bam",sample=sample,loc=loc,map=get_mapping(mode))
    output:
        innit=expand("{loc}/{sample}/tmp/stats/map.tsv",sample=sample,loc=loc),
        filt=expand("{loc}/{sample}/tmp/stats/map_filt.tsv",sample=sample,loc=loc),
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/stats" ]; then mkdir {loc}/{sample}/tmp/stats; fi
        cramino --threads {threads} {input} > {output.innit}
        awk 'BEGIN {{FS = "\t"}} {{col1[NR] = $1; col2[NR] = $2}} END {{print col1[4] "\t" col1[11] "\t" col1[12]; print col2[4] "\t" col2[11] "\t" col2[12];}}' {output.innit} > {output.filt}
        if ! grep -q "cramino" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "cramino " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule ExtraQC_aggregate:
    conda: config["env_9"]
    input:
        variants=expand("{loc}/{sample}/tmp/extra_qc/variants/variants.tsv",sample=sample,loc=loc),
        gunc=expand("{loc}/{sample}/tmp/extra_qc/gunc.tsv",sample=sample,loc=loc),
        links=expand("{loc}/{sample}/tmp/binning/contig_bin.tsv",sample=sample,loc=loc),
    output:
        contigs=expand("{loc}/{sample}/tmp/extra_qc/contigs_extraqc.tsv",sample=sample,loc=loc),
        bins=expand("{loc}/{sample}/tmp/extra_qc/bins_extraqc.tsv",sample=sample,loc=loc)
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
        
        gunc=read.delim("{input.gunc}", sep="\t", header=T)
        colnames(gunc) <- c("bin", "gunc_css", "gunc_rrs", "gunc_pass")
        
        variants <- merge(variants,links,by="contig",all=TRUE)
        variants[is.na(variants$var_n),]$var_n <- 0
        
        bins <- aggregate(variants$var_n, by=list(variants$bin), FUN=sum)
        colnames(bins) <- c("bin", "var_n")
        
        bins <- merge(gunc,bins,by="bin")
        write.table(bins,"{output.bins}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule ExtraQC_variants_prep:
    conda: config["env_13"]
    input:
        expand("{loc}/{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample,loc=loc)
    output:
        list=expand("{loc}/{sample}/tmp/extra_qc/variants/contigs.txt",sample=sample,loc=loc),
        ids=expand("{loc}/{sample}/tmp/extra_qc/variants/contigs_{split}.txt",sample=sample,loc=loc,split=[f"{i:02}" for i in range(1, get_splits(proc,proc_sub,config["split_max"])+1)]),
        contigs=expand("{loc}/{sample}/tmp/extra_qc/variants/contigs.fasta",sample=sample,loc=loc),
    params:
        splits=get_splits(proc,proc_sub,config["split_max"]),
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
    conda: config["env_13"]
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
    resources: usage=config["longshot_usage"]
    threads: proc_sub
    shell:
        """
        samtools view --threads $(({threads} / 2)) -bh {loc}/{sample}/tmp/binning/mapping/1_{params.map}.bam $(cat {input}) | samtools sort --threads $(({threads} / 2)) --write-index -o {output.map}##idx##{output.map}.bai
        longshot --bam {output.map} --ref {loc}/{sample}/tmp/extra_qc/variants/contigs.fasta --out {output.var} -n --min_cov {params.min_cov} --min_mapq {params.min_mapq} --min_alt_frac {params.min_frac} --min_alt_count {params.min_count}
        """

rule ExtraQC_variants_sum:
    conda: config["env_13"]
    input:
        expand("{loc}/{sample}/tmp/extra_qc/variants/contigs_{split}.vcf",sample=sample,loc=loc,split=[f"{i:02}" for i in range(1, get_splits(proc,proc_sub,config["split_max"])+1)]),
    output:
        var_sort=expand("{loc}/{sample}/tmp/extra_qc/variants/variants_sorted.vcf",sample=sample,loc=loc),
        var_sum=expand("{loc}/{sample}/tmp/extra_qc/variants/variants.tsv",sample=sample,loc=loc),
    threads: proc
    shell:
        """
        grep -hv '^#' {input} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' > {output.var_sort}
        awk -F "\t" '{{ if ($7 == "PASS") {{print $1}} }}' {output.var_sort} | uniq -c - | awk '{{$1=$1;print}}' | awk -F" " '{{print $2,$1}}' OFS='\t' > {output.var_sum}
        if ! grep -q "longshot" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "longshot " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule ExtraQC_contamination:
    conda: config["env_12"]
    input:
        expand("{loc}/{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample,loc=loc)
    output:
        raw=expand("{loc}/{sample}/tmp/extra_qc/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv",sample=sample,loc=loc),
        trunc=expand("{loc}/{sample}/tmp/extra_qc/gunc.tsv",sample=sample,loc=loc)
    params:
        db=config["db_gunc"]
    resources: usage=100
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/extra_qc" ]; then mkdir {loc}/{sample}/tmp/extra_qc; fi
        if [ ! -d "{loc}/{sample}/tmp/extra_qc/tmp" ]; then mkdir {loc}/{sample}/tmp/extra_qc/tmp; fi
        if [ ! -d "{loc}/{sample}/tmp/extra_qc/gunc" ]; then rm -r {loc}/{sample}/tmp/extra_qc/gunc; fi

        gunc run --db_file {params.db} --input_dir {loc}/{sample}/results/bins --out_dir {loc}/{sample}/tmp/extra_qc/gunc --threads {threads} --temp_dir {loc}/{sample}/tmp/extra_qc/tmp
        cut -f1,8,12,13 {output.raw} > {output.trunc}
        if ! grep -q "gunc" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "gunc " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi
        """

rule Annotation_aggregate:
    conda: config["env_9"]
    input:
        bakta=expand("{loc}/{sample}/tmp/annotation/bakta_stats.csv",sample=sample,loc=loc),
        barrnap=expand("{loc}/{sample}/tmp/annotation/barrnap_stats.csv",sample=sample,loc=loc),
        trnascan=expand("{loc}/{sample}/tmp/annotation/trna_stats.csv",sample=sample,loc=loc)
    output:
        df=expand("{loc}/{sample}/tmp/annotation/bins_annotation.tsv",sample=sample,loc=loc)
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
    conda: config["env_9"]
    input:
        expand("{loc}/{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample,loc=loc)
    output:
        df=expand("{loc}/{sample}/tmp/annotation/barrnap_stats.csv",sample=sample,loc=loc)
    params:
        db=config["db_barrnap"]
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/annotation" ]; then mkdir {loc}/{sample}/tmp/annotation; fi
        if [ ! -d "{loc}/{sample}/tmp/annotation/tmp" ]; then mkdir {loc}/{sample}/tmp/annotation/tmp; fi
        if [ ! -d "{loc}/{sample}/tmp/annotation/trna" ]; then mkdir {loc}/{sample}/tmp/annotation/trna; fi
        if [ ! -d "{loc}/{sample}/tmp/annotation/bins" ]; then mkdir {loc}/{sample}/tmp/annotation/bins; fi
        
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
    conda: config["env_11"]
    input:
        "{loc}/{sample}/results/bins/{bin}.fa",
        "{loc}/{sample}/tmp/annotation/barrnap_stats.csv",
    output:
        "{loc}/{sample}/tmp/annotation/trna/trna_{bin}.txt"
    params:
        db=config["db_trnascan"],
        bin="{bin}"
    resources: usage=config["trnascan_usage"]
    threads: proc_sub
    shell:
        """
        if [ -f "{loc}/{sample}/tmp/annotation/trna/trna_{params.bin}.txt" ]; then rm {loc}/{sample}/tmp/annotation/trna/trna_{params.bin}.txt; fi
        if [ -f "{loc}/{sample}/tmp/annotation/trna/stats_{params.bin}.txt" ]; then rm {loc}/{sample}/tmp/annotation/trna/stats_{params.bin}.txt; fi
        tRNAscan-SE -{params.db} -o {loc}/{sample}/tmp/annotation/trna/trna_{params.bin}.txt -m {loc}/{sample}/tmp/annotation/trna/stats_{params.bin}.txt -d {loc}/{sample}/results/bins/{params.bin}.fa --thread {threads}
        sed -i -e '1,3'd -e "s/$/\t{params.bin}/g" {loc}/{sample}/tmp/annotation/trna/trna_{params.bin}.txt
        """

rule Annotation_trna_sum:
    conda: config["env_11"]
    input:
        expand("{loc}/{sample}/tmp/annotation/trna/trna_{bin}.txt",sample=sample,loc=loc,bin=get_bins())
    output:
        df=expand("{loc}/{sample}/tmp/annotation/trna_stats.csv",sample=sample,loc=loc)
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
    conda: config["env_11"]
    input:
        "{loc}/{sample}/results/bins/{bin}.fa",
        "{loc}/{sample}/tmp/annotation/barrnap_stats.csv",
    output:
        "{loc}/{sample}/tmp/annotation/bins/{bin}/{bin}.tsv"
    params:
        db=config["db_bakta"],
        extra=get_bakta_params(config["bakta_extra"]),
        bin="{bin}"	
    resources: usage=config["bakta_usage"]
    threads: proc_sub
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/annotation/bins/{params.bin}" ]; then rm -r {loc}/{sample}/tmp/annotation/bins/{params.bin}; fi
        export TMPDIR={loc}/{sample}/tmp/annotation/tmp
        bakta --compliant --db {params.db} --prefix {params.bin} --output {loc}/{sample}/tmp/annotation/bins/{params.bin} --keep-contig-headers --tmp-dir {loc}/{sample}/tmp/annotation/tmp --threads {threads} {loc}/{sample}/results/bins/{params.bin}.fa --meta --force {params.extra}
        bakta_plot --prefix {params.bin} --output {loc}/{sample}/tmp/annotation/bins/{params.bin} --force {loc}/{sample}/tmp/annotation/bins/{params.bin}/{params.bin}.json
        """

rule Annotation_bins_sum:
    conda: config["env_11"]
    input:
        expand("{loc}/{sample}/tmp/annotation/bins/{bin}/{bin}.tsv",sample=sample,loc=loc,bin=get_bins())
    output:
        df=expand("{loc}/{sample}/tmp/annotation/bakta_stats.csv",sample=sample,loc=loc)
    threads: proc
    shell:
        """
        echo "bin,bakta_cds_all,bakta_cds_hyp,bakta_trna_all,bakta_trna_uniq,bakta_16s,bakta_23s,bakta_5s" > {output.df}
        for file in {loc}/{sample}/tmp/annotation/bins/*; do
        name=$(basename $file )
        CDS_all=$(awk -F "\t" '{{ if ($2 == "cds") {{print $2}} }}' $file/${{name}}.tsv | grep -c "cds" -) || true
        CDS_hyp=$(awk -F "\t" '{{ if ($2 == "cds") {{print $8}} }}' $file/${{name}}.tsv | grep -c "hypothetical protein" -) || true
        tRNA_all=$(awk -F "\t" '{{ if ($2 == "tRNA") {{print $7}} }}' $file/${{name}}.tsv | grep -c "trn" -) || true
        tRNA_uniq=$(awk -F "\t" '{{ if ($2 == "tRNA") {{print $7}} }}' $file/${{name}}.tsv | sort -u - | grep -c "trn" -) || true
        rRNA_16S=$(awk -F "\t" '{{ if ($7 == "rrs") {{print $2}} }}' $file/${{name}}.tsv | grep -c "rRNA" -) || true
        rRNA_23S=$(awk -F "\t" '{{ if ($7 == "rrl") {{print $2}} }}' $file/${{name}}.tsv | grep -c "rRNA" -) || true
        rRNA_5S=$(awk -F "\t" '{{ if ($7 == "rrf") {{print $2}} }}' $file/${{name}}.tsv | grep -c "rRNA" -) || true
        echo "$name,$CDS_all,$CDS_hyp,$tRNA_all,$tRNA_uniq,$rRNA_16S,$rRNA_23S,$rRNA_5S" >> {output.df}; done

        if [ -d "{loc}/{sample}/results/bakta" ]; then rm -r "{loc}/{sample}/results/bakta"; fi
        mkdir {loc}/{sample}/results/bakta
        rsync -r {loc}/{sample}/tmp/annotation/bins/* {loc}/{sample}/results/bakta/.
        if ! grep -q "bakta" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "bakta " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi  
        """

rule Taxonomy_aggregate:
    conda: config["env_9"]
    input:
        rrna=expand("{loc}/{sample}/tmp/taxa/rrna/rrna.tsv",sample=sample,loc=loc),
        kaiju=expand("{loc}/{sample}/tmp/taxa/contigs/kaiju_tax_filt.tsv",sample=sample,loc=loc),
        gtdb=expand("{loc}/{sample}/tmp/taxa/gtdb/gtdb.tsv",sample=sample,loc=loc),
    output:
        contigs=expand("{loc}/{sample}/tmp/taxa/contigs_taxonomy.tsv",sample=sample,loc=loc),
        bins=expand("{loc}/{sample}/tmp/taxa/bins_taxonomy.tsv",sample=sample,loc=loc)
    shell:
        """
        grep -v '^#' {loc}/{sample}/tmp/assembly/assembly_info.txt > {loc}/{sample}/tmp/taxa/assembly_info.txt
        R --no-echo --silent --args << 'df' > /dev/null 2>&1
        # Load and wrangle
        assembly_info=read.delim("{loc}/{sample}/tmp/taxa/assembly_info.txt", sep="\t", header=F)
        assembly_info=assembly_info[,c(1,2,3,4)]
        colnames(assembly_info) <- c("contig","len_bp","cov","status_circular")
        
        kaiju=read.delim("{input.kaiju}", sep="\t", header=F)
        colnames(kaiju) <- c("contig","tax_kaiju")
        kaiju=as.data.frame(kaiju, stringsAsFactors = FALSE)
	
        rrna <- read.delim("{input.rrna}", sep="\t", header=F)
        rrna=cbind(rrna, data.frame(do.call('rbind', strsplit(as.character(rrna$V1),':',fixed=TRUE))))
        rrna=rrna[c("X3","V2","V3")]
        colnames(rrna) <- c("contig","tax_rrna","tophit_rrna")
        
        links <- read.delim("{loc}/{sample}/tmp/binning/contig_bin.tsv", sep="\t", header=F)
        colnames(links) <- c("contig","bin")
        
        gtdb <- read.delim("{input.gtdb}", sep="\t", header=T)
        gtdb <- gtdb[,c(1,2,8,11,12,17,19,20)]
        colnames(gtdb)[1] <- "bin"
        colnames(gtdb)[2] <- "tax_gtdb"
        colnames(gtdb)[3] <- "ref_gtdb"
        colnames(gtdb)[4] <- "ani_gtdb"
        colnames(gtdb)[5] <- "af_gtdb"
        colnames(gtdb)[6] <- "msa_gtdb"
        colnames(gtdb)[7] <- "red_gtdb"
        colnames(gtdb)[8] <- "warning_gtdb"
        
        # SSU rRNA taxonomy for contigs
        rrna_tophit <- rrna
        rrna_tophit$id <- paste(rrna_tophit$contig,rrna_tophit$tax_rrna,sep="_")
        rrna_tophit <- rrna_tophit[c("id","tophit_rrna")]
        rrna_tophit <- do.call(rbind, lapply(split(rrna_tophit,rrna_tophit$id), function(x) {{return(x[which.max(x$tophit_rrna),])}}))
        
        rrna_contig <- aggregate(rrna$contig, by=list(rrna$contig, rrna$tax_rrna), FUN=length)
        rrna_contig$id <- paste(rrna_contig$Group.1,rrna_contig$Group.2,sep="_")
        rrna_contig <- merge(rrna_contig,rrna_tophit, by="id")
        
        rrna_contig <- rrna_contig[order(rrna_contig$tophit_rrna,decreasing=TRUE),]
        rrna_contig <- rrna_contig[order(rrna_contig$Group.1,decreasing=TRUE),]
        rrna_contig <- do.call(rbind, lapply(split(rrna_contig,rrna_contig$Group.1), function(x) {{return(x[which.max(x$x),])}}))
        
        rrna_contig <- rrna_contig[,c(2,3,5)]
        colnames(rrna_contig) <- c("contig","tax_rrna","tophit_rrna")
        
        # Make contigs dataframe
        contigs <- merge(assembly_info,kaiju, by="contig", all=TRUE)
        contigs <- merge(contigs,rrna_contig, by="contig", all=TRUE)
        contigs <- merge(contigs,links, by="contig", all=TRUE)
        write.table(contigs,"{output.contigs}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        
        # SSU rRNA taxonomy for bins
        rrna_bin <- merge(rrna,links, by="contig")
        rrna_bin$bin <- as.factor(rrna_bin$bin)
        
        rrna_bin_tophit <- rrna_bin
        rrna_bin_tophit$id <- paste(rrna_bin_tophit$bin,rrna_bin_tophit$tax_rrna,sep="_")
        rrna_bin_tophit <- rrna_bin_tophit[c("id","tophit_rrna")]
        rrna_bin_tophit <- do.call(rbind, lapply(split(rrna_bin_tophit,rrna_bin_tophit$id), function(x) {{return(x[which.max(x$tophit_rrna),])}}))
        
        rrna_bin <- aggregate(rrna_bin$bin, by=list(rrna_bin$bin,	rrna_bin$tax_rrna), FUN=length)
        rrna_bin$id <- paste(rrna_bin$Group.1,rrna_bin$Group.2,sep="_")
        rrna_bin <- merge(rrna_bin,rrna_bin_tophit, by="id")	
        
        rrna_bin <- rrna_bin[order(rrna_bin$tophit_rrna,decreasing=TRUE),]
        rrna_bin <- rrna_bin[order(rrna_bin$Group.1,decreasing=TRUE),]
        rrna_bin <- do.call(rbind, lapply(split(rrna_bin,rrna_bin$Group.1), function(x) {{return(x[which.max(x$x),])}}))
        
        rrna_bin <- rrna_bin[,c(2,3,5)]	
        colnames(rrna_bin) <- c("bin","tax_rrna","tophit_rrna")

        # Make bins dataframe
        bins <- merge(gtdb,rrna_bin, by="bin", all=TRUE)
        write.table(bins,"{output.bins}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Taxonomy_bins:
    conda: config["env_10"]
    input:
        expand("{loc}/{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample,loc=loc)
    output:
        expand("{loc}/{sample}/tmp/taxa/gtdb/gtdb.tsv",sample=sample,loc=loc)
    params:
        db=config["db_gtdb"]
    resources: usage=100
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/taxa" ]; then mkdir {loc}/{sample}/tmp/taxa; fi
        if [ -d "{loc}/{sample}/tmp/taxa/gtdb" ]; then rm -r {loc}/{sample}/tmp/taxa/gtdb; fi
        export GTDBTK_DATA_PATH="{params.db}"
		
        gtdbtk classify_wf --cpus {threads} --genome_dir {loc}/{sample}/results/bins --extension .fa --out_dir {loc}/{sample}/tmp/taxa/gtdb --tmpdir {loc}/{sample}/tmp/taxa --skip_ani_screen
        cat {loc}/{sample}/tmp/taxa/gtdb/classify/gtdbtk.bac120.summary.tsv > {output}
        if [ -f "{loc}/{sample}/tmp/taxa/gtdb/classify/gtdbtk.ar53.summary.tsv" ]; then tail -n+2 {loc}/{sample}/tmp/taxa/gtdb/classify/gtdbtk.ar53.summary.tsv >> {output}; fi
        if ! grep -q "gtdbtk" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "gtdbtk " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi  
        """

rule Taxonomy_contigs:
    conda: config["env_9"]
    input:
        expand("{loc}/{sample}/results/{sample}_assembly.fasta",sample=sample,loc=loc)
    output:
        res=expand("{loc}/{sample}/tmp/taxa/contigs/kaiju.tsv",sample=sample,loc=loc),
        sum=expand("{loc}/{sample}/tmp/taxa/contigs/kaiju_summary_phylum.tsv",sample=sample,loc=loc),
        names=expand("{loc}/{sample}/tmp/taxa/contigs/kaiju_tax.tsv",sample=sample,loc=loc),
        names_filt=expand("{loc}/{sample}/tmp/taxa/contigs/kaiju_tax_filt.tsv",sample=sample,loc=loc),
    params:
        db=config["db_kaiju"]
    resources: usage=100
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/taxa" ]; then mkdir {loc}/{sample}/tmp/taxa; fi
        if [ -d "{loc}/{sample}/tmp/taxa/contigs" ]; then rm -r {loc}/{sample}/tmp/taxa/contigs; fi
        mkdir {loc}/{sample}/tmp/taxa/contigs
        
        kaiju -z {threads} -t {params.db}/nodes.dmp -f {params.db}/*.fmi -i {input} -o {output.res}
        kaiju2table -t {params.db}/nodes.dmp -n {params.db}/names.dmp -r phylum -o {output.sum} {output.res}
        kaiju-addTaxonNames -r superkingdom,phylum,class,order,family,genus,species -t {params.db}/nodes.dmp -n {params.db}/names.dmp -i {output.res} -o {output.names}
        awk -F "\t" -v OFS="\t" '{{ if ($1 == "C") {{print $2, $4}} }}' {output.names} > {output.names_filt}
        if ! grep -q "kaiju" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "kaiju " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi  
        """

rule Taxonomy_rrna:
    conda: config["env_9"]
    input:
        contigs=expand("{loc}/{sample}/results/{sample}_assembly.fasta",sample=sample,loc=loc),
        links=expand("{loc}/{sample}/tmp/binning/contig_bin.tsv",sample=sample,loc=loc),
    output:
        rrna=expand("{loc}/{sample}/tmp/taxa/rrna/rRNA.fa",sample=sample,loc=loc),
        ssu=expand("{loc}/{sample}/tmp/taxa/rrna/rRNA_ssu.fa",sample=sample,loc=loc),
        id=expand("{loc}/{sample}/tmp/taxa/rrna/ssu_id.txt",sample=sample,loc=loc),
        ssu2=expand("{loc}/{sample}/tmp/taxa/rrna/rRNA_ssu_renamed.fa",sample=sample,loc=loc),
        tax=expand("{loc}/{sample}/tmp/taxa/rrna/rrna.tsv",sample=sample,loc=loc),
    params:
        min_len_ssu=config["min_len_ssu"],
        min_id_ssu=config["min_id_ssu"],
        ssu=config["ssu"],
        db=config["db_rrna"],
    resources: usage=100
    threads: proc
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/taxa" ]; then mkdir {loc}/{sample}/tmp/taxa; fi
        if [ -d "{loc}/{sample}/tmp/taxa/rrna" ]; then rm -r {loc}/{sample}/tmp/taxa/rrna; fi
        mkdir {loc}/{sample}/tmp/taxa/rrna
        if [ -f "{loc}/{sample}/results/{sample}_{params.ssu}.fa" ]; then rm {loc}/{sample}/results/{sample}_{params.ssu}.fa; fi
        
        barrnap {input.contigs} --threads {threads} --outseq {output.rrna}
        grep -A1 ">{params.ssu}" {output.rrna} > {output.ssu}
        grep ">" {output.ssu} | cut -c2- > {output.id}
        cat {output.ssu} > {output.ssu2}
        
        cat {output.id} | while read id; do
        contig=$(echo "$id" | cut -f3 -d":" -)
        if grep -wq $contig {input.links}; then bin=$(grep -w $contig {input.links} | cut -f2); else bin="{sample}.unbinned"; fi
        sed -i "s/${{id}}/${{bin}}:${{id}}/g" {output.ssu2}; done
        
        sed -i "s/::/:/g" {output.ssu2}
        rsync {output.ssu2} {loc}/{sample}/results/{sample}_{params.ssu}.fa
        
        vsearch --usearch_global {output.ssu} --db {params.db} --threads {threads} --minseqlength {params.min_len_ssu} --strand both --id {params.min_id_ssu} --top_hits_only --maxhits 1 --blast6out {output.tax}
        if ! grep -q "vsearch" {loc}/{sample}/tmp/dep_mmlong2-proc.csv; then conda list | grep -w "vsearch " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-proc.csv; fi  
        """
