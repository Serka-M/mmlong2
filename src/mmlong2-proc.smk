# DESCRIPTION: Snakemake workflow for analysing metagenome assembled genomes
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

sing=config["sing"]
wf_v=config["version"]
mode=config["mode"]
modes=['Nanopore-simplex', 'PacBio-HiFi']
proc=config["proc"]
fastq=config["fastq"]
sample=config["sample"]
minimap_ram=config["minimap_ram"]
medaka_split=config["medaka_split"]
medak_mod_pol=config["medak_mod_pol"]
bakta_split=config["bakta_split"]
db_bakta=config["db_bakta"]
db_gtdb=config["db_gtdb"]
db_kaiju=config["db_kaiju"]
db_midas=config["db_midas"]
db_silva=config["db_silva"]
db_gunc=config["db_gunc"]
db_barrnap=config["db_barrnap"]
dbs_barrnap=['bac', 'arc', 'euk', 'mito']

singularity: sing
shell.executable("/bin/bash")

onstart:
    from snakemake.utils import min_version
    import os
    import sys
    min_version("7.0.0")
    if not os.path.exists(db_kaiju): sys.exit(print("Bakta database input (((",db_kaiju,"))) not found. Aborting..."))
    if not os.path.exists(db_bakta): sys.exit(print("Bakta database input (((",db_bakta,"))) not found. Aborting..."))
    if not os.path.exists(db_gtdb): sys.exit(print("GTDB-tk database input (((",db_gtdb,"))) not found. Aborting..."))
    if not os.path.exists(db_silva): sys.exit(print("SILVA database input (((",db_silva,"))) not found. Aborting..."))
    if not os.path.exists(db_midas): sys.exit(print("MiDAS database input (((",db_midas,"))) not found. Aborting..."))
    if not os.path.exists(db_gunc): sys.exit(print("GUNC database input (((",db_gunc,"))) not found. Aborting..."))
    if not os.path.exists(fastq): sys.exit(print("Read input (((",fastq,"))) not found. Aborting..."))
    if not db_barrnap in dbs_barrnap: sys.exit(print("Provided Barrnap database (((",db_barrnap,"))) not recognised. Aborting..."))
    if not mode in modes: sys.exit(print("Provided workflow mode (((",mode,"))) not recognised. Aborting..."))

onsuccess:
    from datetime import datetime
    now = datetime.now()
    print("MAG processing with version",wf_v,"completed at",now.strftime("%Y/%m/%d %H:%M:%S"))
    print("Thank you for using mmlong2")

onerror:
    print("An error has occurred. Inspect Snakemake log files for troubleshooting.")

rule Finalise:
    conda: "env_0"
    input:
        expand("{sample}/tmp/taxonomy/contigs_taxonomy.tsv",sample=sample),
        expand("{sample}/tmp/taxonomy/bins_taxonomy.tsv",sample=sample),
        expand("{sample}/tmp/annotation/bins_annotation.tsv",sample=sample),
        expand("{sample}/tmp/stats/nanoq.ssf",sample=sample),
	expand("{sample}/tmp/stats/assembly.csv",sample=sample),
	expand("{sample}/tmp/stats/gc.tsv",sample=sample),
        expand("{sample}/tmp/extra_qc/var.tsv",sample=sample),
        expand("{sample}/tmp/extra_qc/gunc.tsv",sample=sample)
    output:
        dep1=expand("{sample}/tmp/dep_mmlong2.csv",sample=sample),
        dep2=expand("{sample}/results/dependencies.csv",sample=sample),
        gen=expand("{sample}/results/{sample}_general.tsv",sample=sample),
        con=expand("{sample}/results/{sample}_contigs.tsv",sample=sample),
        bin=expand("{sample}/results/{sample}_bins.tsv",sample=sample)
    shell:
        """
	if [ -f "$(pwd)/{sample}/results/{sample}_bins.tsv" ]; then rm {sample}/results/{sample}_bins.tsv; fi
	if [ -f "$(pwd)/{sample}/results/dependencies.csv" ]; then rm {sample}/results/dependencies.csv; fi

        printf "software,version\n" > {output.dep1}
	printf "mmlong2,{wf_v}\n" >> {output.dep1}
	printf "Kaiju,1.9.2\n" >> {output.dep1}
	printf "Bakta,1.6.1\n" >> {output.dep1}
	printf "Gunc,1.0.5\n" >> {output.dep1}
	printf "NanoQ,0.9.0\n" >> {output.dep1}
	printf "Vsearch,2.22.1\n" >> {output.dep1}
	printf "Barrnap,0.9\n" >> {output.dep1}
	printf "GTDBTk,2.1.1\n" >> {output.dep1}
	cat {output.dep1} > {output.dep2}
	awk '{{if(NR>2)print}}' {sample}/tmp/binning/dep_mmlong2-lite.csv >> {output.dep2}

	R --slave --silent --args << 'df'
	# Load data
	bins <- read.delim("{sample}/tmp/binning/bins_mmlong2-lite.tsv", sep="\t", header=T)
	annot <- read.delim("{sample}/tmp/annotation/bins_annotation.tsv", sep="\t", header=T)
	gunc <- read.delim("{sample}/tmp/extra_qc/gunc.tsv", sep="\t", header=T)
	taxa <- read.delim("{sample}/tmp/taxonomy/bins_taxonomy.tsv", sep="\t", header=T)

	contigs <- read.delim("{sample}/tmp/taxonomy/contigs_taxonomy.tsv", sep="\t", header=T)
	var <- read.delim("{sample}/tmp/extra_qc/var.tsv", sep="\t", header=F)
	gc <- read.delim("{sample}/tmp/stats/gc.tsv", sep="\t", header=F)

	nanoq <- read.delim("{sample}/tmp/stats/nanoq.ssf", sep=" ", header=T)
	flye <- read.delim("{sample}/tmp/stats/assembly.csv", sep=",", header=T)

	# Combine data
	colnames(gc) <- c("contig","GC")
	colnames(var) <- c("contig","var_n")
	contigs <- merge(contigs,gc,by="contig", all=TRUE)
	contigs <- merge(contigs,var,by="contig", all=TRUE)
	contigs[is.na(contigs$var_n),]$var_n <- 0
	contigs$var_perc <- round(contigs$var_n/contigs$len_bp*100,3)
	contigs$wf_name <- "{sample}"
	contigs$wf_mode <- "{mode}"
	contigs$wf_v <- "{wf_v}"
	contigs$wf_date <- Sys.Date()

	bins <- merge(bins,annot, by="bin")
	bins$MAG_status <- ifelse((bins$Completeness >= 90 & bins$Contamination <= 5 & bins$bakta_tRNA_uniq >= 18 &
              (bins$bakta_5S >= 1 | bins$barrnap_5S >= 1) & (bins$bakta_16S >= 1 | bins$barrnap_16S >= 1) &
              (bins$bakta_23S >= 1 | bins$barrnap_23S >= 1)),"HQ",
              ifelse(bins$Completeness >= 50 & bins$Contamination <= 10, "MQ",
              ifelse(bins$Completeness <= 50 & bins$Contamination <= 10, "LQ", "Contaminated")))
	bins$cMAG_status <- ifelse(grepl("bin.c",bins$bin),"Y","N")
	colnames(gunc) <- c("bin","gunc_contamination","gunc_status")
	bins <- merge(bins,gunc, by="bin")
	bins <- merge(bins,taxa, by="bin")
	var_bins <- contigs[!is.na(contigs$bin),]
	var_bins <- aggregate(var_bins$var_n, by=list(var_bins$bin), FUN=sum)
	colnames(var_bins) <- c("bin","var_n")
	bins <- merge(bins,var_bins, by="bin")
	bins$var_perc <- round(bins$var_n/bins$Genome_Size*100,3)
	bins$wf_name <- NULL
	bins$wf_mode <- NULL
	bins$wf_v <- NULL
	bins$wf_date <- NULL
	bins$wf_name <- "{sample}"
	bins$wf_mode <- "{mode}"
	bins$wf_v <- "{wf_v}"
	bins$wf_date <- Sys.Date()

	gen <- cbind(nanoq,flye)
	gen$contigs_circ <- nrow(contigs[(contigs$status_circular == "Y"), ])
	gen$all_mags <- nrow(bins)
	gen$circ_mags <- nrow(bins[(bins$cMAG_status == "Y"), ])
	gen$hq_mags <- nrow(bins[(bins$MAG_status == "HQ"), ])
	gen$mq_mags <- nrow(bins[(bins$MAG_status == "MQ"), ])
	gen$lq_mags <- nrow(bins[(bins$MAG_status == "LQ"), ])
	gen$contaminated_mags <- nrow(bins[(bins$MAG_status == "Contaminated"), ])
	gen$asm_binned <- sum(bins$Genome_Size)/gen$assembly_size_bp*100
	gen$r_abund_bins <- sum(bins$r_abund)
	gen$mag_cov_med_hq_mq <- median(bins[(bins$MAG_status == "MQ" | bins$MAG_status == "HQ"), ]$cov)
	gen$mag_cov_med_hq <- median(bins[(bins$MAG_status == "HQ"), ]$cov)
	gen$mag_cov_med_mq <- median(bins[(bins$MAG_status == "MQ"), ]$cov)
	gen$mag_cov_mad_hq_mq <- mad(bins[(bins$MAG_status == "MQ" | bins$MAG_status == "HQ"), ]$cov)
	gen$mag_cov_mad_hq <- mad(bins[(bins$MAG_status == "HQ"), ]$cov)
	gen$mag_cov_mad_mq <- mad(bins[(bins$MAG_status == "MQ"), ]$cov)
	gen$contigs_med_hq_mq <- median(bins[(bins$MAG_status == "MQ" | bins$MAG_status == "HQ"), ]$contigs)
	gen$contigs_med_hq <- median(bins[(bins$MAG_status == "HQ"), ]$contigs)
	gen$contigs_med_mq <- median(bins[(bins$MAG_status == "MQ"), ]$contigs)
	gen$wf_name <- "{sample}"
	gen$wf_mode <- "{mode}"
	gen$wf_v <- "{wf_v}"
	gen$wf_date <- Sys.Date()

	# Save data
	write.table(gen,"{output.gen}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
	write.table(contigs,"{output.con}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
	write.table(bins,"{output.bin}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule ExtraQC_microdiversity:
    conda: "env_5"
    input:
        expand("{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample)
    output:
        contigs=expand("{sample}/tmp/binning/contigs_binned.fasta",sample=sample),
        ids=expand("{sample}/tmp/binning/contigs_binned.txt",sample=sample),
        var=expand("{sample}/tmp/extra_qc/var.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/extra_qc" ]; then mkdir {sample}/tmp/extra_qc; fi
	cat {sample}/tmp/binning/bins/*.fa > {output.contigs}
	grep ">" {output.contigs} | cut -c 2- | awk '{{print $1}}' > {output.ids}

	if [ {mode} == "Nanopore-simplex" ]; then

	if [ ! -d "$(pwd)/{sample}/tmp/extra_qc/medaka" ]; then mkdir {sample}/tmp/extra_qc/medaka; fi
	n_contigs=$(grep ">" {output.contigs} | wc -l)	
	splits=$(({proc} / 2))
	if [ $splits -gt {medaka_split} ]; then splits={medaka_split}; fi;
	if [ $splits -gt $n_contigs ]; then splits=$n_contigs; fi;
	split -n l/$splits --numeric-suffixes=1 --additional-suffix=.txt -d {output.ids} {sample}/tmp/extra_qc/medaka/contigs_id_
	for file in {sample}/tmp/extra_qc/medaka/contigs_id_*.txt; do filename=lin_$(basename $file)
	awk 'BEGIN {{ ORS = " " }} {{ print }}' $file > {sample}/tmp/extra_qc/medaka/$filename; done
	mini_align -I {minimap_ram}G -i {fastq} -r {sample}/tmp/polishing/asm_pol_lenfilt.fasta -m -p {sample}/tmp/extra_qc/medaka/calls_to_draft -t {proc}
	find {sample}/tmp/extra_qc/medaka -type f -name "lin_*" | xargs -i --max-procs=$splits bash -c 'list=$(head -1 {{}}) && file=$(basename {{}} | sed 's/\.txt//') && medaka consensus {sample}/tmp/extra_qc/medaka/calls_to_draft.bam {sample}/tmp/extra_qc/medaka/$file.hdf --model {medak_mod_pol} --batch 200 --threads 2 --region $list'
	medaka variant {sample}/tmp/polishing/asm_pol_lenfilt.fasta {sample}/tmp/extra_qc/medaka/*.hdf {sample}/tmp/extra_qc/medaka/medaka_var.vcf
	cat {sample}/tmp/extra_qc/medaka/medaka_var.vcf | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' > {sample}/tmp/extra_qc/medaka/medaka_var_sort.vcf
	awk -F "\t" '{{ if ($7 == "PASS") {{print $1}} }}' {sample}/tmp/extra_qc/medaka/medaka_var_sort.vcf | uniq -c - | tr -d " \t" | sed 's/contig/,contig/g' | awk -F, '{{print $2,$1}}' OFS='\t' > {output.var}

	else

	awk -F, '{{print $1,0}}' OFS='\t' {output.ids} > {output.var}
	
	fi	
        """

rule ExtraQC_contamination:
    conda: "env_4"
    input:
        expand("{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample)
    output:
        raw=expand("{sample}/tmp/extra_qc/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv",sample=sample),
        trunc=expand("{sample}/tmp/extra_qc/gunc.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/extra_qc" ]; then mkdir {sample}/tmp/extra_qc; fi
	if [ ! -d "$(pwd)/{sample}/tmp/extra_qc/tmp" ]; then mkdir {sample}/tmp/extra_qc/tmp; fi

	gunc run --db_file {db_gunc} --input_dir {sample}/tmp/binning/bins --out_dir {sample}/tmp/extra_qc/gunc --threads {proc} --temp_dir {sample}/tmp/extra_qc/tmp
	cut -f1,9,13 {output.raw} > {output.trunc}
        """

rule Stats_summary:
    conda: "env_4"
    input:
        expand("{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample),
        expand("{sample}/tmp/flye/flye.log",sample=sample),
        {fastq}
    output:
        reads=expand("{sample}/tmp/stats/nanoq.ssf",sample=sample),
	asm=expand("{sample}/tmp/stats/assembly.csv",sample=sample),
	gc=expand("{sample}/tmp/stats/gc.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/stats" ]; then mkdir {sample}/tmp/stats; fi

	echo "reads_n reads_size_bp reads_N50_bp reads_len_max_bp reads_len_min_bp reads_mean_len_bp reads_median_len_bp reads_mean_q reads_median_q" > {output.reads}
	nanoq -s -i {fastq} >> {output.reads}

	contig_sum=$(grep '	Total length:	' {sample}/tmp/flye/flye.log | sed 's/Total length://g') 
	contig_count=$(grep '	Fragments:	' {sample}/tmp/flye/flye.log | sed 's/Fragments://g') 
	contig_n50=$(grep '	Fragments N50:	' {sample}/tmp/flye/flye.log | sed 's/Fragments N50://g')
	contig_cov=$(grep '	Mean coverage:	' {sample}/tmp/flye/flye.log | sed 's/Mean coverage://g') 
	contig_aln=$(grep 'Aligned read sequence' {sample}/tmp/flye/flye.log | sed 's/.*(//g' | sed 's/)//g') 

	echo "assembly_size_bp,assembly_contigs_n,assembly_contigs_N50_bp,assembly_contigs_meanCOV,assembly_reads_mapped" > {output.asm}
	echo "$contig_sum,$contig_count,$contig_n50,$contig_cov,$contig_aln" >> {output.asm}

	seqkit fx2tab -n --gc {sample}/tmp/polishing/asm_pol_lenfilt.fasta > {output.gc}
        """

rule Annotation_aggregate:
    conda: "env_0"
    input:
        expand("{sample}/tmp/annotation/bakta_stats.csv",sample=sample),
        expand("{sample}/tmp/annotation/barrnap_stats.csv",sample=sample)
    output:
        df=expand("{sample}/tmp/annotation/bins_annotation.tsv",sample=sample)
    shell:
        """
	R --slave --silent --args << 'df'
	bakta=read.delim("{sample}/tmp/annotation/bakta_stats.csv", sep=",", header=T)
	barrnap=read.delim("{sample}/tmp/annotation/barrnap_stats.csv", sep=",", header=T)
	bakta=merge(bakta,barrnap,by="bin")
	write.table(bakta,"{output.df}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Annotation_16S:
    conda: "env_1"
    input:
        expand("{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample)
    output:
        df=expand("{sample}/tmp/annotation/barrnap_stats.csv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/annotation" ]; then mkdir {sample}/tmp/annotation; fi

	function rrna_stats {{
	FILE=$1
	FILE_NAME=${{FILE##*/}}
	FILE_NAME=${{FILE_NAME%%.*}}
	BIN=$(basename $FILE | sed 's/\.fa//')
	S=`barrnap --threads {proc} --kingdom {db_barrnap} --quiet $FILE |\
	awk -F "=" -v bin=$BIN '
	NR == FNR {{a[$1]=0; next}}
	/^[^#]/{{gsub(/ .*/, "", $3); a[$3]++}}
	END{{OFS = ","; print bin, a["16S"], a["23S"], a["5S"]}}
	' <(printf "%s\n" 16S 23S 5S) -` 
	echo "$S"; }}

    	echo "bin, barrnap_16S, barrnap_23S, barrnap_5S" > {output.df}
	find {sample}/tmp/binning/bins -name "*.fa" | while read file; do rrna_stats "$file" >> {output.df}; done
        """

rule Annotation_MAGs:
    conda: "env_3"
    input:
        expand("{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample)
    output:
        df=expand("{sample}/tmp/annotation/bakta_stats.csv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/annotation" ]; then mkdir {sample}/tmp/annotation; fi
	if [ ! -d "$(pwd)/{sample}/tmp/annotation/tmp" ]; then mkdir {sample}/tmp/annotation/tmp; fi
	splits=$(({proc} / 3)) && if [ $splits -gt {bakta_split} ]; then splits={bakta_split}; fi;
	find {sample}/tmp/binning/bins -type f -name "*.fa" -printf '%f\n' | sed 's/...$//' | xargs -i --max-procs=$splits bash -c 'bakta --compliant --db {db_bakta} --prefix {{}} --output {sample}/tmp/annotation/bakta/{{}} --keep-contig-headers --tmp-dir {sample}/tmp/annotation/tmp --threads 3 {sample}/tmp/binning/bins/{{}}.fa --skip-crispr'

	echo "bin,bakta_CDS_all,bakta_CDS_hyp,bakta_tRNA_all,bakta_tRNA_uniq,bakta_16S,bakta_23S,bakta_5S" > {output.df}
	for file in {sample}/tmp/annotation/bakta/*; do
	name=$(basename $file )
	CDS_all=$(awk -F "\t" '{{ if ($2 == "cds") {{print $2}} }}' $file/${{name}}.tsv | grep -c "cds" -) || true
	CDS_hyp=$(awk -F "\t" '{{ if ($2 == "cds") {{print $8}} }}' $file/${{name}}.tsv | grep -c "hypothetical protein" -) || true
	tRNA_all=$(awk -F "\t" '{{ if ($2 == "tRNA") {{print $7}} }}' $file/${{name}}.tsv | sed 's/fMet_trna/Met_trna/' - | sed 's/Ile2_trna/Ile_trna/' - | sed 's/SeC_trna/Cys_trna/' - | grep -c "trna" -) || true
	tRNA_uniq=$(awk -F "\t" '{{ if ($2 == "tRNA") {{print $7}} }}' $file/${{name}}.tsv | sed 's/fMet_trna/Met_trna/' - | sed 's/Ile2_trna/Ile_trna/' - | sed 's/SeC_trna/Cys_trna/' - | sort -u - | grep -c "trna" -) || true
	rRNA_16S=$(awk -F "\t" '{{ if ($7 == "16S_rrna") {{print $2}} }}' $file/${{name}}.tsv | grep -c "rRNA" -) || true
	rRNA_23S=$(awk -F "\t" '{{ if ($7 == "23S_rrna") {{print $2}} }}' $file/${{name}}.tsv | grep -c "rRNA" -) || true
	rRNA_5S=$(awk -F "\t" '{{ if ($7 == "5S_rrna") {{print $2}} }}' $file/${{name}}.tsv | grep -c "rRNA" -) || true
	echo "$name,$CDS_all,$CDS_hyp,$tRNA_all,$tRNA_uniq,$rRNA_16S,$rRNA_23S,$rRNA_5S" >> {output.df}; done
	
	cp -avr {sample}/tmp/annotation/bakta {sample}/results/.
        """

rule Taxonomy_aggregate:
    conda: "env_0"
    input:
        expand("{sample}/tmp/taxonomy/midas_16s.tsv",sample=sample),
        expand("{sample}/tmp/taxonomy/silva_16s.tsv",sample=sample),
        expand("{sample}/tmp/taxonomy/kaiju.out.names",sample=sample),
        expand("{sample}/tmp/taxonomy/gtdbtk.tsv",sample=sample)
    output:
        contigs=expand("{sample}/tmp/taxonomy/contigs_taxonomy.tsv",sample=sample),
        bins=expand("{sample}/tmp/taxonomy/bins_taxonomy.tsv",sample=sample)
    shell:
        """
	R --slave --silent --args << 'df'
	# Load and wrangle
	assembly_info=read.delim("{sample}/tmp/flye/assembly_info.txt", sep="\t", header=T)
	assembly_info=assembly_info[,c(1,2,3,4,5,8)]
	colnames(assembly_info) <- c("contig","len_bp","cov","status_circular","status_repeat","graph_path")
	
	kaiju=read.delim("{sample}/tmp/taxonomy/kaiju.out.names", sep="\t", header=F)
	kaiju=kaiju[(kaiju$V1 == "C"), ]
	kaiju=kaiju[c("V2","V4")]
	colnames(kaiju) <- c("contig","tax_kaiju")
	kaiju=as.data.frame(kaiju, stringsAsFactors = FALSE)
	
	silva <- read.delim("{sample}/tmp/taxonomy/silva_16s.tsv", sep="\t", header=F)
	silva=cbind(silva, data.frame(do.call('rbind', strsplit(as.character(silva$V1),':',fixed=TRUE))))
	silva=silva[c("X3","V2","V3")]
	colnames(silva) <- c("contig","tax_silva","tophit_silva")

	midas <- read.delim("{sample}/tmp/taxonomy/midas_16s.tsv", sep="\t", header=F)
	midas=cbind(midas, data.frame(do.call('rbind', strsplit(as.character(midas$V1),':',fixed=TRUE))))
	midas=midas[c("X3","V2","V3")]
	colnames(midas) <- c("contig","tax_midas","tophit_midas")

	links <- read.delim("{sample}/tmp/binning/contig_bin.tsv", sep="\t", header=F)
	colnames(links) <- c("contig","bin")

	gtdb <- read.delim("{sample}/tmp/taxonomy/gtdbtk.tsv", sep="\t", header=T)
	gtdb <- gtdb[,c(1,2,6,7,8,11,12,17,19,20)]
	colnames(gtdb)[1] <- "bin"
	colnames(gtdb)[2] <- "tax_gtdb"
	colnames(gtdb)[10] <- "gtdb_warning"
	
	# Silva taxonomy for contigs
	silva_tmp <- aggregate(silva$contig, by=list(silva$contig, silva$tax_silva), FUN=length)
	silva_tmp$sub <- paste(silva_tmp$Group.1,silva_tmp$Group.2,sep="_")
	silva_tmp <- do.call(rbind, lapply(split(silva_tmp,silva_tmp$Group.1), function(x) {{return(x[which.max(x$x),])}}))
	colnames(silva_tmp) <- c("contig","tax_silva","n","id")

	silva_tmp2 <- silva
	silva_tmp2$id <- paste(silva_tmp2$contig,silva_tmp2$tax_silva,sep="_")
	silva_tmp2 <- silva_tmp2[c("id","tophit_silva")]
	silva_tmp2 <- do.call(rbind, lapply(split(silva_tmp2,silva_tmp2$id), function(x) {{return(x[which.max(x$tophit_silva),])}}))
	
	silva_tmp <- merge(silva_tmp,silva_tmp2, by="id")
	silva_tmp$id <- NULL
	silva_tmp$n <- NULL
	
	# Midas taxonomy for contigs
	midas_tmp <- aggregate(midas$contig, by=list(midas$contig, midas$tax_midas), FUN=length)
	midas_tmp$sub <- paste(midas_tmp$Group.1,midas_tmp$Group.2,sep="_")
	midas_tmp <- do.call(rbind, lapply(split(midas_tmp,midas_tmp$Group.1), function(x) {{return(x[which.max(x$x),])}}))
	colnames(midas_tmp) <- c("contig","tax_midas","n","id")

	midas_tmp2 <- midas
	midas_tmp2$id <- paste(midas_tmp2$contig,midas_tmp2$tax_midas,sep="_")
	midas_tmp2 <- midas_tmp2[c("id","tophit_midas")]
	midas_tmp2 <- do.call(rbind, lapply(split(midas_tmp2,midas_tmp2$id), function(x) {{return(x[which.max(x$tophit_midas),])}}))
	
	midas_tmp <- merge(midas_tmp,midas_tmp2, by="id")
	midas_tmp$id <- NULL
	midas_tmp$n <- NULL
	
	# Make contigs dataframe
	contigs <- merge(assembly_info,kaiju, by="contig")
	contigs <- merge(contigs,silva_tmp, by="contig", all=TRUE)
	contigs <- merge(contigs,midas_tmp, by="contig", all=TRUE)
	contigs <- merge(contigs,links, by="contig", all=TRUE)
	write.table(contigs,"{output.contigs}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
	
	# Silva taxonomy for bins
	silva_tmp <- merge(silva,links, by="contig")
	silva_tmp$bin <- as.factor(silva_tmp$bin)
	
	silva_tmp2 <- aggregate(silva_tmp$bin, by=list(silva_tmp$bin, silva_tmp$tax_silva), FUN=length)
	silva_tmp2$sub <- paste(silva_tmp2$Group.1,silva_tmp2$Group.2,sep="_")
	silva_tmp2 <- do.call(rbind, lapply(split(silva_tmp2,silva_tmp2$Group.1), function(x) {{return(x[which.max(x$x),])}}))
	colnames(silva_tmp2) <- c("bin","tax_silva","n","id")
	
	silva_tmp3 <- silva_tmp
	silva_tmp3$id <- paste(silva_tmp3$bin,silva_tmp3$tax_silva,sep="_")
	silva_tmp3 <- silva_tmp3[c("id","tophit_silva")]
	silva_tmp3 <- do.call(rbind, lapply(split(silva_tmp3,silva_tmp3$id), function(x) {{return(x[which.max(x$tophit_silva),])}}))
	
	silva_tmp <- merge(silva_tmp2,silva_tmp3, by="id")
	silva_tmp$id <- NULL
	silva_tmp$n <- NULL
	
	# Midas taxonomy for bins
	midas_tmp <- merge(midas,links, by="contig")
	midas_tmp$bin <- as.factor(midas_tmp$bin)
	
	midas_tmp2 <- aggregate(midas_tmp$bin, by=list(midas_tmp$bin, midas_tmp$tax_midas), FUN=length)
	midas_tmp2$sub <- paste(midas_tmp2$Group.1,midas_tmp2$Group.2,sep="_")
	midas_tmp2 <- do.call(rbind, lapply(split(midas_tmp2,midas_tmp2$Group.1), function(x) {{return(x[which.max(x$x),])}}))
	colnames(midas_tmp2) <- c("bin","tax_midas","n","id")
	
	midas_tmp3 <- midas_tmp
	midas_tmp3$id <- paste(midas_tmp3$bin,midas_tmp3$tax_midas,sep="_")
	midas_tmp3 <- midas_tmp3[c("id","tophit_midas")]
	midas_tmp3 <- do.call(rbind, lapply(split(midas_tmp3,midas_tmp3$id), function(x) {{return(x[which.max(x$tophit_midas),])}}))
	
	midas_tmp <- merge(midas_tmp2,midas_tmp3, by="id")
	midas_tmp$id <- NULL
	midas_tmp$n <- NULL
	
	# Make bins dataframe
	bins <- merge(gtdb,silva_tmp, by="bin", all=TRUE)
	bins <- merge(bins,midas_tmp, by="bin", all=TRUE)
	write.table(bins,"{output.bins}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Taxonomy_MAGs:
    conda: "env_2"
    input:
        expand("{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample)
    output:
        expand("{sample}/tmp/taxonomy/gtdbtk.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/taxonomy" ]; then mkdir {sample}/tmp/taxonomy; fi
	export GTDBTK_DATA_PATH="{db_gtdb}"
	gtdbtk classify_wf --cpus {proc} --genome_dir {sample}/tmp/binning/bins --extension .fa --out_dir {sample}/tmp/taxonomy/gtdb --tmpdir {sample}/tmp/taxonomy 
	cp {sample}/tmp/taxonomy/gtdb/classify/gtdbtk.bac120.summary.tsv {output}
	if [ -f "$(pwd)/{sample}/tmp/taxonomy/gtdb/classify/gtdbtk.ar53.summary.tsv" ]; then tail -n+2 {sample}/tmp/taxonomy/gtdb/classify/gtdbtk.ar53.summary.tsv >> {output}; fi
        """

rule Taxonomy_contigs:
    conda: "env_1"
    input:
        expand("{sample}/tmp/polishing/asm_pol_lenfilt.fasta",sample=sample)
    output:
        res=expand("{sample}/tmp/taxonomy/kaiju.out",sample=sample),
        names=expand("{sample}/tmp/taxonomy/kaiju.out.names",sample=sample),
        sum=expand("{sample}/tmp/taxonomy/kaiju_summary_phylum.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/taxonomy" ]; then mkdir {sample}/tmp/taxonomy; fi
	kaiju -z {proc} -t {db_kaiju}/nodes.dmp -f {db_kaiju}/kaiju_db.fmi -i {input} -o {output.res}
	kaiju-addTaxonNames -r superkingdom,phylum,class,order,family,genus,species -t {db_kaiju}/nodes.dmp -n {db_kaiju}/names.dmp -i {output.res} -o {output.names}
	kaiju2table -t {db_kaiju}/nodes.dmp -n {db_kaiju}/names.dmp -r phylum -o {output.sum} {output.res}
        """

rule Taxonomy_16S:
    conda: "env_1"
    input:
        expand("{sample}/tmp/polishing/asm_pol_lenfilt.fasta",sample=sample)
    output:
        rRNA=expand("{sample}/tmp/taxonomy/rRNA.fa",sample=sample),
        SSU=expand("{sample}/tmp/taxonomy/rRNA_16S.fa",sample=sample),
        midas=expand("{sample}/tmp/taxonomy/midas_16s.tsv",sample=sample),
        silva=expand("{sample}/tmp/taxonomy/silva_16s.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/taxonomy" ]; then mkdir {sample}/tmp/taxonomy; fi
	barrnap {input} --threads {proc} --outseq {output.rRNA}
	grep -A1 ">16S" {output.rRNA} > {output.SSU}
	cp {output.rRNA} {sample}/results/.
	cp {output.SSU} {sample}/results/.
	vsearch -usearch_global {output.SSU} -db {db_midas} --threads {proc} --minseqlength 1000 -strand both -id 0.70 -top_hits_only -maxaccepts 1 -maxrejects 0 -blast6out {output.midas}
	vsearch -usearch_global {output.SSU} -db {db_silva} --threads {proc} --minseqlength 1000 -strand both -id 0.70 -top_hits_only -maxaccepts 1 -maxrejects 0 -blast6out {output.silva}
        """
