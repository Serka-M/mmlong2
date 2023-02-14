# DESCRIPTION: Snakemake workflow for recovering metagenome assembled genomes with long reads
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

sing=config["sing"]
wf_v=config["version"]
mode=config["mode"]
modes=['Nanopore-simplex', 'PacBio-HiFi']
proc=config["proc"]
fastq=config["fastq"]
sample=config["sample"]
flye_cov=config["flye_cov"]
flye_ovlp=config["flye_ovlp"]
min_contig_len=config["min_contig_len"]
medaka_split=config["medaka_split"]
medak_mod_pol=config["medak_mod_pol"]
minimap_ram=config["minimap_ram"]
min_mag_len=config["min_mag_len"]
semibin_mod=config["semibin_mod"]
semibin_mods=['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'chicken_caecum', 'global']
semibin_prot=config["semibin_prot"]
semibin_prots=['prodigal', 'fraggenescan']
reads_diffcov=config["reads_diffcov"]
np_map_ident=config["np_map_ident"]
pb_map_ident=config["pb_map_ident"]
il_map_ident=config["il_map_ident"]
min_compl_1=config["min_compl_1"]
min_compl_2=config["min_compl_2"]
min_compl_3=config["min_compl_3"]
min_compl_4=config["min_compl_4"]
min_cont_1=config["min_cont_1"]
min_cont_2=config["min_cont_2"]
min_cont_3=config["min_cont_3"]
min_cont_4=config["min_cont_4"]
das_tool_score_1=config["das_tool_score_1"]
das_tool_score_2=config["das_tool_score_2"]
das_tool_score_3=config["das_tool_score_3"]

singularity: sing
shell.executable("/bin/bash")

onstart:
    from snakemake.utils import min_version
    import os
    import sys
    min_version("7.0.0")
    if not os.path.exists(fastq): sys.exit(print("Read input (((",fastq,"))) not found. Aborting..."))
    if reads_diffcov != "none" and not os.path.exists(reads_diffcov): sys.exit(print("Dataframe for differential coverage binning (((",reads_diffcov,"))) not found. Aborting..."))
    if not mode in modes: sys.exit(print("Provided workflow mode (((",mode,"))) not recognised. Aborting..."))
    if not semibin_mod in semibin_mods: sys.exit(print("Provided model for SemiBin (((",semibin_mod,"))) not recognised. Aborting..."))
    if not semibin_prot in semibin_prots: sys.exit(print("Provided gene predictor for SemiBin (((",semibin_prot,"))) not recognised. Aborting..."))

onsuccess:
    from datetime import datetime
    now = datetime.now()
    print("MAG recovery with version",wf_v,"completed at",now.strftime("%Y/%m/%d %H:%M:%S"))
    print("Thank you for using mmlong2-lite")

onerror:
    print("An error has occurred. Inspect Snakemake log files for troubleshooting.")

rule Finalise:
    conda: "env_0"
    input:
        expand("{sample}/tmp/binning/checkm2.tsv",sample=sample)
    output:
        dep=expand("{sample}/tmp/binning/dep_mmlong2-lite.csv",sample=sample),
        df1=expand("{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample),
        df2=expand("{sample}/results/{sample}_bins.tsv",sample=sample)
    shell:
        """
	printf "software,version\n" > {output.dep}
	printf "mmlong2-lite,{wf_v}\n" >> {output.dep}
	printf "Snakemake,7.19.1\n" >> {output.dep}
	printf "Singularity,3.8.6\n" >> {output.dep}
	printf "R,4.2.2\n" >> {output.dep}
	printf "Minimap2,2.24\n" >> {output.dep}
	printf "SAMtools,1.16.1\n" >> {output.dep}
	printf "SeqKit,2.3.1\n" >> {output.dep}
	printf "Flye,2.9.1\n" >> {output.dep}
	printf "Medaka,1.7.2\n" >> {output.dep}
	printf "Tiara,1.0.3\n" >> {output.dep}
	printf "MetaBAT2,2.15\n" >> {output.dep}
	printf "SemiBin,1.5\n" >> {output.dep}
	printf "GraphMB,0.1.5\n" >> {output.dep}
	printf "DAS_Tool,1.1.3\n" >> {output.dep}
	printf "CheckM2,1.0.0\n" >> {output.dep}
	printf "CoverM,0.6.1\n" >> {output.dep}
	printf "QUAST,5.2.0\n" >> {output.dep}
	cp {output.dep} {sample}/results/dependencies.csv
	cp -r {sample}/tmp/binning/bins {sample}/results/.
	
	quast.py {sample}/tmp/binning/bins/*.fa -o {sample}/tmp/binning/quast -t {proc}
	cut -f1,14,15,19,20,23 {sample}/tmp/binning/quast/transposed_report.tsv | sed 1d - > {sample}/tmp/binning/quast.tsv
	coverm genome -b {sample}/tmp/binning/mapping_tmp/1_cov.bam -d {sample}/tmp/binning/bins -x fa -m mean -o {sample}/tmp/binning/bin_cov.tsv
	coverm genome -b {sample}/tmp/binning/mapping_tmp/1_cov.bam -d {sample}/tmp/binning/bins -x fa -m relative_abundance -o {sample}/tmp/binning/bin_abund.tsv

	R --slave --silent --args << 'df'
	quast <- read.delim("{sample}/tmp/binning/quast.tsv", sep="\t", header=F)
	colnames(quast) <- c("bin","contigs","longest_contig","N90","auN","N_per_100kb")
	abund <- read.delim("{sample}/tmp/binning/bin_abund.tsv", sep="\t", header=T)
	colnames(abund) <- c("bin","r_abund")
	cov <- read.delim("{sample}/tmp/binning/bin_cov.tsv", sep="\t", header=T)
	colnames(cov) <- c("bin","cov")
	checkm2 <- read.delim("{sample}/tmp/binning/checkm2.tsv", sep="\t", header=T)
	checkm2$Translation_Table_Used <- NULL
	checkm2$Additional_Notes <- NULL
	colnames(checkm2)[1] <- "bin"
	bins <- merge(checkm2,merge(quast,merge(cov,abund,by="bin"),by="bin"), by="bin")
	bins$wf_name <- "{sample}"
	bins$wf_mode <- "{mode}"
	bins$wf_v <- "{wf_v}"
	bins$wf_date <- Sys.Date()
	write.table(bins,"{output.df1}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
	write.table(bins,"{output.df2}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Assembly:
    conda: "env_1"
    output:
        expand("{sample}/tmp/flye/assembly.fasta",sample=sample),
        expand("{sample}/tmp/flye/assembly_graph.gfa",sample=sample),
	expand("{sample}/tmp/flye/assembly_info.txt",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}" ]; then mkdir {sample}; fi
	if [ ! -d "$(pwd)/{sample}/results" ]; then mkdir {sample}/results; fi
	if [ ! -d "$(pwd)/{sample}/tmp" ]; then mkdir {sample}/tmp; fi
	if [ {mode} == "Nanopore-simplex" ]; then flye_opt="--nano-hq"; fi
	if [ {mode} == "PacBio-HiFi" ]; then flye_opt="--read-error 0.01 --pacbio-hifi"; fi
        if [ {flye_ovlp} -eq 0 ]; then flye_ovlp=""; else flye_ovlp="--min-overlap {flye_ovlp}"; fi
        flye $flye_opt {fastq} --out-dir {sample}/tmp/flye --threads {proc} --meta $flye_ovlp --extra-params min_read_cov_cutoff={flye_cov}     
        """

rule Polishing:
    conda: "env_1"
    input:
        expand("{sample}/tmp/flye/assembly.fasta",sample=sample)
    output:
        org=expand("{sample}/tmp/polishing/asm_pol.fasta",sample=sample),
        filt=expand("{sample}/tmp/polishing/asm_pol_lenfilt.fasta",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/polishing" ]; then mkdir {sample}/tmp/polishing; fi
	if [ {mode} == "PacBio-HiFi" ]; then cp {input} {output.org}; else
	grep ">" {input} | cut -c 2- | awk '{{print $1}}' > {sample}/tmp/polishing/ids.txt
	splits=$(({proc} / 2)) && if [ $splits -gt {medaka_split} ]; then splits={medaka_split}; fi;
	split -n l/$splits --numeric-suffixes=1 --additional-suffix=.txt -d {sample}/tmp/polishing/ids.txt {sample}/tmp/polishing/contigs_id_
	for file in {sample}/tmp/polishing/contigs_id_*.txt; do filename=lin_$(basename $file)
	awk 'BEGIN {{ ORS = " " }} {{ print }}' $file > {sample}/tmp/polishing/$filename; done
	mini_align -I {minimap_ram}G -i {fastq} -r {input} -m -p {sample}/tmp/polishing/calls_to_draft -t {proc}
	find {sample}/tmp/polishing -type f -name "lin_*" | xargs -i --max-procs=$splits bash -c 'list=$(head -1 {{}}) && file=$(basename {{}} | sed 's/\.txt//') && medaka consensus {sample}/tmp/polishing/calls_to_draft.bam {sample}/tmp/polishing/$file.hdf --model {medak_mod_pol} --batch 200 --threads 2 --region $list'
	medaka stitch {sample}/tmp/polishing/*.hdf {input} {output.org} --threads $splits; fi
	cp {output.org}  {sample}/results/{sample}_assembly.fasta
	seqkit seq -m {min_contig_len} {output.org} > {output.filt}
        """

rule Eukaryote_removal:
    conda: "env_1"
    input:
        expand("{sample}/tmp/polishing/asm_pol_lenfilt.fasta",sample=sample)
    output:
        expand("{sample}/tmp/eukfilt/asm_pol_lenfilt_eukfilt.fasta",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/eukfilt" ]; then mkdir {sample}/tmp/eukfilt; fi
	tiara -i {input} -t {proc} -m {min_contig_len} -o {sample}/tmp/eukfilt/tiara
	cut -f1,2 {sample}/tmp/eukfilt/tiara | grep -e "prokarya" -e "bacteria" -e "archaea" -e "unknown" - | cut -f1 | sort > {sample}/tmp/eukfilt/contigs_filt.txt
	seqkit grep -f {sample}/tmp/eukfilt/contigs_filt.txt {input} > {output}
        """

rule cMAG_extraction:
    conda: "env_1"
    input:
        expand("{sample}/tmp/polishing/asm_pol_lenfilt.fasta",sample=sample)
    output:
        df=expand("{sample}/tmp/binning/bins_c/contig_cbin.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning" ]; then mkdir {sample}/tmp/binning; fi
	if [ ! -d "$(pwd)/{sample}/tmp/binning/bins" ]; then mkdir {sample}/tmp/binning/bins; fi
	if [ ! -d "$(pwd)/{sample}/tmp/binning/bins_c" ]; then mkdir {sample}/tmp/binning/bins_c; fi
	if [ ! -d "$(pwd)/{sample}/tmp/binning/bins_c/cbins_init" ]; then mkdir {sample}/tmp/binning/bins_c/cbins_init; fi
	awk -F "\t" '{{ if (($2 >= {min_mag_len}) && ($4 == "Y")) {{print $1 "\t" "bin.c" ++i; next}} }}' {sample}/tmp/flye/assembly_info.txt > {output.df}
	if [ $(cat {output.df} | wc -l) -ge 1 ]; then cat {output.df} | xargs -i --max-procs=1 -n 2 bash -c 'samtools faidx {sample}/tmp/polishing/asm_pol_lenfilt.fasta $0 > {sample}/tmp/binning/bins_c/cbins_init/$1.fa'; fi;
	find {sample}/tmp/binning/bins_c/cbins_init -name "*.fa" -type f -size -2k -delete
        """

rule Coverage_calculation:
    conda: "env_1"
    input:
        expand("{sample}/tmp/polishing/asm_pol_lenfilt.fasta",sample=sample)
    output:
        reads=expand("{sample}/tmp/binning/reads.csv",sample=sample),
        cov=expand("{sample}/tmp/binning/metabat_cov.tsv",sample=sample),
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning" ]; then mkdir {sample}/tmp/binning; fi
	if [ ! -d "$(pwd)/{sample}/tmp/binning/mapping_tmp" ]; then mkdir {sample}/tmp/binning/mapping_tmp; fi
	if [ ! -d "$(pwd)/{sample}/tmp/binning/cov_tmp" ]; then mkdir {sample}/tmp/binning/cov_tmp; fi
	if [ {mode} == "PacBio-HiFi" ]; then reads="PB";fi
	if [ {mode} == "Nanopore-simplex" ]; then reads="NP";fi
	printf "$reads,{fastq}\n" > {output.reads}
	if [ -f {reads_diffcov} ]; then cat {reads_diffcov} >> {output.reads}; fi
	cov_init_=$(pwd)/{sample}/tmp/binning/cov_tmp

	count=0 && cat {output.reads} | while read line || [ -n "$line" ]; do
	count=$(($count + 1)) && type="$(echo $line | cut -f1 -d",")" && reads="$(echo $line | cut -f2 -d",")"
	if [ $type == "IL" ]; then mapping="sr" && ident={il_map_ident};
	elif [ $type == "NP" ]; then mapping="map-ont" && ident={np_map_ident};
	elif [ $type == "PB" ]; then mapping="map-hifi" && ident={pb_map_ident}; fi

	if [ ! -f $(pwd)/{sample}/tmp/binning/metabat_tmp/"${{count}}_metabat.tsv" ]; then 
	minimap2  -I {minimap_ram}G -K 10G -t {proc} -ax $mapping {sample}/tmp/polishing/asm_pol_lenfilt.fasta $reads | samtools view --threads $(({proc} / 2)) -Sb -F 2308 - | samtools sort --threads $(({proc} / 2)) - > $(pwd)/{sample}/tmp/binning/mapping_tmp/"${{count}}_cov.bam"
	jgi_summarize_bam_contig_depths $(pwd)/{sample}/tmp/binning/mapping_tmp/"${{count}}_cov.bam" --percentIdentity $ident --outputDepth $(pwd)/{sample}/tmp/binning/cov_tmp/"${{count}}_cov.tsv"	
	cov_init=$(pwd)/{sample}/tmp/binning/cov_tmp/${{count}}_cov.tsv
	cov_metabat=$(pwd)/{sample}/tmp/binning/cov_tmp/${{count}}_metabat.tsv
	if [ "$cov_init" == "$cov_init_/1_cov.tsv" ]; then cp $cov_init $cov_metabat; else cut -f4,5 $cov_init > $cov_metabat; fi; fi; done

	paste -d "\t" $(ls -v $cov_init_/*_metabat.tsv) > {output.cov}
        """

rule Binning_prep_1:
    conda: "env_1"
    params: 1
    input:
        expand("{sample}/tmp/eukfilt/asm_pol_lenfilt_eukfilt.fasta",sample=sample),
        expand("{sample}/tmp/binning/metabat_cov.tsv",sample=sample),
        expand("{sample}/tmp/binning/bins_c/contig_cbin.tsv",sample=sample)
    output:
        contigs=expand("{sample}/tmp/binning/round_1/contigs_lin.fasta",sample=sample),
        id1=expand("{sample}/tmp/binning/round_1/contigs_lin_unfilt.txt",sample=sample),
        id2=expand("{sample}/tmp/binning/round_1/contigs_lin.txt",sample=sample),
        cov=expand("{sample}/tmp/binning/round_1/metabat_cov_filt.tsv",sample=sample),
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning/round_{params}" ]; then mkdir {sample}/tmp/binning/round_{params}; fi
	awk -F "\t" '{{ if (($2 >= {min_contig_len}) && ($4 == "N")) {{print $1}} }}' {sample}/tmp/flye/assembly_info.txt | sort > {output.id1}
	seqkit grep -f {output.id1} {sample}/tmp/eukfilt/asm_pol_lenfilt_eukfilt.fasta | seqkit sort --by-name - > {output.contigs}
	grep ">" {output.contigs} | cut -c2- > {output.id2}
	head -n1 {sample}/tmp/binning/metabat_cov.tsv > {output.cov}
	cat {output.id2} | while read line || [ -n "$line" ]; do
	grep -w $line {sample}/tmp/binning/metabat_cov.tsv >> {output.cov}; done
        """

rule Binning_MetaBat2_1:
    conda: "env_1"
    params: 1
    input:
        expand("{sample}/tmp/binning/round_1/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_1/metabat2/bins_metabat2",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning/round_{params}/metabat2" ]; then mkdir {sample}/tmp/binning/round_{params}/metabat2; fi
	metabat2 -i {sample}/tmp/binning/round_{params}/contigs_lin.fasta -a {sample}/tmp/binning/round_{params}/metabat_cov_filt.tsv -o {output} -t {proc} -m {min_contig_len} -s {min_mag_len} --saveCls
        """

rule Binning_SemiBin_1:
    conda: "env_2"
    params: 1
    input:
        expand("{sample}/tmp/binning/round_1/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_1/semibin/bins_info.tsv",sample=sample)
    shell:
        """
	SemiBin2 single_easy_bin -i {input} -b {sample}/tmp/binning/mapping_tmp/*_cov.bam -o {sample}/tmp/binning/round_{params}/semibin -p {proc} -m {min_contig_len} --minfasta-kbs $(({min_mag_len}/1000)) --self-supervised --sequencing-type long_read --engine cpu --orf-finder {semibin_prot} --compression none
        """

rule Binning_GraphMB_1:
    conda: "env_3"
    params: 1
    input:
        expand("{sample}/tmp/binning/round_1/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_1/graphmb/_best_contig2bin.tsv",sample=sample)
    shell:
        """
	cp {sample}/tmp/flye/assembly_graph.gfa {sample}/tmp/binning/round_{params}/.
	cd {sample}
	graphmb --assembly tmp/binning/round_{params} --outdir tmp/binning/round_{params}/graphmb --assembly_name contigs_lin.fasta --depth metabat_cov_filt.tsv --contignodes --numcores {proc} --assembly_type flye --vamb --minbin {min_mag_len} --mincontig {min_contig_len}
	rm features.tsv
	cd ..
        """

rule Binning_DASTool_1:
    conda: "env_1"
    params: 1
    input:
        expand("{sample}/tmp/binning/round_1/metabat2/bins_metabat2",sample=sample),
        expand("{sample}/tmp/binning/round_1/semibin/bins_info.tsv",sample=sample),
        expand("{sample}/tmp/binning/round_1/graphmb/_best_contig2bin.tsv",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_1/das_tool/output_DASTool_scaffolds2bin.txt",sample=sample),
        metabat=expand("{sample}/tmp/binning/round_1/das_tool/metabat2.tsv",sample=sample),
        semibin=expand("{sample}/tmp/binning/round_1/das_tool/semibin.tsv",sample=sample),
        graphmb=expand("{sample}/tmp/binning/round_1/das_tool/graphmb.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning/round_{params}/das_tool" ]; then mkdir {sample}/tmp/binning/round_{params}/das_tool; fi
	Fasta_to_Scaffolds2Bin.sh -i {sample}/tmp/binning/round_{params}/metabat2 -e fa > {output.metabat}
	Fasta_to_Scaffolds2Bin.sh -i {sample}/tmp/binning/round_{params}/semibin/output_bins -e fa > {output.semibin}
	Fasta_to_Scaffolds2Bin.sh -i {sample}/tmp/binning/round_{params}/graphmb/_bins -e fa > {output.graphmb}
	DAS_Tool -i {output.metabat},{output.semibin},{output.graphmb} -l MetaBAT2,SemiBin,GraphMB -c {sample}/tmp/binning/round_{params}/contigs_lin.fasta -o {sample}/tmp/binning/round_{params}/das_tool/output -t {proc} --score_threshold {das_tool_score_1} --search_engine diamond --write_bins
        """

rule Binning_QC_1:
    conda: "env_4"
    params: 1 
    input:
        expand("{sample}/tmp/binning/round_1/das_tool/output_DASTool_scaffolds2bin.txt",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_2/contigs_lin.fasta",sample=sample),
        expand("{sample}/tmp/binning/round_2/metabat_cov_filt.tsv",sample=sample)
    shell:
        """
	next={sample}/tmp/binning/round_$(({params}+1))
	export prev={sample}/tmp/binning/round_{params}
	if [ ! -d "$(pwd)/$next" ]; then mkdir $next; fi
	if [ ! -d "$(pwd)/$prev/bins_innit" ]; then mkdir $prev/bins_innit; fi
	n=1 && for i in $prev/das_tool/output_DASTool_bins/*.fa; do cp $i $prev/bins_innit/{sample}.bin.{params}.${{n}}.fa && n=$(($n+1)); done
	/checkm2/bin/checkm2 predict -x .fa -i $prev/bins_innit -o $prev/checkm2 -t {proc}
	awk -F "\t" '{{ if (($2 >= {min_compl_1}) && ($3 <= {min_cont_1})) {{print $1}} }}' $prev/checkm2/quality_report.tsv > $prev/bins_keep.txt
	if [ $(cat $prev/bins_keep.txt | wc -l) -ge 1 ]; then
	cat $prev/bins_keep.txt | xargs -i --max-procs={proc} bash -c 'cp ${{prev}}/bins_innit/{{}}.fa {sample}/tmp/binning/bins/.'
	cat {sample}/tmp/binning/bins/{sample}.bin.{params}.*.fa > $prev/binned.fasta
	grep ">" $prev/binned.fasta | cut -c 2- | awk '{{print $1}}' | sort > $prev/binned.txt
	grep -v -w -f $prev/binned.txt $prev/contigs_lin.txt | sort > $prev/unbinned.txt
	seqkit grep -f $prev/unbinned.txt $prev/contigs_lin.fasta | seqkit sort --by-name - > $next/contigs_lin.fasta
	head -n1 $prev/metabat_cov_filt.tsv > $next/metabat_cov_filt.tsv
	cat $prev/unbinned.txt | while read line || [ -n "$line" ]; do
	grep -w $line $prev/metabat_cov_filt.tsv >> $next/metabat_cov_filt.tsv; done
	else cp $prev/contigs_lin.fasta $next/contigs_lin.fasta && cp $prev/metabat_cov_filt.tsv $next/metabat_cov_filt.tsv; fi
	grep ">" $next/contigs_lin.fasta | cut -c 2- | awk '{{print $1}}' | sort > $next/contigs_lin.txt
        """

rule Binning_MetaBat2_2:
    conda: "env_1"
    params: 2
    input:
        expand("{sample}/tmp/binning/round_2/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_2/metabat2/bins_metabat2",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning/round_{params}/metabat2" ]; then mkdir {sample}/tmp/binning/round_{params}/metabat2; fi
	metabat2 -i {sample}/tmp/binning/round_{params}/contigs_lin.fasta -a {sample}/tmp/binning/round_{params}/metabat_cov_filt.tsv -o {output} -t {proc} -m {min_contig_len} -s {min_mag_len} --saveCls
        """

rule Binning_SemiBin_2:
    conda: "env_2"
    params: 2
    input:
        expand("{sample}/tmp/binning/round_2/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_2/semibin/bins_info.tsv",sample=sample)
    shell:
        """
	SemiBin2 single_easy_bin -i {input} -b {sample}/tmp/binning/mapping_tmp/*_cov.bam -o {sample}/tmp/binning/round_{params}/semibin -p {proc} -m {min_contig_len} --minfasta-kbs $(({min_mag_len}/1000)) --self-supervised --sequencing-type long_read --engine cpu --orf-finder {semibin_prot} --compression none
        """

rule Binning_GraphMB_2:
    conda: "env_3"
    params: 2
    input:
        expand("{sample}/tmp/binning/round_2/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_2/graphmb/_best_contig2bin.tsv",sample=sample)
    shell:
        """
	cp {sample}/tmp/flye/assembly_graph.gfa {sample}/tmp/binning/round_{params}/.
	cd {sample}
	graphmb --assembly tmp/binning/round_{params} --outdir tmp/binning/round_{params}/graphmb --assembly_name contigs_lin.fasta --depth metabat_cov_filt.tsv --contignodes --numcores {proc} --assembly_type flye --vamb --minbin {min_mag_len} --mincontig {min_contig_len}
	rm features.tsv
	cd ..
        """

rule Binning_DASTool_2:
    conda: "env_1"
    params: 2
    input:
        expand("{sample}/tmp/binning/round_2/metabat2/bins_metabat2",sample=sample),
        expand("{sample}/tmp/binning/round_2/semibin/bins_info.tsv",sample=sample),
        expand("{sample}/tmp/binning/round_2/graphmb/_best_contig2bin.tsv",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_2/das_tool/output_DASTool_scaffolds2bin.txt",sample=sample),
        metabat=expand("{sample}/tmp/binning/round_2/das_tool/metabat2.tsv",sample=sample),
        semibin=expand("{sample}/tmp/binning/round_2/das_tool/semibin.tsv",sample=sample),
        graphmb=expand("{sample}/tmp/binning/round_2/das_tool/graphmb.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning/round_{params}/das_tool" ]; then mkdir {sample}/tmp/binning/round_{params}/das_tool; fi
	Fasta_to_Scaffolds2Bin.sh -i {sample}/tmp/binning/round_{params}/metabat2 -e fa > {output.metabat}
	Fasta_to_Scaffolds2Bin.sh -i {sample}/tmp/binning/round_{params}/semibin/output_bins -e fa > {output.semibin}
	Fasta_to_Scaffolds2Bin.sh -i {sample}/tmp/binning/round_{params}/graphmb/_bins -e fa > {output.graphmb}
	DAS_Tool -i {output.metabat},{output.semibin},{output.graphmb} -l MetaBAT2,SemiBin,GraphMB -c {sample}/tmp/binning/round_{params}/contigs_lin.fasta -o {sample}/tmp/binning/round_{params}/das_tool/output -t {proc} --score_threshold {das_tool_score_2} --search_engine diamond --write_bins 
        """

rule Binning_QC_2:
    conda: "env_4"
    params: 2 
    input:
        expand("{sample}/tmp/binning/round_2/das_tool/output_DASTool_scaffolds2bin.txt",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_3/contigs_lin.fasta",sample=sample),
        expand("{sample}/tmp/binning/round_3/metabat_cov_filt.tsv",sample=sample)
    shell:
        """
	next={sample}/tmp/binning/round_$(({params}+1))
	export prev={sample}/tmp/binning/round_{params}
	if [ ! -d "$(pwd)/$next" ]; then mkdir $next; fi
	if [ ! -d "$(pwd)/$prev/bins_innit" ]; then mkdir $prev/bins_innit; fi
	n=1 && for i in $prev/das_tool/output_DASTool_bins/*.fa; do cp $i $prev/bins_innit/{sample}.bin.{params}.${{n}}.fa && n=$(($n+1)); done
	/checkm2/bin/checkm2 predict -x .fa -i $prev/bins_innit -o $prev/checkm2 -t {proc}
	awk -F "\t" '{{ if (($2 >= {min_compl_2}) && ($3 <= {min_cont_2})) {{print $1}} }}' $prev/checkm2/quality_report.tsv > $prev/bins_keep.txt
	if [ $(cat $prev/bins_keep.txt | wc -l) -ge 1 ]; then
	cat $prev/bins_keep.txt | xargs -i --max-procs={proc} bash -c 'cp ${{prev}}/bins_innit/{{}}.fa {sample}/tmp/binning/bins/.'
	cat {sample}/tmp/binning/bins/{sample}.bin.{params}.*.fa > $prev/binned.fasta
	grep ">" $prev/binned.fasta | cut -c 2- | awk '{{print $1}}' | sort > $prev/binned.txt
	grep -v -w -f $prev/binned.txt $prev/contigs_lin.txt | sort > $prev/unbinned.txt
	seqkit grep -f $prev/unbinned.txt $prev/contigs_lin.fasta | seqkit sort --by-name - > $next/contigs_lin.fasta
	head -n1 $prev/metabat_cov_filt.tsv > $next/metabat_cov_filt.tsv
	cat $prev/unbinned.txt | while read line || [ -n "$line" ]; do
	grep -w $line $prev/metabat_cov_filt.tsv >> $next/metabat_cov_filt.tsv; done
	else cp $prev/contigs_lin.fasta $next/contigs_lin.fasta && cp $prev/metabat_cov_filt.tsv $next/metabat_cov_filt.tsv; fi
	grep ">" $next/contigs_lin.fasta | cut -c 2- | awk '{{print $1}}' | sort > $next/contigs_lin.txt
        """

rule Binning_MetaBat2_3:
    conda: "env_1"
    params: 3
    input:
        expand("{sample}/tmp/binning/round_3/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_3/metabat2/bins_metabat2",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning/round_{params}/metabat2" ]; then mkdir {sample}/tmp/binning/round_{params}/metabat2; fi
	metabat2 -i {sample}/tmp/binning/round_{params}/contigs_lin.fasta -a {sample}/tmp/binning/round_{params}/metabat_cov_filt.tsv -o {output} -t {proc} -m {min_contig_len} -s {min_mag_len} --saveCls
        """

rule Binning_SemiBin_3:
    conda: "env_2"
    params: 3
    input:
        expand("{sample}/tmp/binning/round_3/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_3/semibin/bins_info.tsv",sample=sample)
    shell:
        """
	SemiBin2 single_easy_bin -i {input} -b {sample}/tmp/binning/mapping_tmp/1_cov.bam -o {sample}/tmp/binning/round_{params}/semibin -p {proc} -m {min_contig_len} --minfasta-kbs $(({min_mag_len}/1000)) --self-supervised --sequencing-type long_read --engine cpu --orf-finder {semibin_prot} --compression none --environment {semibin_mod}
        """

rule Binning_GraphMB_3:
    conda: "env_3"
    params: 3
    input:
        expand("{sample}/tmp/binning/round_3/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_3/graphmb/_best_contig2bin.tsv",sample=sample)
    shell:
        """
	cp {sample}/tmp/flye/assembly_graph.gfa {sample}/tmp/binning/round_{params}/.
	cd {sample}
	graphmb --assembly tmp/binning/round_{params} --outdir tmp/binning/round_{params}/graphmb --assembly_name contigs_lin.fasta --depth metabat_cov_filt.tsv --contignodes --numcores {proc} --assembly_type flye --vamb --minbin {min_mag_len} --mincontig {min_contig_len}
	rm features.tsv
	cd ..
        """

rule Binning_DASTool_3:
    conda: "env_1"
    params: 3
    input:
        expand("{sample}/tmp/binning/round_3/metabat2/bins_metabat2",sample=sample),
        expand("{sample}/tmp/binning/round_3/semibin/bins_info.tsv",sample=sample),
        expand("{sample}/tmp/binning/round_3/graphmb/_best_contig2bin.tsv",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_3/das_tool/output_DASTool_scaffolds2bin.txt",sample=sample),
        metabat=expand("{sample}/tmp/binning/round_3/das_tool/metabat2.tsv",sample=sample),
        semibin=expand("{sample}/tmp/binning/round_3/das_tool/semibin.tsv",sample=sample),
        graphmb=expand("{sample}/tmp/binning/round_3/das_tool/graphmb.tsv",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning/round_{params}/das_tool" ]; then mkdir {sample}/tmp/binning/round_{params}/das_tool; fi
	Fasta_to_Scaffolds2Bin.sh -i {sample}/tmp/binning/round_{params}/metabat2 -e fa > {output.metabat}
	Fasta_to_Scaffolds2Bin.sh -i {sample}/tmp/binning/round_{params}/semibin/output_bins -e fa > {output.semibin}
	Fasta_to_Scaffolds2Bin.sh -i {sample}/tmp/binning/round_{params}/graphmb/_bins -e fa > {output.graphmb}
	DAS_Tool -i {output.metabat},{output.semibin},{output.graphmb} -l MetaBAT2,SemiBin,GraphMB -c {sample}/tmp/binning/round_{params}/contigs_lin.fasta -o {sample}/tmp/binning/round_{params}/das_tool/output -t {proc} --score_threshold {das_tool_score_3} --search_engine diamond --write_bins
        """

rule Binning_QC_3:
    conda: "env_4"
    params: 3
    input:
        expand("{sample}/tmp/binning/round_3/das_tool/output_DASTool_scaffolds2bin.txt",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_4/contigs_lin.fasta",sample=sample),
        expand("{sample}/tmp/binning/round_4/metabat_cov_filt.tsv",sample=sample)
    shell:
        """
	next={sample}/tmp/binning/round_$(({params}+1))
	export prev={sample}/tmp/binning/round_{params}
	if [ ! -d "$(pwd)/$next" ]; then mkdir $next; fi
	if [ ! -d "$(pwd)/$prev/bins_innit" ]; then mkdir $prev/bins_innit; fi
	n=1 && for i in $prev/das_tool/output_DASTool_bins/*.fa; do cp $i $prev/bins_innit/{sample}.bin.{params}.${{n}}.fa && n=$(($n+1)); done
	/checkm2/bin/checkm2 predict -x .fa -i $prev/bins_innit -o $prev/checkm2 -t {proc}
	awk -F "\t" '{{ if (($2 >= {min_compl_3}) && ($3 <= {min_cont_3})) {{print $1}} }}' $prev/checkm2/quality_report.tsv > $prev/bins_keep.txt
	if [ $(cat $prev/bins_keep.txt | wc -l) -ge 1 ]; then
	cat $prev/bins_keep.txt | xargs -i --max-procs={proc} bash -c 'cp ${{prev}}/bins_innit/{{}}.fa {sample}/tmp/binning/bins/.'
	cat {sample}/tmp/binning/bins/{sample}.bin.{params}.*.fa > $prev/binned.fasta
	grep ">" $prev/binned.fasta | cut -c 2- | awk '{{print $1}}' | sort > $prev/binned.txt
	grep -v -w -f $prev/binned.txt $prev/contigs_lin.txt | sort > $prev/unbinned.txt
	seqkit grep -f $prev/unbinned.txt $prev/contigs_lin.fasta | seqkit sort --by-name - > $next/contigs_lin.fasta
	head -n1 $prev/metabat_cov_filt.tsv > $next/metabat_cov_filt.tsv
	cat $prev/unbinned.txt | while read line || [ -n "$line" ]; do
	grep -w $line $prev/metabat_cov_filt.tsv >> $next/metabat_cov_filt.tsv; done
	else cp $prev/contigs_lin.fasta $next/contigs_lin.fasta && cp $prev/metabat_cov_filt.tsv $next/metabat_cov_filt.tsv; fi
	grep ">" $next/contigs_lin.fasta | cut -c 2- | awk '{{print $1}}' | sort > $next/contigs_lin.txt
        """

rule Binning_MetaBat2_4:
    conda: "env_1"
    params: 4
    input:
        expand("{sample}/tmp/binning/round_4/contigs_lin.fasta",sample=sample)
    output:
        expand("{sample}/tmp/binning/round_4/metabat2/bins_metabat2",sample=sample)
    shell:
        """
	if [ ! -d "$(pwd)/{sample}/tmp/binning/round_{params}/metabat2" ]; then mkdir {sample}/tmp/binning/round_{params}/metabat2; fi
	metabat2 -i {sample}/tmp/binning/round_{params}/contigs_lin.fasta -a {sample}/tmp/binning/round_{params}/metabat_cov_filt.tsv -o {output} -t {proc} -m {min_contig_len} -s {min_mag_len} --saveCls
        """

rule Binning_QC_4:
    conda: "env_4"
    params: 4
    input:
        expand("{sample}/tmp/binning/round_4/metabat2/bins_metabat2",sample=sample)
    output:
        expand("{sample}/tmp/binning/checkm2.tsv",sample=sample)
    shell:
        """
	export prev={sample}/tmp/binning/round_{params}
	if [ ! -d "$(pwd)/$prev/bins_innit" ]; then mkdir $prev/bins_innit; fi
	n=1 && for i in $prev/metabat2/*.fa; do cp $i $prev/bins_innit/{sample}.bin.{params}.${{n}}.fa && n=$(($n+1)); done
	n=1 && for i in {sample}/tmp/binning/bins_c/cbins_init/*.fa; do cp $i $prev/bins_innit/{sample}.bin.c.${{n}}.fa && n=$(($n+1)); done
	/checkm2/bin/checkm2 predict -x .fa -i $prev/bins_innit -o $prev/checkm2 -t {proc}
	awk -F "\t" '{{ if (($2 >= {min_compl_4}) && ($3 <= {min_cont_4})) {{print $1}} }}' $prev/checkm2/quality_report.tsv > $prev/bins_keep.txt
	if [ $(cat $prev/bins_keep.txt | wc -l) -ge 1 ]; then
	cat $prev/bins_keep.txt | xargs -i --max-procs={proc} bash -c 'cp ${{prev}}/bins_innit/{{}}.fa {sample}/tmp/binning/bins/.'; fi
	grep ">" {sample}/tmp/binning/bins/*.fa | xargs -L 1 basename | awk 'gsub(".fa:>", "\t")' | awk -F'\t' -v OFS='\t' '{{print $2,$1}}' > {sample}/tmp/binning/contig_bin.tsv
	ls {sample}/tmp/binning/bins/*.fa | xargs -L 1 basename | awk 'gsub(".fa", "")' > {sample}/tmp/binning/bins.txt
	head -n1 {sample}/tmp/binning/round_1/checkm2/quality_report.tsv > {output}
	if [ $(cat {sample}/tmp/binning/round_1/bins_keep.txt | wc -l) -ge 1 ]; then grep -w -f {sample}/tmp/binning/bins.txt {sample}/tmp/binning/round_1/checkm2/quality_report.tsv >> {output}; fi
	if [ $(cat {sample}/tmp/binning/round_2/bins_keep.txt | wc -l) -ge 1 ]; then grep -w -f {sample}/tmp/binning/bins.txt {sample}/tmp/binning/round_2/checkm2/quality_report.tsv >> {output}; fi
	if [ $(cat {sample}/tmp/binning/round_3/bins_keep.txt | wc -l) -ge 1 ]; then grep -w -f {sample}/tmp/binning/bins.txt {sample}/tmp/binning/round_3/checkm2/quality_report.tsv >> {output}; fi
	if [ $(cat {sample}/tmp/binning/round_4/bins_keep.txt | wc -l) -ge 1 ]; then grep -w -f {sample}/tmp/binning/bins.txt {sample}/tmp/binning/round_4/checkm2/quality_report.tsv >> {output}; fi
        """
