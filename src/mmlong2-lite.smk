# DESCRIPTION: Snakemake workflow for recovering metagenome assembled genomes (MAGs) with long reads from PacBio or Nanopore sequencing
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

import os
import sys
import glob
import pandas as pd
import numpy as np

shell.executable("/bin/bash")
singularity: config["sing"]

wf_v=config["version"]
mode=config["mode"]
modes=['Nanopore-simplex', 'PacBio-HiFi']

proc=config["proc"]
proc_sub=config["proc_sub"]

loc=config["loc"]
sample=config["sample"]

fastq=config["fastq"]
reads_diffcov=config["reads_diffcov"]

semibin_mod=config["semibin_mod"]
semibin_mods=['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'chicken_caecum', 'global']

semibin_prot=config["semibin_prot"]
semibin_prots=['prodigal', 'fraggenescan', 'fast-naive']

min_compl_n={'1': config["min_compl_1"],'2': config["min_compl_2"],'3': config["min_compl_3"],'4': config["min_compl_4"]}
min_cont_n={'1': config["min_cont_1"],'2': config["min_cont_2"],'3': config["min_cont_3"],'4': config["min_cont_4"]}
das_tool_score_n={'1': config["das_tool_score_1"],'2': config["das_tool_score_2"],'3': config["das_tool_score_3"]}

def get_reads(reads_diffcov,fastq,mode):
    reads_diffcov_main = pd.DataFrame([get_read_type(mode),fastq]).T

    if reads_diffcov != "none":
        reads_diffcov_sup = pd.read_csv(reads_diffcov, index_col=None, header=None)
        reads_diffcov_main = pd.concat([reads_diffcov_main,reads_diffcov_sup])

    reads_diffcov_main.columns = ['type', 'reads']
    reads_diffcov_main['pos'] = np.arange(1, len(reads_diffcov_main) + 1).astype(str)
    reads_diffcov_main['id'] = reads_diffcov_main[["pos","type"]].agg("_".join, axis=1)
    reads = reads_diffcov_main['id'].tolist()
    return reads

def get_read_type(mode):
    if mode == "PacBio-HiFi": return "PB"
    else: return "NP"

def get_splits(proc,medaka_split):
    if int(proc/4) < medaka_split: return int(proc/4)
    else: return medaka_split

def get_contigs():	
    contigs_all=expand("{loc}/{sample}/tmp/polishing/lin_contigs_id_{split}",sample=sample,loc=loc,split=[f"{i:02}" for i in range(1, get_splits(proc,config["medaka_split"]) + 1)])
    return contigs_all

def get_assembly(mode,mdbg_status,medaka_status):
    if mode == "PacBio-HiFi" and mdbg_status == "FALSE": return "assembly/assembly.fasta"
    elif mode == "PacBio-HiFi" and mdbg_status == "TRUE": return "assembly/contigs.fasta.gz"
    elif medaka_status == "TRUE": return "polishing/asm_pol.fasta"
    else: return "assembly/assembly.fasta"

def get_map_mode(reads,np_map_mode,pb_map_mode,il_map_mode):
    if "NP" in reads: return np_map_mode
    if "PB" in reads: return pb_map_mode
    if "IL" in reads: return il_map_mode

def get_map_ident(reads,np_map_ident,pb_map_ident,il_map_ident):
    if "NP" in reads: return np_map_ident
    if "PB" in reads: return pb_map_ident
    if "IL" in reads: return il_map_ident

def get_map_cut(reads):
    if reads.startswith("1_"): return "1,2,3,4,5"
    else: return "4,5"

def get_semibin_cov(loc,sample,reads_diffcov,round,semibin_mod):
    if reads_diffcov == "none": return ''.join([str(loc), "/", str(sample), "/tmp/binning/mapping/*.bam", " --environment ", str(semibin_mod), " "])
    if round == "3": return ''.join([str(loc), "/", str(sample), "/tmp/binning/mapping/1_*.bam", " --environment ", str(semibin_mod), " "])
    else: return ''.join([str(loc), "/", str(sample), "/tmp/binning/mapping/*.bam "])

def get_qc_in(loc,sample,round):
    if round == "2": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/comebin/comebin_res/comebin_res.tsv"])
    elif round == "4": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/metabat2/bins_metabat2"])
    else: return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/dastool/output_DASTool_scaffolds2bin.txt"])

def get_qc_dir(loc,sample,round):
    if round == "2": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/comebin/comebin_res/comebin_res_bins"])
    elif round == "4": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/metabat2"])
    else: return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/dastool/output_DASTool_bins"])

def get_vae_in(vamb_status):
    if vamb_status == "TRUE": return "vamb/clusters.tsv"
    else: return "graphmb/_best_contig2bin.tsv"

def get_vae_bins(vamb_status):
    if vamb_status == "TRUE": return "vamb/bins"
    else: return "graphmb/_bins"

onstart:
    from snakemake.utils import min_version
    min_version("8.0.0")
    if not os.path.exists(fastq): sys.exit(print("Read input (((",fastq,"))) not found. Aborting..."))
    if not os.path.isdir(loc): sys.exit(print("Provided path for output (((",loc,"))) not found. Aborting..."))
    if reads_diffcov != "none" and not os.path.exists(reads_diffcov): sys.exit(print("Dataframe for differential coverage binning (((",reads_diffcov,"))) not found. Aborting..."))
    if not mode in modes: sys.exit(print("Provided workflow mode (((",mode,"))) not recognised. Aborting..."))
    if not semibin_mod in semibin_mods: sys.exit(print("Provided model for SemiBin (((",semibin_mod,"))) not recognised. Aborting..."))
    if not semibin_prot in semibin_prots: sys.exit(print("Provided gene predictor for SemiBin (((",semibin_prot,"))) not recognised. Aborting..."))
    if len(os.path.join(loc, sample)) > 85: sys.exit(print("Path for provided output too long: (((",os.path.join(loc, sample),")))\nPlease re-run with different output location."))
    if not os.path.exists(os.path.join(loc, sample)): os.makedirs(os.path.join(loc, sample))
    if not os.path.exists(os.path.join(loc, sample, "results")): os.makedirs(os.path.join(loc, sample, "results"))
    if not os.path.exists(os.path.join(loc, sample, "tmp")): os.makedirs(os.path.join(loc, sample, "tmp"))
    if mode == "PacBio-HiFi": ruleorder: Assembly_metaFlye > Filtering_length > Polishing_stitch > Binning_prep_innit > Binning_prep_main
    else: ruleorder: Assembly_metaFlye > Polishing_prep > Polishing_consensus > Polishing_stitch > Filtering_length > Binning_prep_innit > Binning_prep_main
    if not os.path.exists(os.path.join(loc, sample, "tmp", "dep_mmlong2-lite.csv")): 
	    with open(os.path.join(loc, sample, "tmp", "dep_mmlong2-lite.csv"), 'w') as f:
		    f.write("dependency,version\n")

onsuccess:
    from datetime import datetime
    now = datetime.now()
    print("MAG recovery with version",wf_v,"completed at",now.strftime("%Y/%m/%d %H:%M:%S"))
    print("Thank you for using mmlong2-lite")

onerror:
    print("An error has occurred. Inspect Snakemake log files for troubleshooting.")

rule Finalise:
    conda: config["env_1"]
    input:
        expand("{loc}/{sample}/tmp/binning/contig_bin.tsv",sample=sample,loc=loc),
    output:
        df1=expand("{loc}/{sample}/tmp/binning/bins_mmlong2-lite.tsv",sample=sample,loc=loc),
        df2=expand("{loc}/{sample}/results/{sample}_bins.tsv",sample=sample,loc=loc)
    shell:
        """
        if [ -f "{loc}/{sample}/results/dependencies.csv" ]; then rm {loc}/{sample}/results/dependencies.csv; fi
        rsync {loc}/{sample}/tmp/dep_mmlong2-lite.csv {loc}/{sample}/results/dependencies.csv
		
        R --no-echo --silent --args << 'make_df'
        quast <- read.delim("{loc}/{sample}/tmp/binning/quast.tsv", sep="\t", header=T)
        colnames(quast) <- c("bin","gc","contig_n50","contig_n90","aun","n_per_100kb")
        abund <- read.delim("{loc}/{sample}/tmp/binning/bin_abund.tsv", sep="\t", header=T)
        colnames(abund) <- c("bin","r_abund")
        cov <- read.delim("{loc}/{sample}/tmp/binning/bin_cov.tsv", sep="\t", header=T)
        colnames(cov) <- c("bin","cov")
        checkm1 <- read.delim("{loc}/{sample}/tmp/binning/checkm1.tsv", sep="\t", header=T)
        checkm1 <- checkm1[, c("Bin.Id","Completeness","Contamination","Strain.heterogeneity")]
        colnames(checkm1) <- c("bin","completeness_checkm1","contamination_checkm1","strain_heterogeneity_checkm1")
        checkm2 <- read.delim("{loc}/{sample}/tmp/binning/checkm2.tsv", sep="\t", header=T)
        checkm2 <- checkm2[, c("Name","Completeness","Contamination","Coding_Density","Genome_Size","Total_Contigs")]
        colnames(checkm2) <- c("bin","completeness_checkm2","contamination_checkm2","coding_density","genome_size","contigs")
        checkm2$coding_density <- checkm2$coding_density * 100
        bins <- merge(checkm1,merge(checkm2,merge(quast,merge(cov,abund,by="bin"),by="bin"), by="bin"),by="bin")
        bins$wf_name <- "{sample}"
        bins$wf_mode <- "{mode}"
        bins$wf_v <- "{wf_v}"
        bins$wf_date <- Sys.Date()
        write.table(bins,"{output.df1}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        write.table(bins,"{output.df2}", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
        """

rule Assembly_metaFlye:
    conda: config["env_2"]
    input:
        fastq
    output:
        expand("{loc}/{sample}/tmp/assembly/assembly.fasta",sample=sample,loc=loc)
    params:
        flye_cov=config["flye_cov"],
        flye_ovlp=config["flye_ovlp"],
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        if [ {mode} == "Nanopore-simplex" ]; then flye_opt="--nano-hq"; fi
        if [ {mode} == "PacBio-HiFi" ]; then flye_opt="--read-error 0.01 --pacbio-hifi"; fi
        if [ {params.flye_ovlp} -eq 0 ]; then flye_ovlp=""; else flye_ovlp="--min-overlap {params.flye_ovlp}"; fi
        
        flye $flye_opt {input} --out-dir {loc}/{sample}/tmp/assembly --threads {threads} --meta $flye_ovlp --extra-params min_read_cov_cutoff={params.flye_cov}
        if ! grep -q "flye" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "flye " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """

rule Assembly_metaMDBG:
    conda: config["env_4"]
    input:
        fastq
    output:
        expand("{loc}/{sample}/tmp/assembly/contigs.fasta.gz",sample=sample,loc=loc)
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        metaMDBG asm {loc}/{sample}/tmp/assembly {input} -t {threads}
        
        zcat {output} | awk '/^>/ {{if (name) {{split(substr(name, 2), parts, "_"); print substr(name, 2) "\t" length(seq) "\t" parts[2] "\t" parts[3]; }} name=$0; seq=""; next}} {{seq = seq $0}}
        END {{if (name) {{split(substr(name, 2), parts, "_"); print substr(name, 2) "\t" length(seq) "\t" parts[2] "\t" parts[3]; }} }}' > {loc}/{sample}/tmp/assembly/assembly.tsv
        
        cat {loc}/{sample}/tmp/assembly/assembly.tsv | sed 's/x\t/\t/g' | awk 'BEGIN {{FS=OFS="\t"}} {{gsub("l", "N", $4); gsub("rc", "N", $4); gsub("c", "Y", $4)}} 1' - > {loc}/{sample}/tmp/assembly/assembly_info.txt
        if ! grep -q "metamdbg" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "metamdbg " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """

rule Polishing_prep:
    conda: config["env_3"]
    input:
        expand("{loc}/{sample}/tmp/assembly/assembly.fasta",sample=sample,loc=loc)
    output:
        expand("{loc}/{sample}/tmp/polishing/lin_contigs_id_{split}.txt",sample=sample,loc=loc,split=[f"{i:02}" for i in range(1, get_splits(proc,config["medaka_split"]) + 1)]),
        expand("{loc}/{sample}/tmp/polishing/calls_to_draft.bam",sample=sample,loc=loc)
    params:
        splits=get_splits(proc,config["medaka_split"]),
        ram_reads=config["minimap_ram"],
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/polishing" ]; then rm -r {loc}/{sample}/tmp/polishing; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fasta.map-ont.mmi" ]; then rm {loc}/{sample}/tmp/assembly/assembly.fasta.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fasta.fai" ]; then rm {loc}/{sample}/tmp/assembly/assembly.fasta.fai; fi
		mkdir {loc}/{sample}/tmp/polishing
        
        grep ">" {input} | cut -c 2- | awk '{{print $1}}' > {loc}/{sample}/tmp/polishing/ids.txt
        split -n l/{params.splits} --numeric-suffixes=1 --additional-suffix=.txt -d {loc}/{sample}/tmp/polishing/ids.txt {loc}/{sample}/tmp/polishing/contigs_id_
        for file in {loc}/{sample}/tmp/polishing/contigs_id_*.txt; do filename=lin_$(basename $file)
        awk 'BEGIN {{ ORS = " " }} {{ print }}' $file > {loc}/{sample}/tmp/polishing/$filename; done
        mini_align -I {params.ram_reads}G -i {fastq} -r {input} -m -p {loc}/{sample}/tmp/polishing/calls_to_draft -t {threads}
        """

rule Polishing_consensus:
    conda: config["env_3"]
    input:
        "{contigs}.txt"
    output:
        "{contigs}.hdf"
    params:
        medaka_model=config["medak_mod_pol"],	
    resources: usage=config["medaka_usage"]
    threads: proc_sub
    shell:
        """
        if [ -f "{output}" ]; then rm {output}; fi
        medaka consensus {loc}/{sample}/tmp/polishing/calls_to_draft.bam {output} --model {params.medaka_model} --batch 200 --threads {threads} --region $(cat {input})
        """	

rule Polishing_stitch:
    conda: config["env_3"]
    input:
        expand("{contigs}.hdf", contigs=get_contigs())
    output:
        expand("{loc}/{sample}/tmp/polishing/asm_pol.fasta",sample=sample,loc=loc)
    threads: proc
    shell:
        """
        medaka stitch {input} {loc}/{sample}/tmp/assembly/assembly.fasta {output} --threads {threads}
        if ! grep -q "medaka" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "medaka " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi        
        """			

rule Filtering_length:
    conda: config["env_2"]
    input:
        expand("{loc}/{sample}/tmp/{assembly}",sample=sample,loc=loc,assembly=get_assembly(mode,config["mdbg_status"],config["medaka_status"]))
    output:
        expand("{loc}/{sample}/tmp/filtering/asm_filt_len.fasta",sample=sample,loc=loc)
    params:
        len=config["min_contig_len"],
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/filtering" ]; then rm -r {loc}/{sample}/tmp/filtering; fi
        mkdir {loc}/{sample}/tmp/filtering
        seqkit seq -m {params.len} {input} > {output}
        rsync {output} {loc}/{sample}/results/{sample}_assembly.fasta
        if ! grep -q "seqkit" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "seqkit " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi        
        """	

rule Filtering_eukaryotes:
    conda: config["env_2"]
    input:
        expand("{loc}/{sample}/tmp/filtering/asm_filt_len.fasta",sample=sample,loc=loc)
    output:
        expand("{loc}/{sample}/tmp/filtering/asm_filt_len_euk.fasta",sample=sample,loc=loc)
    params:
        len=config["min_contig_len"],
    threads: proc
    shell:
        """
        tiara -i {input} -t {threads} -m {params.len} -o {loc}/{sample}/tmp/filtering/tiara.tsv
        cut -f1,2 {loc}/{sample}/tmp/filtering/tiara.tsv | grep -e "prokarya" -e "bacteria" -e "archaea" -e "unknown" - | cut -f1 | sort > {loc}/{sample}/tmp/filtering/contigs_filt_len_euk.txt
        seqkit grep -f {loc}/{sample}/tmp/filtering/contigs_filt_len_euk.txt {input} > {output}
        if ! grep -q "tiara" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "tiara " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi     
        """

rule Singletons_circ:
    conda: config["env_1"]
    input:
        expand("{loc}/{sample}/tmp/filtering/asm_filt_len_euk.fasta",sample=sample,loc=loc)
    output:
        expand("{loc}/{sample}/tmp/binning/singl/contig_c.tsv",sample=sample,loc=loc)
    params:
        mag_len=config["min_mag_len"],
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/binning" ]; then mkdir {loc}/{sample}/tmp/binning; fi
        if [ ! -d "{loc}/{sample}/tmp/binning/singl" ]; then mkdir {loc}/{sample}/tmp/binning/singl; fi
        if [ ! -d "{loc}/{sample}/tmp/binning/singl/innit" ]; then mkdir {loc}/{sample}/tmp/binning/singl/innit; fi
        
        awk -F "\t" '{{ if (($2 >= {params.mag_len}) && ($4 == "Y")) {{print $1 "\t" "{sample}.bin.c." ++i; next}} }}' {loc}/{sample}/tmp/assembly/assembly_info.txt > {output}
        if [ -f "{loc}/{sample}/tmp/binning/singl/innit/{sample}.bin.c.*.fa" ]; then rm {loc}/{sample}/tmp/binning/singl/innit/{sample}.bin.c.*.fa; fi
        if [ $(cat {output} | wc -l) -ge 1 ]; then cat {output} | xargs -i --max-procs=1 -n 2 bash -c 'samtools faidx {loc}/{sample}/tmp/filtering/asm_filt_len.fasta $0 > {loc}/{sample}/tmp/binning/singl/innit/$1.fa'; fi;
        
        find {loc}/{sample}/tmp/binning/singl/innit -name "*.fa" -type f -size -2k -delete
        if ! grep -q "samtools" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "samtools " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Singletons_lin:
    conda: config["env_1"]
    input:
        expand("{loc}/{sample}/tmp/binning/singl/contig_c.tsv",sample=sample,loc=loc)
    output:
        expand("{loc}/{sample}/tmp/binning/singl/contig_l.tsv",sample=sample,loc=loc)
    params:
        mag_len=config["min_smag_len"],
    shell:
        """
        awk -F "\t" '{{ if (($2 >= {params.mag_len}) && ($4 == "N")) {{print $1 "\t" "{sample}.bin.s." ++i; next}} }}' {loc}/{sample}/tmp/assembly/assembly_info.txt > {output}
        if [ -f "{loc}/{sample}/tmp/binning/singl/innit/{sample}.bin.s.*.fa" ]; then rm {loc}/{sample}/tmp/binning/singl/innit/{sample}.bin.s.*.fa; fi
        if [ $(cat {output} | wc -l) -ge 1 ]; then cat {output} | xargs -i --max-procs=1 -n 2 bash -c 'samtools faidx {loc}/{sample}/tmp/filtering/asm_filt_len.fasta $0 > {loc}/{sample}/tmp/binning/singl/innit/$1.fa'; fi;
        find {loc}/{sample}/tmp/binning/singl/innit -name "*.fa" -type f -size -2k -delete
        """

rule Singletons_qc:
    conda: config["env_7"]
    input:
        expand("{loc}/{sample}/tmp/binning/singl/contig_l.tsv",sample=sample,loc=loc)
    output:
        expand("{loc}/{sample}/tmp/binning/singl/checkm2/quality_report.tsv",sample=sample,loc=loc)
    params:
        compl_c=config["min_compl_circ"],
        compl_l=config["min_compl_lin"],
        cont=config["min_cont_singl"],
        apptainer=config["apptainer_status"],
    threads: proc
    shell:
        """
        if [ "{params.apptainer}" == "FALSE" ]; then if [ ! -f $CONDA_DEFAULT_ENV/CheckM2_database/uniref100.KO.1.dmnd ]; then checkm2 database --download --path $CONDA_DEFAULT_ENV; fi; fi 
        
        if [ -d "{loc}/{sample}/tmp/binning/singl/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/singl/checkm2; fi
        checkm2 predict -x .fa -i {loc}/{sample}/tmp/binning/singl/innit -o {loc}/{sample}/tmp/binning/singl/checkm2 -t {threads}
        if ! grep -q "checkm2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then echo "checkm2,$(checkm2 --version)" >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi 
        
        awk -F "\t" '{{ if ((($1 ~ /bin.c/) && ($2 >= {params.compl_c}) && ($3 <= {params.cont})) || (($1 ~ /bin.s/) && ($2 >= {params.compl_l}) && ($3 <= {params.cont}))) {{print $1}} }}' {output} > {loc}/{sample}/tmp/binning/singl/bins_keep.txt
        if [ -d "{loc}/{sample}/tmp/binning/singl/bins" ]; then rm -r {loc}/{sample}/tmp/binning/singl/bins; fi
        mkdir {loc}/{sample}/tmp/binning/singl/bins
		
        if [ $(cat {loc}/{sample}/tmp/binning/singl/bins_keep.txt | wc -l) -ge 1 ]; then
        cat {loc}/{sample}/tmp/binning/singl/bins_keep.txt | xargs -i --max-procs=1 bash -c 'rsync {loc}/{sample}/tmp/binning/singl/innit/{{}}.fa {loc}/{sample}/tmp/binning/singl/bins/.'
        grep ">" {loc}/{sample}/tmp/binning/singl/bins/*.fa | cut -f2 -d">" | sort > {loc}/{sample}/tmp/binning/singl/binned.txt; fi
        """

rule Coverage_prep:
    conda: config["env_2"]
    input:
        expand("{loc}/{sample}/tmp/binning/singl/checkm2/quality_report.tsv",sample=sample,loc=loc)
    output:
        reads=expand("{loc}/{sample}/tmp/binning/mapping/reads.csv",sample=sample,loc=loc),
        links=expand("{loc}/{sample}/tmp/binning/mapping/{reads}.lnk",sample=sample,loc=loc,reads=get_reads(reads_diffcov,fastq,mode)),			
    params:
        reads=fastq,
        type=get_read_type(mode),
        extra=reads_diffcov,
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/mapping" ]; then rm -r {loc}/{sample}/tmp/binning/mapping; fi
        mkdir {loc}/{sample}/tmp/binning/mapping
        
        printf "{params.type},{params.reads}\n" > {output.reads}
        if [ -f {params.extra} ]; then cat {params.extra} >> {output.reads}; fi
        
        count=0
        cat {output.reads} | while read line || [ -n "$line" ]; do
        count=$(($count + 1))
        type="$(echo $line | cut -f1 -d",")"
        reads="$(echo $line | cut -f2 -d",")"
        ln -s $reads {loc}/{sample}/tmp/binning/mapping/${{count}}_${{type}}.lnk; done
        """

rule Coverage_map:
    conda: config["env_2"]
    input:
        "{loc}/{sample}/tmp/binning/mapping/{reads}.lnk"
    output:
        bam="{loc}/{sample}/tmp/binning/mapping/{reads}.bam",
        cov="{loc}/{sample}/tmp/binning/mapping/{reads}.tsv",
        tr="{loc}/{sample}/tmp/binning/mapping/{reads}_tr.tsv",
    params:
        map=lambda wildcards: get_map_mode(wildcards.reads,config["minimap_np"],config["minimap_pb"],config["minimap_il"]),
        ram_ref=config["minimap_ref"],
        ram_reads=config["minimap_ram"],
        ident=lambda wildcards: get_map_ident(wildcards.reads,config["np_map_ident"],config["pb_map_ident"],config["il_map_ident"]),
        tr=lambda wildcards: get_map_cut(wildcards.reads),
    resources: usage=config["mapping_usage"]
    threads: proc * (config["mapping_usage"]/100)
    shell:
        """
        if [ -f "{output.bam}" ]; then rm {output.bam}; fi
        minimap2 -I {params.ram_reads}G -K {params.ram_ref}G -t {threads} -ax {params.map} {loc}/{sample}/tmp/filtering/asm_filt_len.fasta {input} | samtools view --threads $(({threads} / 2)) -Sb -F 2308 - | samtools sort --threads $(({threads} / 2)) --write-index -o {output.bam}##idx##{output.bam}.bai -
        jgi_summarize_bam_contig_depths {output.bam} --percentIdentity {params.ident} --outputDepth {output.cov}
        cut -f{params.tr} {output.cov} > {output.tr}
        """

rule Coverage_aggregate:
    conda: config["env_2"]
    input:
        expand("{loc}/{sample}/tmp/binning/mapping/{reads}_tr.tsv",sample=sample,loc=loc,reads=get_reads(reads_diffcov,fastq,mode))
    output:
        expand("{loc}/{sample}/tmp/binning/mapping/cov_all.tsv",sample=sample,loc=loc)
    shell:
        """
        paste -d "\t" $(ls -v {input}) > {output}
        if ! grep -q "minimap2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "minimap2 " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if ! grep -q "metabat2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "metabat2 " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """

rule Comebin_prep_innit:
    conda: config["env_2"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.txt"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/mapping/header.sam"
    params:
        round="{round}",
        type=get_read_type(mode)
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/mapping" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/mapping; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/mapping
        awk 'FNR==NR{{contigs[$1]; next}} /^@SQ/ {{contig=$2; sub(/^SN:/, "", contig); if (contig in contigs) print}}' "{input}" <(samtools view --threads {threads} -H "{loc}/{sample}/tmp/binning/mapping/1_{params.type}.bam") > {output}
        """

rule Comebin_prep_main:
    conda: config["env_2"]
    input:
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.txt",
        header="{loc}/{sample}/tmp/binning/round_{round}/mapping/header.sam",
        bam="{loc}/{sample}/tmp/binning/mapping/{reads}.bam"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/mapping/{reads}.bam"
    params:
        round="{round}",
    resources: usage=config["mapping_usage"]
    threads: proc * (config["mapping_usage"]/100)
    shell:
        """
        grep -wFf "{input.contigs}" <(samtools view --threads {threads} "{input.bam}") | cat {input.header} - | samtools view --threads {threads} -b -o {output} -
        """

rule Binning_prep_innit:
    conda: config["env_2"]
    input:
        expand("{loc}/{sample}/tmp/binning/mapping/cov_all.tsv",sample=sample,loc=loc)
    output:
        id1=expand("{loc}/{sample}/tmp/binning/round_1/contigs_w_singl.txt",sample=sample,loc=loc),
        id2=expand("{loc}/{sample}/tmp/binning/round_1/contigs_wo_singl.txt",sample=sample,loc=loc),
        id3=expand("{loc}/{sample}/tmp/binning/round_1/contigs.txt",sample=sample,loc=loc),
        contigs=expand("{loc}/{sample}/tmp/binning/round_1/contigs.fasta",sample=sample,loc=loc),
        cov=expand("{loc}/{sample}/tmp/binning/round_1/cov.tsv",sample=sample,loc=loc)
    params:
        skip_con_filter=config["vamb_status"],
        len=config["min_contig_len"],
        connections=config["max_contig_con"],
        inc_len=config["inc_contig_len"],
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_1" ]; then rm -r {loc}/{sample}/tmp/binning/round_1; fi
        mkdir {loc}/{sample}/tmp/binning/round_1
        
        if [ {params.skip_con_filter} == "TRUE" ]
        then cut -f1 {loc}/{sample}/tmp/assembly/assembly_info.txt | sort > {output.id1}
        else awk -F "\t" '{{ if (($2 >= {params.len})) {{print $1, $2, gsub(/,/,"",$8)}} }}' {loc}/{sample}/tmp/assembly/assembly_info.txt | awk -F " " '{{ if (($2 >= {params.inc_len}) || ($3 <= {params.connections})) {{print $1}} }}' | sed 1d | sort > {output.id1}; fi
        
        comm -23 <(sort {output.id1}) <(sort {loc}/{sample}/tmp/binning/singl/binned.txt) > {output.id2}
        seqkit grep -f {output.id2} {loc}/{sample}/tmp/filtering/asm_filt_len_euk.fasta | seqkit sort --by-name - > {output.contigs}
        grep ">" {output.contigs} | cut -c2- > {output.id3}
        
        head -n1 {input} > {output.cov}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id3} {input} | sort >> {output.cov}
        """

rule Binning_prep_main:
    conda: config["env_2"]
    input:
        lambda wildcards: ''.join([str(wildcards.loc), "/", str(wildcards.sample), "/tmp/binning/round_",str(int(wildcards.round)-1),"/checkm2/quality_report.tsv"])
    output:
        id="{loc}/{sample}/tmp/binning/round_{round}/contigs.txt",
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
        cov="{loc}/{sample}/tmp/binning/round_{round}/cov.tsv",
    params:
        round_curr="{round}",
        round_prev=lambda wildcards: str(int(wildcards.round)-1),
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round_curr}" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round_curr}; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round_curr}
        
        comm -23 <(sort {loc}/{sample}/tmp/binning/round_{params.round_prev}/contigs.txt) <(sort {loc}/{sample}/tmp/binning/round_{params.round_prev}/binned.txt) > {output.id}
        seqkit grep -f {output.id} {loc}/{sample}/tmp/binning/round_{params.round_prev}/contigs.fasta | seqkit sort --by-name - > {output.contigs}
        head -n1 {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov.tsv > {output.cov}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id} {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov.tsv | sort >> {output.cov}
        """

rule Binning_metabat2:
    conda: config["env_2"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/metabat2/bins_metabat2"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        mag_len=config["min_mag_len"],
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/metabat2" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/metabat2; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/metabat2
        metabat2 -i {input} -a {loc}/{sample}/tmp/binning/round_{params.round}/cov.tsv -o {output} -t {threads} -m {params.len} -s {params.mag_len} --saveCls --seed {params.seed}
        if ! grep -q "metabat2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "metabat2 " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi	
        """

rule Binning_graphmb:
    conda: config["env_6"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/graphmb/_best_contig2bin.tsv"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        mag_len=config["min_mag_len"],
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/graphmb" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/graphmb; fi
        if [ -f "{loc}/{sample}/features.tsv" ]; then rm {loc}/{sample}/features.tsv; fi
        rsync {loc}/{sample}/tmp/assembly/assembly_graph.gfa {loc}/{sample}/tmp/binning/round_{params.round}/.
        loc_main=$(pwd) && cd {loc}/{sample}
        graphmb --assembly {loc}/{sample}/tmp/binning/round_{params.round} --outdir {loc}/{sample}/tmp/binning/round_{params.round}/graphmb --assembly_name contigs.fasta --depth cov.tsv --contignodes --numcores {threads} --assembly_type flye --vamb --minbin {params.mag_len} --mincontig {params.len} --seed {params.seed}
        rm features.tsv && cd $loc_main
        if ! grep -q "graphmb" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "graphmb " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_vamb:
    conda: config["env_6"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/vamb/clusters.tsv"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        mag_len=config["min_mag_len"],
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/vamb" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/vamb; fi
        loc_main=$(pwd) && cd {loc}/{sample}
        vamb --fasta {input} --outdir {loc}/{sample}/tmp/binning/round_{params.round}/vamb -m {params.len} --minfasta {params.mag_len} -p {threads} --jgi {loc}/{sample}/tmp/binning/round_{params.round}/cov.tsv
        rm embs.tsv && cd $loc_main
        for bin in {loc}/{sample}/tmp/binning/round_{params.round}/vamb/bins/*.fna; do mv "$bin" "${{bin%.fna}}.fa"; done
        if ! grep -q "vamb" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "vamb " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi	
        """

rule Binning_semibin2:
    conda: config["env_4"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/semibin/bins_info.tsv"
    params:
        cov=lambda wildcards: get_semibin_cov(wildcards.loc,wildcards.sample,reads_diffcov,wildcards.round,semibin_mod),
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        mag_len=config["min_mag_len"],
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/semibin" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/semibin; fi
        SemiBin2 single_easy_bin -i {input} -b {params.cov} -o {loc}/{sample}/tmp/binning/round_{params.round}/semibin -p {threads} -m {params.len} --minfasta-kbs $(({params.mag_len}/1000)) --self-supervised --sequencing-type long_read --engine cpu --orf-finder {semibin_prot} --compression none --random-seed {params.seed}
        if ! grep -q "semibin" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "semibin " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_dastool:
    conda: config["env_2"]
    input:
        expand("{{loc}}/{{sample}}/tmp/binning/round_{{round}}/{vae}", vae=get_vae_in(config["vamb_status"])),
        "{loc}/{sample}/tmp/binning/round_{round}/metabat2/bins_metabat2",
        "{loc}/{sample}/tmp/binning/round_{round}/semibin/bins_info.tsv"
    output:
        metabat="{loc}/{sample}/tmp/binning/round_{round}/dastool/metabat2.tsv",
        semibin="{loc}/{sample}/tmp/binning/round_{round}/dastool/semibin.tsv",
        vae="{loc}/{sample}/tmp/binning/round_{round}/dastool/vae.tsv",
        dastool="{loc}/{sample}/tmp/binning/round_{round}/dastool/output_DASTool_scaffolds2bin.txt"
    params:
        round="{round}",
        score=lambda wildcards: das_tool_score_n.get(str(wildcards.round), None),
        vae_bins=get_vae_bins(config["vamb_status"]),
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/dastool" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/dastool; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/dastool
        
        Fasta_to_Scaffolds2Bin.sh -i {loc}/{sample}/tmp/binning/round_{params.round}/metabat2 -e fa > {output.metabat}
        Fasta_to_Scaffolds2Bin.sh -i {loc}/{sample}/tmp/binning/round_{params.round}/semibin/output_bins -e fa > {output.semibin}
        Fasta_to_Scaffolds2Bin.sh -i {loc}/{sample}/tmp/binning/round_{params.round}/{params.vae_bins} -e fa > {output.vae}
        
        DAS_Tool -i {output.metabat},{output.semibin},{output.vae} -l MetaBAT2,SemiBin2,VAE -c {loc}/{sample}/tmp/binning/round_{params.round}/contigs.fasta -o {loc}/{sample}/tmp/binning/round_{params.round}/dastool/output --score_threshold {params.score} -t {threads} --search_engine diamond --write_bins
        if ! grep -q "das_tool" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "das_tool " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_comebin:
    conda: config["env_5"]
    input:
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
        bams=expand("{{loc}}/{{sample}}/tmp/binning/round_{{round}}/mapping/{reads}.bam", reads=get_reads(reads_diffcov,fastq,mode))
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/comebin/comebin_res/comebin_res.tsv"
    params:
        round="{round}",
        len=config["min_contig_len"],
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/comebin" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/comebin; fi
        n_contigs=$(grep ">" {input.contigs} | wc -l)
        if [ $n_contigs -lt 1024 ]; then batch_size=$n_contigs; else batch_size=1024;fi
        run_comebin.sh -a {input.contigs} -o {loc}/{sample}/tmp/binning/round_{params.round}/comebin -t {threads} -p {loc}/{sample}/tmp/binning/round_{params.round}/mapping -b $batch_size
        if ! grep -q "comebin" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "comebin " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi	
        """

rule Binning_qc:
    conda: config["env_7"]
    input:
        lambda wildcards: get_qc_in(wildcards.loc,wildcards.sample,wildcards.round)
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/checkm2/quality_report.tsv"
    params:
        round="{round}",
        dir=lambda wildcards: get_qc_dir(wildcards.loc,wildcards.sample,wildcards.round),
        compl=lambda wildcards: min_compl_n.get(str(wildcards.round), None),
        cont=lambda wildcards: min_cont_n.get(str(wildcards.round), None),
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/bins" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/bins; fi
        
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit
        n=1 && for i in {params.dir}/*.fa; do rsync $i {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/{sample}.bin.{params.round}.${{n}}.fa && n=$(($n+1)); done
        checkm2 predict -x .fa -i {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit -o {loc}/{sample}/tmp/binning/round_{params.round}/checkm2 -t {threads}
        
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/bins
        awk -F "\t" '{{ if (($2 >= {params.compl}) && ($3 <= {params.cont})) {{print $1}} }}' {output} > {loc}/{sample}/tmp/binning/round_{params.round}/bins_keep.txt
        
        if [ $(cat {loc}/{sample}/tmp/binning/round_{params.round}/bins_keep.txt | wc -l) -ge 1 ]; then
        cat {loc}/{sample}/tmp/binning/round_{params.round}/bins_keep.txt | xargs -i --max-procs=1 bash -c 'rsync {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/{{}}.fa {loc}/{sample}/tmp/binning/round_{params.round}/bins/.'
        grep ">" {loc}/{sample}/tmp/binning/round_{params.round}/bins/*.fa | cut -f2 -d">" | sort > {loc}/{sample}/tmp/binning/round_{params.round}/binned.txt; fi
        """

rule Binning_aggregate:
    conda: config["env_2"]
    input:
        expand("{loc}/{sample}/tmp/binning/round_{round}/checkm2/quality_report.tsv",sample=sample,loc=loc,round=list(range(1, 5))),
    output:
        expand("{loc}/{sample}/tmp/binning/checkm2.tsv",sample=sample,loc=loc)
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/bins_innit; fi
        mkdir {loc}/{sample}/tmp/binning/bins_innit
        rsync {loc}/{sample}/tmp/binning/*/bins/*.fa {loc}/{sample}/tmp/binning/bins_innit/.
        
        head -n1 {loc}/{sample}/tmp/binning/round_1/checkm2/quality_report.tsv > {output}
        for stage in "singl" "round_1" "round_2" "round_3" "round_4" 
        do if [ $(cat {loc}/{sample}/tmp/binning/$stage/bins_keep.txt | wc -l) -ge 1 ]; then grep -w -f {loc}/{sample}/tmp/binning/$stage/bins_keep.txt {loc}/{sample}/tmp/binning/$stage/checkm2/quality_report.tsv >> {output}; fi; done
        """

rule Binning_qc2:
    conda: config["env_8"]
    input:
        expand("{loc}/{sample}/tmp/binning/checkm2.tsv",sample=sample,loc=loc)
    output:
        checkm1=expand("{loc}/{sample}/tmp/binning/checkm1.tsv",sample=sample,loc=loc),
        bins=expand("{loc}/{sample}/tmp/binning/bins_keep.txt",sample=sample,loc=loc),
    params:
        compl=config["min_compl_checkm1"],
        cont=config["min_cont_checkm1"],
        score=config["min_score_checkm1"],
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/checkm1" ]; then rm -r {loc}/{sample}/tmp/binning/checkm1; fi
        checkm lineage_wf -x fa -t {threads} -f {output.checkm1} --pplacer_threads $(({threads} / 2)) --reduced_tree --tab_table {loc}/{sample}/tmp/binning/bins_innit {loc}/{sample}/tmp/binning/checkm1
        if ! grep -q "checkm-" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "checkm-genome " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        
        if [ -d "{loc}/{sample}/tmp/binning/bins" ]; then rm -r {loc}/{sample}/tmp/binning/bins; fi
        mkdir {loc}/{sample}/tmp/binning/bins
        awk -F "\t" '{{ if (($12 >= {params.compl}) && ($13 <= {params.cont}) && ($12 - $13 * 5 > {params.score})) {{print $1}} }}' {output.checkm1} > {output.bins}
        cat {output.bins} | xargs -i --max-procs=1 bash -c 'rsync {loc}/{sample}/tmp/binning/bins_innit/{{}}.fa {loc}/{sample}/tmp/binning/bins/.'
        
        if [ -d "{loc}/{sample}/results/bins" ]; then rm -r {loc}/{sample}/results/bins; fi
        rsync -r {loc}/{sample}/tmp/binning/bins {loc}/{sample}/results/.
        """

rule Summary_coverage:
    conda: config["env_1"]
    input:
        expand("{loc}/{sample}/tmp/binning/checkm1.tsv",sample=sample,loc=loc),
    output:
        cov=expand("{loc}/{sample}/tmp/binning/bin_cov.tsv",sample=sample,loc=loc),
        abund=expand("{loc}/{sample}/tmp/binning/bin_abund.tsv",sample=sample,loc=loc),
    params:
        cov_method=config["cov_method"],
    threads: proc
    shell:
        """
        coverm genome -t {threads} -b {loc}/{sample}/tmp/binning/mapping/1_*.bam -d {loc}/{sample}/tmp/binning/bins -x fa -o {output.cov} -m {params.cov_method}
        coverm genome -t {threads} -b {loc}/{sample}/tmp/binning/mapping/1_*.bam -d {loc}/{sample}/tmp/binning/bins -x fa -o {output.abund} -m relative_abundance 
        if ! grep -q "coverm" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "coverm " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Summary_stats:
    conda: config["env_1"]
    input:
        expand("{loc}/{sample}/tmp/binning/bin_cov.tsv",sample=sample,loc=loc),
    output:
        links=expand("{loc}/{sample}/tmp/binning/contig_bin.tsv",sample=sample,loc=loc),
        stats=expand("{loc}/{sample}/tmp/binning/quast.tsv",sample=sample,loc=loc),
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/quast" ]; then rm -r {loc}/{sample}/tmp/binning/quast; fi
        grep ">" {loc}/{sample}/tmp/binning/bins/*.fa | xargs -n1 basename | sed "s/.fa:>/\t/" | awk -F'\t' '{{print $2 "\t" $1}}' > {output.links}
        quast.py {loc}/{sample}/tmp/binning/bins/*.fa -o {loc}/{sample}/tmp/binning/quast -t {threads} --no-plots --no-html --no-icarus
        cut -f1,17,18,19,20,23 {loc}/{sample}/tmp/binning/quast/transposed_report.tsv > {output.stats}
        if ! grep -q "quast" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "quast " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """
