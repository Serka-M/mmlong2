# DESCRIPTION: Snakemake workflow for recovering metagenome assembled genomes (MAGs) with long reads from PacBio or Nanopore sequencing
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

import os
import sys
import glob
import re
import pandas as pd
import numpy as np

shell.executable("/bin/bash")
singularity: config["sing"]
wf_v=config["version"]

mode=config["mode"]
modes=['Nanopore-simplex', 'PacBio-HiFi']
assemblers=['myloasm', 'metaflye', 'metamdbg', 'custom']

binmode=config["binmode"]
binmodes=['fast', 'default', 'extended']

proc=config["proc"]
proc_sub=config["proc_sub"]

loc=config["loc"]
sample=config["sample"]

fastq=config["fastq"]
reads_diffcov=config["reads_diffcov"]

semibin_mod=config["semibin_mod"]
semibin_mods=['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'chicken_caecum', 'global']

semibin_prot=config["semibin_prot"]
semibin_prots=['prodigal', 'fast-naive']

def get_reads(reads_diffcov,fastq,mode):
    reads_diffcov_main = pd.DataFrame([get_read_type(mode),fastq]).T

    if reads_diffcov != "none":
        reads_diffcov_sup = pd.read_csv(reads_diffcov, index_col=None, header=None)
        reads_diffcov_main = pd.concat([reads_diffcov_main,reads_diffcov_sup])

    reads_diffcov_main.columns = ['type', 'reads']
    reads_diffcov_main['pos'] = np.arange(1, len(reads_diffcov_main) + 1).astype(str)
    reads_diffcov_main['id'] = reads_diffcov_main[["pos","type"]].agg("-".join, axis=1)
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

def get_assembler(assembler):
    if assembler == "metaflye": return "tmp/assembly/assembly.fasta"
    if assembler == "metamdbg": return "tmp/assembly/contigs.fasta"
    if assembler == "myloasm": return "tmp/assembly/assembly.fa"
    else: return "tmp/assembly/assembly_custom.fa"

def get_assembly(mode,assembler,medaka_status,curation):
    if curation == "TRUE": return "tmp/curation/asm_curated.fasta"
    elif mode == "PacBio-HiFi": return get_assembler(assembler)
    elif medaka_status == "TRUE": return "tmp/polishing/asm_pol.fasta"
    else: return get_assembler(assembler)

def get_asm_info(assembler,assembler_config,curation):
    if curation == "TRUE": return "tmp/curation/assembly_info.tsv"
    elif assembler == assembler_config: return "tmp/assembly/assembly_info.tsv"
    else: return assembler

def get_myloasm_params(myloasm_extra):
    if myloasm_extra == "FALSE": return " "
    else: return re.sub(","," ",myloasm_extra)

def get_prok(tiara_status):
    if tiara_status == "TRUE": return "contigs_filt_prok.txt"
    else: return "whokaryote/prokaryote_contig_headers.txt"

def get_euk(tiara_status):
    if tiara_status == "TRUE": return "contigs_filt_euk.txt"
    else: return "whokaryote/eukaryote_contig_headers.txt"

def get_map_mode(reads,np_map_mode,pb_map_mode,il_map_mode):
    if "NP" in reads: return np_map_mode
    if "PB" in reads: return pb_map_mode
    if "IL" in reads: return il_map_mode

def get_map_ident(reads,np_map_ident,pb_map_ident,il_map_ident):
    if "NP" in reads: return np_map_ident
    if "PB" in reads: return pb_map_ident
    if "IL" in reads: return il_map_ident

def get_map_cut1(reads):
    if reads.startswith("1-"): return "1,2,3,4,5"
    else: return "4,5"

def get_map_cut2(reads):
    if reads.startswith("1-"): return "1,4"
    else: return "4"

def get_semibin_cov(loc,sample,reads_diffcov,round,semibin_mod):
    if reads_diffcov == "none": return ''.join([str(loc), "/", str(sample), "/tmp/binning/mapping/*.bam", " --environment ", str(semibin_mod), " "])
    if round == "3": return ''.join([str(loc), "/", str(sample), "/tmp/binning/mapping/1-*.bam", " --environment ", str(semibin_mod), " "])
    else: return ''.join([str(loc), "/", str(sample), "/tmp/binning/mapping/*.bam "])

def get_rounds(binmode):
    if binmode == "extended": return list(range(1, 5))
    else: return 1

def get_qc_in(loc,sample,round, binmode):
    if binmode == "fast": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_1/semibin/bins_info.tsv"])
    if round == "2": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/comebin/comebin_res/comebin_res.tsv"])
    elif round == "4": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/metabat2/bins_metabat2.MemberMatrix.txt"])
    else: return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/binette/final_bins_quality_reports.tsv"])

def get_qc_dir(loc,sample,round,binmode):
    if binmode == "fast": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_1/semibin/output_bins"])
    if round == "2": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/comebin/comebin_res/comebin_res_bins"])
    elif round == "4": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/metabat2"])
    else: return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/binette/final_bins"])

def get_compl(round,binmode,compl1,compl2):
    if round == "1" and binmode == "extended": return compl1
    else: return compl2

def get_cont(round,binmode,cont1,cont2):
    if round == "1" and binmode == "extended": return cont1
    else: return cont2

onstart:
    from snakemake.utils import min_version
    min_version("9.9.0")
    if not os.path.exists(fastq): sys.exit(print("Read input (((",fastq,"))) not found. Aborting..."))
    if not os.path.isdir(loc): sys.exit(print("Provided path for output (((",loc,"))) not found. Aborting..."))
    if reads_diffcov != "none" and not os.path.exists(reads_diffcov): sys.exit(print("Dataframe for differential coverage binning (((",reads_diffcov,"))) not found. Aborting..."))
    if not mode in modes: sys.exit(print("Provided workflow mode (((",mode,"))) not recognised. Aborting..."))
    if not config["assembler"] in assemblers: sys.exit(print("Provided assembler (((",config["assembler"],"))) not recognised. Aborting..."))
    if not binmode in binmodes: sys.exit(print("Provided binning mode (((",binmode,"))) not recognised. Aborting..."))
    if not semibin_mod in semibin_mods: sys.exit(print("Provided model for SemiBin (((",semibin_mod,"))) not recognised. Aborting..."))
    if not semibin_prot in semibin_prots: sys.exit(print("Provided gene predictor for SemiBin (((",semibin_prot,"))) not recognised. Aborting..."))
    if len(os.path.join(loc, sample)) > 85: sys.exit(print("Path for provided output too long: (((",os.path.join(loc, sample),")))\nPlease re-run with different output location or shorter name."))
    if not os.path.exists(os.path.join(loc, sample)): os.makedirs(os.path.join(loc, sample))
    if not os.path.exists(os.path.join(loc, sample, "results")): os.makedirs(os.path.join(loc, sample, "results"))
    if not os.path.exists(os.path.join(loc, sample, "tmp")): os.makedirs(os.path.join(loc, sample, "tmp"))
    if not os.path.exists(os.path.join(loc, sample, "tmp", "logs")): os.makedirs(os.path.join(loc, sample, "tmp", "logs"))
    if not os.path.exists(os.path.join(loc, sample, "tmp", "dep_mmlong2-lite.csv")): 
	    with open(os.path.join(loc, sample, "tmp", "dep_mmlong2-lite.csv"), 'w') as f:
		    f.write("dependency,version\n")

onsuccess:
    from datetime import datetime
    now = datetime.now()
    print("MAG recovery with pipeline version",wf_v,"completed at",now.strftime("%Y/%m/%d %H:%M:%S"))
    print("Thank you for using mmlong2-lite")

onerror:
    print("An error has occurred. Inspect Snakemake log files for troubleshooting.")

rule Finalise:
    conda: config["env_1"]
    input:
        quast=os.path.join(loc, sample, "tmp/binning/quast.tsv"),
        abund=os.path.join(loc, sample, "tmp/binning/bin_abund.tsv"),
        cov=os.path.join(loc, sample, "tmp/binning/bin_cov.tsv"),
        checkm1=os.path.join(loc, sample, "tmp/binning/checkm1.tsv"),
        checkm2=os.path.join(loc, sample, "tmp/binning/checkm2.tsv"),
        clean=os.path.join(loc, sample, "tmp/logs/cleanup.txt"),
    output:
        df1=os.path.join(loc, sample, "tmp/binning/bins_mmlong2-lite.tsv"),
        df2=os.path.join(loc, sample, "results", f"{sample}_bins.tsv")
    shell:
        """
        if ! grep -q "r-base" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "r-base " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if [ -f "{loc}/{sample}/results/dependencies.csv" ]; then rm {loc}/{sample}/results/dependencies.csv; fi
        rsync {loc}/{sample}/tmp/dep_mmlong2-lite.csv {loc}/{sample}/results/dependencies.csv
		
        R --no-echo --silent --args << 'make_df'
        quast <- read.delim("{input.quast}", sep="\t", header=T)
        colnames(quast) <- c("bin","contigs","longest_contig","genome_size","gc","contig_n50","contig_n90","aun","n_per_100kb")
        abund <- read.delim("{input.abund}", sep="\t", header=T)
        colnames(abund) <- c("bin","r_abund")
        cov <- read.delim("{input.cov}", sep="\t", header=T)
        colnames(cov) <- c("bin","cov")
        checkm1 <- read.delim("{input.checkm1}", sep="\t", header=T)
        checkm1 <- checkm1[, c("Bin.Id","Completeness","Contamination","Strain.heterogeneity")]
        colnames(checkm1) <- c("bin","completeness_checkm1","contamination_checkm1","strain_heterogeneity_checkm1")
        checkm2 <- read.delim("{input.checkm2}", sep="\t", header=T)
        colnames(checkm2) <- c("bin","completeness_checkm2","contamination_checkm2")
        bins <- merge(checkm1,merge(checkm2,merge(quast,merge(cov,abund,by="bin"),by="bin"), by="bin"),by="bin")
        bins$wf_name <- "{sample}"
        bins$wf_read_mode <- "{mode}"
        bins$wf_binning_mode <- "{binmode}"
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
        asm=os.path.join(loc, sample, "tmp/assembly/assembly.fasta"),
        info=os.path.join(loc, sample, get_asm_info("metaflye",config["assembler"],"FALSE")),
    params:
        cov=config["flye_cov"],
        ovlp=config["flye_ovlp"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_assembly_metaflye.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        if [ {mode} == "Nanopore-simplex" ]; then flye_opt="--nano-hq"; fi
        if [ {mode} == "PacBio-HiFi" ]; then flye_opt="--read-error 0.01 --pacbio-hifi"; fi
        if [ {params.ovlp} -eq 0 ]; then flye_ovlp=""; else flye_ovlp="--min-overlap {params.ovlp}"; fi
        flye $flye_opt {input} --out-dir {loc}/{sample}/tmp/assembly --threads {threads} --meta $flye_ovlp --extra-params min_read_cov_cutoff={params.cov}
        sed 1d {loc}/{sample}/tmp/assembly/assembly_info.txt | cut -f1,2,3,4 > {output.info}
        if ! grep -q "flye" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "flye " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """

rule Assembly_metaMDBG:
    conda: config["env_4"]
    input:
        fastq
    output:
        asm=os.path.join(loc, sample, "tmp/assembly/contigs.fasta"),
        info=os.path.join(loc, sample, get_asm_info("metamdbg",config["assembler"],"FALSE")),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_assembly_metamdbg.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        if [ {mode} == "Nanopore-simplex" ]; then mdbg_opt="--in-ont"; else mdbg_opt="--in-hifi";fi
        metaMDBG asm $mdbg_opt {input} --out-dir {loc}/{sample}/tmp/assembly --threads {threads}
        unpigz {output.asm}.gz
        grep ">" {output.asm} | cut -c2- > {loc}/{sample}/tmp/assembly/contigs.txt
        sed -e "s/length=//" -e "s/coverage=//" -e "s/circular=//" -e "s/yes/Y/" -e "s/no/N/" -e 's/ /\t/g' {loc}/{sample}/tmp/assembly/contigs.txt > {output.info}
        if ! grep -q "metamdbg" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "metamdbg " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """

rule Assembly_myloasm:
    conda: config["env_4"]
    input:
        fastq
    output:
        asm=os.path.join(loc, sample, "tmp/assembly/assembly.fa"),
        info=os.path.join(loc, sample, get_asm_info("myloasm",config["assembler"],"FALSE")),
    params:
        cov=config["myloasm_cov"],
        ovlp=config["myloasm_ovlp"],
        bloom=config["myloasm_bloom"],
        circ_prob=config["circ_prob"],
        extra=get_myloasm_params(config["myloasm_extra"]),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_assembly_myloasm.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        if [ {mode} == "Nanopore-simplex" ]; then mylo_opt=""; else mylo_opt="--hifi";fi
        if [ {params.ovlp} -eq 500 ]; then mylo_ovlp=""; else mylo_ovlp="--min-ol {params.ovlp}"; fi
        if [ {params.cov} -eq 1 ]; then mylo_cov=""; else mylo_cov="--absolute-coverage-threshold {params.cov}"; fi
        myloasm {input} -o {loc}/{sample}/tmp/assembly -t {threads} -b {params.bloom} --clean-dir $mylo_cov $mylo_ovlp $mylo_opt {params.extra}
        sed '/^>/ s/_/ /g' {loc}/{sample}/tmp/assembly/assembly_primary.fa > {output.asm}
        grep ">" {loc}/{sample}/tmp/assembly/assembly_primary.fa | cut -c2- > {loc}/{sample}/tmp/assembly/contigs.txt
        sed -e 's/_/\t/g' -e "s/len-//" -e "s/circular-//" -e "s/yes/Y/" -e "s/no/N/" -e "s/possibly/{params.circ_prob}/" -e "s/depth-//" {loc}/{sample}/tmp/assembly/contigs.txt | sed 's/-.*//' | awk -F'\t' '{{t=$3; $3=$4; $4=t; print}}' OFS='\t' > {output.info}
        if ! grep -q "myloasm" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "myloasm " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """

rule Assembly_custom:
    conda: config["env_2"]
    input:
        config["custom_assembly"],
    output:
        asm=os.path.join(loc, sample, "tmp/assembly/assembly_custom.fa"),
        info=os.path.join(loc, sample, get_asm_info("custom",config["assembler"],"FALSE")),
    params:
        info=config["assembly_info"],
        circ_prob=config["circ_prob"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_assembly_custom.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        mkdir {loc}/{sample}/tmp/assembly

        if head -c2 {input} | grep -q $'\\x1f\\x8b'
        then unpigz -c {input} > {output.asm}
        else rsync {input} {output.asm}; fi
        grep ">" {output.asm} | cut -c2- > {loc}/{sample}/tmp/assembly/contigs.txt

        if grep -q "circular=" {output.asm}
        then sed  -e "s/length=//" -e "s/coverage=//" -e "s/circular=//" -e "s/yes/Y/" -e "s/no/N/" -e 's/ /\t/g' {loc}/{sample}/tmp/assembly/contigs.txt > {output.info}; fi

        if grep -q "_circular-" {output.asm}; then
        sed -e 's/_/\t/g' -e "s/len-//" -e "s/circular-//" -e "s/yes/Y/" -e "s/no/N/" -e "s/possibly/{params.circ_prob}/" -e "s/depth-//" {loc}/{sample}/tmp/assembly/contigs.txt | sed 's/-.*//' | awk -F'\t' '{{t=$3; $3=$4; $4=t; print}}' OFS='\t' > {output.info}
        sed -i '/^>/ s/_/ /g' {output.asm}; fi

        if [ ! "{params.info}" == "FALSE" ]; then 
        if grep -q '^#' {params.info}; then sed 1d {params.info} | cut -f1-4 > {output.info}; else rsync {params.info} {output.info}; fi; fi

        if [ ! -f "{output.info}" ]; then sed -i 's/^\\(>[^ ]*\\).*/\\1/' {output.asm} && seqkit fx2tab -l -n {output.asm} | awk '{{print $1"\t"$2"\t0\tN"}}' > {output.info}; fi
        """

rule Polishing_prep:
    conda: config["env_9"]
    input:
        os.path.join(loc, sample, get_assembler(config["assembler"]))
    output:
        expand("{loc}/{sample}/tmp/polishing/lin_contigs_id_{split}.txt",sample=sample,loc=loc,split=[f"{i:02}" for i in range(1, get_splits(proc,config["medaka_split"]) + 1)]),
        os.path.join(loc, sample, "tmp/polishing/calls_to_draft.bam")
    params:
        splits=get_splits(proc,config["medaka_split"]),
        ram_reads=config["minimap_ram"],
        asm=get_assembler(config["assembler"]),
        seed=config["seed"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_polishing_prep.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/polishing" ]; then rm -r {loc}/{sample}/tmp/polishing; fi
        if [ -f "{loc}/{sample}/{params.asm}.map-ont.mmi" ]; then rm {loc}/{sample}/{params.asm}.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/{params.asm}.fai" ]; then rm {loc}/{sample}/{params.asm}.fai; fi
		mkdir {loc}/{sample}/tmp/polishing
        
        grep ">" {input} | cut -c 2- | awk '{{print $1}}' | shuf --random-source=<(yes {params.seed}) > {loc}/{sample}/tmp/polishing/ids.txt
        split -n l/{params.splits} --numeric-suffixes=1 --additional-suffix=.txt -d {loc}/{sample}/tmp/polishing/ids.txt {loc}/{sample}/tmp/polishing/contigs_id_
        for file in {loc}/{sample}/tmp/polishing/contigs_id_*.txt; do filename=lin_$(basename $file)
        awk 'BEGIN {{ ORS = " " }} {{ print }}' $file > {loc}/{sample}/tmp/polishing/$filename; done
        mini_align -I {params.ram_reads}G -i {fastq} -r {input} -m -p {loc}/{sample}/tmp/polishing/calls_to_draft -t {threads}
        """

rule Polishing_consensus:
    conda: config["env_9"]
    input:
        "{contigs}.txt"
    output:
        "{contigs}.hdf"
    params:
        medaka_model=config["medak_mod_pol"],
        medaka_batch=config["medaka_batch"],
    benchmark:
        "{contigs}.tsv"
    resources: usage=config["medaka_usage"]
    threads: proc_sub
    shell:
        """
        if [ -f "{output}" ]; then rm {output}; fi
        medaka inference {loc}/{sample}/tmp/polishing/calls_to_draft.bam {output} --model {params.medaka_model} --batch {params.medaka_batch} --threads {threads} --region $(cat {input})
        """	

rule Polishing_stitch:
    conda: config["env_9"]
    input:
        hdf=expand("{contigs}.hdf", contigs=get_contigs()),
        asm=os.path.join(loc, sample, get_assembler(config["assembler"])),
    output:
        os.path.join(loc, sample, "tmp/polishing/asm_pol.fasta")
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_polishing_stitch.tsv")
    threads: proc
    shell:
        """
        medaka sequence {input.hdf} {input.asm} {output} --threads {threads}
        if [ -f "{loc}/{sample}/tmp/polishing/lin_contigs_id_01.tsv" ]; then
        mv {loc}/{sample}/tmp/polishing/*.tsv {loc}/{sample}/tmp/logs
        for f in {loc}/{sample}/tmp/logs/lin_contigs_id_*.tsv; do mv "$f" "$(dirname "$f")/usage_polishing_consensus-${{f##*/lin_contigs_id_}}"; done; fi
        if ! grep -q "medaka" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then pip list | grep -w "medaka" | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi        
        """			

rule Curation_map:
    conda: config["env_2"]
    input:
        reads=fastq,
        asm=os.path.join(loc, sample, get_assembly(mode,config["assembler"],config["medaka_status"],"FALSE")),
    output:
        os.path.join(loc, sample, "tmp/curation/mapping.bam")
    params:
        map=lambda wildcards: get_map_mode(get_read_type(mode),config["minimap_np"],config["minimap_pb"],config["minimap_il"]),
        ram_ref=config["minimap_ref"],
        ram_reads=config["minimap_ram"],
        ratio=config["minimap_ratio"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_curation_map.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/curation" ]; then rm -r {loc}/{sample}/tmp/curation; fi
        mkdir {loc}/{sample}/tmp/curation
        minimap2 -I {params.ram_reads}G -K {params.ram_ref}G -t {threads} --secondary-seq -p {params.ratio} -ax {params.map} {input.asm} {input.reads} | samtools view --threads $(({threads} / 2)) -Sb -F 4 - | samtools sort --threads $(({threads} / 2)) --write-index -o {output}##idx##{output}.bai -
        if ! grep -q "minimap2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "minimap2 " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if ! grep -q "samtools" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "samtools " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi        
        """	

rule Curation_screening:
    conda: config["env_2"]
    input:
        os.path.join(loc, sample, "tmp/curation/mapping.bam"),
    output:
        cov=os.path.join(loc, sample, "tmp/curation/errors-zero_cov.txt"),
        cov_filt=os.path.join(loc, sample, "tmp/curation/errors-zero_cov.tsv"),
        clip=os.path.join(loc, sample, "tmp/curation/errors-clipping.txt"),
        clip_filt=os.path.join(loc, sample, "tmp/curation/errors-clipping.tsv"),
    params:
        clip_ratio=config["clip_ratio"],
        min_dist=config["min_dist"],
        min_len_zerocov=config["min_len_zerocov"],
        min_cov_clip=config["min_cov_clip"],
        min_cov_clip_all=config["min_cov_clip_all"],
        apptainer=config["apptainer_status"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_curation_screening.tsv")
    shell:
        """
        if [ "{params.apptainer}" == "FALSE" ]; then if ! command -v bam_error_detector &> /dev/null; then cargo install --git https://github.com/bluenote-1577/rust-anvio-mis --root $CONDA_DEFAULT_ENV; fi; fi

        org_dir=$(pwd) && cd {loc}/{sample}/tmp/curation
        bam_error_detector --clipping-ratio {params.clip_ratio} --min-dist-to-end {params.min_dist} {input} errors
        cd $org_dir

        awk -F'\t' '$4 > {params.min_len_zerocov}' {output.cov} > {output.cov_filt}
        awk -F'\t' '$5 >= {params.min_cov_clip_all} && $6 >= {params.min_cov_clip}' {output.clip} > {output.clip_filt}
        if ! grep -q "rust" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "rust " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi      
        if ! grep -q "bam_error_detector" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then bam_error_detector --version | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi    
        """	

rule Curation_aggregate:
    conda: config["env_1"]
    input:
        asm=os.path.join(loc, sample, "tmp/assembly/assembly_info.tsv"),
        cov=os.path.join(loc, sample, "tmp/curation/errors-zero_cov.tsv"),
        clip=os.path.join(loc, sample, "tmp/curation/errors-clipping.tsv"),
    output:
        df=os.path.join(loc, sample, "tmp/curation/contigs.tsv"),
        con=os.path.join(loc, sample, "tmp/curation/contigs.txt"),
        asm=os.path.join(loc, sample, "tmp/curation/assembly_info.tsv"),
    params:
        len=config["min_contig_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_curation_aggregate.tsv")
    shell:
        """
        if ! grep -q "r-base" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "r-base " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if ! grep -q "r-tidyverse" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "r-tidyverse " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  

        R --no-echo --silent --args << 'curation'
        library(tidyverse)
        cov <- read.delim("{input.cov}", sep="\t", header=T)
        if (nrow(cov) == 0) {{cov <- tibble(contig="contig_placeholder", length=1000, range="900-1000", range_size=100)}}

        cov <- cov[order(cov$range,decreasing=FALSE),]
        cov <- cov[order(cov$contig,decreasing=FALSE),]
        cov_ <- cov %>% separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>% select(contig, length, start, end)

        cov_ <- cov_ %>% arrange(contig, start, end) %>% group_by(contig, length) %>% mutate(group = cumsum(coalesce(start > lag(end) + 1, TRUE))) %>%
                group_by(contig, length, group) %>% summarise(start = min(start), end = max(end), .groups = "drop")

        cov_ <- cov_ %>% group_by(contig, length) %>% arrange(start) %>% summarise(safe = list({{
                        starts <- c(0, end)
                        ends <- c(start, unique(length))
                        tibble(start = starts, end = ends) %>% filter(start < end) }}), .groups = "drop") %>% unnest(safe) %>% select(contig, start, end)

        clip <- read.delim("{input.clip}", sep="\t", header=T)
        if (nrow(clip) == 0) {{clip <- tibble(contig="contig_placeholder", length=1000, pos=500, relative_pos=0.5, cov=1, clipping=1, clipping_ratio=1)}}

        clip <- clip[order(clip$pos,decreasing=FALSE),]
        clip <- clip[order(clip$contig,decreasing=FALSE),]
        clip_ <- clip %>% group_by(contig) %>% summarise(length = max(length), breakpoints = list(sort(unique(pos))))

        clip_ <- clip_ %>% rowwise() %>% mutate(chunks = list({{
                        bp <- breakpoints
                        starts <- c(0, bp)
                        ends <- c(bp - 1, length)
                        data.frame(start = starts, end = ends) }})) %>% tidyr::unnest(chunks) %>% select(contig, start, end)

        sub <- bind_rows(cov_, clip_) %>% pivot_longer(cols = c(start, end), names_to = "type", values_to = "bp") %>% distinct() %>% arrange(contig, bp) %>%
               group_by(contig) %>% summarise(bps = list(bp), .groups = "drop") %>% mutate(intervals = map(bps, ~{{ tibble(
                        start = head(.x, -1), end = tail(.x, -1)) }})) %>% select(contig, intervals) %>% unnest(intervals)

        regions <- bind_rows(cov_, clip_)
        sub <- sub %>% semi_join(regions, by = "contig") %>% rowwise() %>% filter(any(
                        regions$contig == contig &
                        regions$start <= end &
                        regions$end >= start )) %>% ungroup() %>% filter(abs(end - start) >= {params.len})
        write.table(sub,"{output.df}",quote=F,row.names=FALSE,col.names=FALSE,sep="\t")

        asm <- read.delim("{input.asm}", sep="\t", header=F)
        asm_ <- asm[! asm$V1 %in% regions$contig,]

        sub_ <- sub
        sub_$len <- sub_$end - sub_$start 
        sub_$start <- sub_$start + 1
        sub_$cov <- NA
        sub_$circ <- "N"
        sub_$contig_ <- paste0(sub_$contig,"_",sub_$start,"-",sub_$end)

        colnames(asm_) <- c("contig_","len","cov","circ")
        asm_ <- rbind(asm_,sub_[,c("contig_","len","cov","circ")])
        write.table(asm_,"{output.asm}",quote=F,row.names=FALSE,col.names=FALSE,sep="\t")
        write.table(asm[! asm$V1 %in% regions$contig,]$V1,"{output.con}",quote=F,row.names=FALSE,col.names=FALSE,sep="\t")  
        """	

rule Curation_selection:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, get_assembly(mode,config["assembler"],config["medaka_status"],"FALSE")),
        df=os.path.join(loc, sample, "tmp/curation/contigs.tsv"),
        con=os.path.join(loc, sample, "tmp/curation/contigs.txt"),
    output:
        os.path.join(loc, sample, "tmp/curation/asm_curated.fasta")
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_curation_selection.tsv")
    shell:
        """
        seqkit grep -f {input.con} {input.asm} > {output}
        seqkit subseq --bed {input.df} {input.asm} | sed '/^>/ s/...$//' >> {output}
        """	

rule Filtering_length:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, get_assembly(mode,config["assembler"],config["medaka_status"],config["curation_status"])),
        info=os.path.join(loc, sample, get_asm_info("TRUE","TRUE",config["curation_status"])),
    output:
        asm=os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta"),
        info=os.path.join(loc, sample, "tmp/filtering/assembly_info.tsv"),
    params:
        len=config["min_contig_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_filtering_length.tsv")
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/filtering" ]; then rm -r {loc}/{sample}/tmp/filtering; fi
        mkdir {loc}/{sample}/tmp/filtering
        seqkit seq -m {params.len} {input.asm} | seqkit replace -p "\\s.+" | seqkit replace -p ^ -r {sample}_ > {output.asm}
        rsync {output.asm} {loc}/{sample}/results/{sample}_assembly.fasta
        sed 's/^/{sample}_/' {input.info} > {output.info}
        if ! grep -q "seqkit" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "seqkit " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi        
        """	

rule Filtering_tiara:
    conda: config["env_3"]
    input:
        os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta")
    output:
        prok=os.path.join(loc, sample, "tmp/filtering/contigs_filt_prok.txt"),
        euk=os.path.join(loc, sample, "tmp/filtering/contigs_filt_euk.txt"),
    params:
        len=config["min_contig_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_filtering_tiara.tsv")
    threads: proc
    shell:
        """
        tiara -i {input} -t {threads} -m {params.len} -o {loc}/{sample}/tmp/filtering/tiara.tsv
        cut -f1,2 {loc}/{sample}/tmp/filtering/tiara.tsv | (grep -e "prokarya" -e "bacteria" -e "archaea" -e "unknown" - || true) | cut -f1 | sort > {output.prok}
        cut -f1,2 {loc}/{sample}/tmp/filtering/tiara.tsv | (grep -e "eukarya" -e "organelle" -e "unknown" - || true) | cut -f1 | sort > {output.euk}
        if ! grep -q "tiara" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "tiara " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi     
        """

rule Filtering_whokaryote:
    conda: config["env_3"]
    input:
        os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta")
    output:
        prok=os.path.join(loc, sample, "tmp/filtering/whokaryote/prokaryote_contig_headers.txt"),
        euk=os.path.join(loc, sample, "tmp/filtering/whokaryote/eukaryote_contig_headers.txt"),
    params:
        len=config["min_contig_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_filtering_whokaryote.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/filtering/whokaryote" ]; then rm -r {loc}/{sample}/tmp/filtering/whokaryote; fi
        whokaryote.py --contigs {input} --threads {threads} --minsize {params.len} --outdir {loc}/{sample}/tmp/filtering/whokaryote 
        if ! grep -q "whokaryote" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "whokaryote " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi     
        """

rule Filtering_eukaryotes:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta"),
        prok=os.path.join(loc, sample, "tmp/filtering", get_prok(config["tiara_status"])),
        euk=os.path.join(loc, sample, "tmp/filtering", get_euk(config["tiara_status"])),
    output:
        prok=os.path.join(loc, sample, "tmp/filtering/asm_filt_prok.fasta"),
        euk=os.path.join(loc, sample, "tmp/filtering/asm_filt_euk.fasta"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_filtering_eukaryotes.tsv")
    shell:
        """
        seqkit grep -f {input.prok} {input.asm} > {output.prok}
        seqkit grep -f {input.euk} {input.asm} > {output.euk} 
        if ! grep -q "seqkit" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "seqkit " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi    
        """

rule Singletons_circ:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta"),
        info=os.path.join(loc, sample, "tmp/filtering/assembly_info.tsv"),
        prok=os.path.join(loc, sample, "tmp/filtering/asm_filt_prok.fasta"),
    output:
        os.path.join(loc, sample, "tmp/binning/singl/contig_c.tsv")
    params:
        mag_len=config["min_mag_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_singletons_circ.tsv")
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/binning" ]; then mkdir {loc}/{sample}/tmp/binning; fi
        if [ ! -d "{loc}/{sample}/tmp/binning/singl" ]; then mkdir {loc}/{sample}/tmp/binning/singl; fi
        if [ ! -d "{loc}/{sample}/tmp/binning/singl/innit" ]; then mkdir {loc}/{sample}/tmp/binning/singl/innit; fi
        
        awk -F "\t" '{{ if (($2 >= {params.mag_len}) && ($4 == "Y")) {{print $1 "\t" "{sample}.bin.c." ++i; next}} }}' {input.info} > {output}
        find {loc}/{sample}/tmp/binning/singl/innit -name "{sample}.bin.c.*.fa" -type f -delete
        if [ $(cat {output} | wc -l) -ge 1 ]; then cat {output} | xargs -i --max-procs=1 -n 2 bash -c 'samtools faidx {input.asm} $0 > {loc}/{sample}/tmp/binning/singl/innit/$1.fa'; fi;
        
        find {loc}/{sample}/tmp/binning/singl/innit -name "*.fa" -type f -size -2k -delete
        if ! grep -q "samtools" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "samtools " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Singletons_lin:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta"),
        info=os.path.join(loc, sample, "tmp/filtering/assembly_info.tsv"),
        circ=os.path.join(loc, sample, "tmp/binning/singl/contig_c.tsv"),
    output:
        os.path.join(loc, sample, "tmp/binning/singl/contig_l.tsv")
    params:
        mag_len=config["min_smag_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_singletons_lin.tsv")
    shell:
        """
        awk -F "\t" '{{ if (($2 >= {params.mag_len}) && ($4 == "N")) {{print $1 "\t" "{sample}.bin.s." ++i; next}} }}' {input.info} > {output}
        find {loc}/{sample}/tmp/binning/singl/innit -name "{sample}.bin.s.*.fa" -type f -delete
        if [ $(cat {output} | wc -l) -ge 1 ]; then cat {output} | xargs -i --max-procs=1 -n 2 bash -c 'samtools faidx {input.asm} $0 > {loc}/{sample}/tmp/binning/singl/innit/$1.fa'; 
        elif [ $(cat {input.circ} | wc -l) -lt 1 ]; then seqkit fx2tab -nl {input.asm} | sort -k2,2nr | head -n1 | cut -f1 | seqkit grep -f - {input.asm} > {loc}/{sample}/tmp/binning/singl/innit/{sample}.bin.s.1.fa || true
        grep ">" {loc}/{sample}/tmp/binning/singl/innit/{sample}.bin.s.1.fa | cut -c 2- | awk '{{print $1 "\t" "{sample}.bin.s.1"}}' > {output}; fi
        find {loc}/{sample}/tmp/binning/singl/innit -name "*.fa" -type f -size -2k -delete
        """

rule Singletons_qc:
    conda: config["env_8"]
    input:
        os.path.join(loc, sample, "tmp/binning/singl/contig_l.tsv")
    output:
        os.path.join(loc, sample, "tmp/binning/singl/checkm2.tsv")
    params:
        compl_c=config["min_compl_circ"],
        compl_l=config["min_compl_lin"],
        cont=config["min_cont_singl"],
        apptainer=config["apptainer_status"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_singletons_qc.tsv")
    threads: proc
    shell:
        """
        if [ "{params.apptainer}" == "FALSE" ]; then if [ ! -f $CONDA_DEFAULT_ENV/CheckM2_database/uniref100.KO.1.dmnd ]; then checkm2 database --download --path $CONDA_DEFAULT_ENV; fi; fi 
        
        if [ -d "{loc}/{sample}/tmp/binning/singl/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/singl/checkm2; fi
        checkm2 predict -x .fa -i {loc}/{sample}/tmp/binning/singl/innit -o {loc}/{sample}/tmp/binning/singl/checkm2 -t {threads}
        cut -f1,2,3 {loc}/{sample}/tmp/binning/singl/checkm2/quality_report.tsv > {output}
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
        os.path.join(loc, sample, "tmp/binning/singl/checkm2.tsv")
    output:
        reads=os.path.join(loc, sample, "tmp/binning/mapping/reads.csv"),
        links=expand("{loc}/{sample}/tmp/binning/mapping/{reads}.lnk",sample=sample,loc=loc,reads=get_reads(reads_diffcov,fastq,mode)),			
    params:
        reads=fastq,
        type=get_read_type(mode),
        extra=reads_diffcov,
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_coverage_prep.tsv")
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
        ln -s $reads {loc}/{sample}/tmp/binning/mapping/${{count}}-${{type}}.lnk; done
        """

rule Coverage_map:
    conda: config["env_2"]
    input:
        "{loc}/{sample}/tmp/binning/mapping/{reads}.lnk"
    output:
        bam="{loc}/{sample}/tmp/binning/mapping/{reads}.bam",
        cov="{loc}/{sample}/tmp/binning/mapping/{reads}.tsv",
        tr1="{loc}/{sample}/tmp/binning/mapping/{reads}_tr1.tsv",
        tr2="{loc}/{sample}/tmp/binning/mapping/{reads}_tr2.tsv",
    params:
        map=lambda wildcards: get_map_mode(wildcards.reads,config["minimap_np"],config["minimap_pb"],config["minimap_il"]),
        ram_ref=config["minimap_ref"],
        ram_reads=config["minimap_ram"],
        ident=lambda wildcards: get_map_ident(wildcards.reads,config["np_map_ident"],config["pb_map_ident"],config["il_map_ident"]),
        tr1=lambda wildcards: get_map_cut1(wildcards.reads),
        tr2=lambda wildcards: get_map_cut2(wildcards.reads),
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_coverage_map-{reads}.tsv"
    resources: usage=config["mapping_usage"]
    threads: proc * (config["mapping_usage"]/100)
    shell:
        """
        if [ -f "{output.bam}" ]; then rm {output.bam}; fi
        minimap2 -I {params.ram_reads}G -K {params.ram_ref}G -t {threads} -ax {params.map} {loc}/{sample}/tmp/filtering/asm_filt_len.fasta {input} | samtools view --threads $(({threads} / 2)) -Sb -F 2308 - | samtools sort --threads $(({threads} / 2)) --write-index -o {output.bam}##idx##{output.bam}.bai -
        jgi_summarize_bam_contig_depths {output.bam} --percentIdentity {params.ident} --outputDepth {output.cov}
        cut -f{params.tr1} {output.cov} > {output.tr1}
        cut -f{params.tr2} {output.cov} > {output.tr2}
        """

rule Coverage_aggregate:
    conda: config["env_2"]
    input:
        cov1=expand("{loc}/{sample}/tmp/binning/mapping/{reads}_tr1.tsv",sample=sample,loc=loc,reads=get_reads(reads_diffcov,fastq,mode)),
        cov2=expand("{loc}/{sample}/tmp/binning/mapping/{reads}_tr2.tsv",sample=sample,loc=loc,reads=get_reads(reads_diffcov,fastq,mode)),
    output:
        cov1=os.path.join(loc, sample, "tmp/binning/mapping/cov_all.tsv"),
        cov2=os.path.join(loc, sample, "tmp/binning/mapping/cov_all_sub.tsv"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_coverage_aggregate.tsv")
    shell:
        """
        paste -d "\t" $(ls -v {input.cov1}) > {output.cov1}
        paste -d "\t" $(ls -v {input.cov2}) | sed "s/contigName/contigname/" > {output.cov2}
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
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_prep-comebin-innit-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/mapping" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/mapping; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/mapping
        awk 'FNR==NR{{contigs[$1]; next}} /^@SQ/ {{contig=$2; sub(/^SN:/, "", contig); if (contig in contigs) print}}' "{input}" <(samtools view --threads {threads} -H "{loc}/{sample}/tmp/binning/mapping/1-{params.type}.bam") > {output}
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
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_prep-comebin-{reads}-R{round}.tsv"
    resources: usage=config["mapping_usage"]
    threads: proc * (config["mapping_usage"]/100)
    shell:
        """
        grep -wFf "{input.contigs}" <(samtools view --threads {threads} "{input.bam}") | cat {input.header} - | samtools view --threads {threads} -b -o {output} -
        """

rule Binning_prep_innit:
    conda: config["env_2"]
    input:
        cov1=os.path.join(loc, sample, "tmp/binning/mapping/cov_all.tsv"),
        cov2=os.path.join(loc, sample, "tmp/binning/mapping/cov_all_sub.tsv"),
    output:
        id1=os.path.join(loc, sample, "tmp/binning/round_1/contigs_w_singl.txt"),
        id2=os.path.join(loc, sample, "tmp/binning/round_1/contigs_wo_singl.txt"),
        id3=os.path.join(loc, sample, "tmp/binning/round_1/contigs.txt"),
        contigs=os.path.join(loc, sample, "tmp/binning/round_1/contigs.fasta"),
        cov1=os.path.join(loc, sample, "tmp/binning/round_1/cov.tsv"),
        cov2=os.path.join(loc, sample, "tmp/binning/round_1/cov_sub.tsv"),
    params:
        cov=config["min_contig_cov"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_prep-R1.tsv")
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_1" ]; then rm -r {loc}/{sample}/tmp/binning/round_1; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_2" ]; then rm -r {loc}/{sample}/tmp/binning/round_2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_3" ]; then rm -r {loc}/{sample}/tmp/binning/round_3; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_4" ]; then rm -r {loc}/{sample}/tmp/binning/round_4; fi
        mkdir {loc}/{sample}/tmp/binning/round_1
        
        sed 1d {loc}/{sample}/tmp/binning/mapping/cov_all_sub.tsv | awk '$2 >= {params.cov}' | cut -f1 | sort > {output.id1}
        comm -23 <(sort {output.id1}) <(sort {loc}/{sample}/tmp/binning/singl/binned.txt) > {output.id2}
        seqkit grep -f {output.id2} {loc}/{sample}/tmp/filtering/asm_filt_prok.fasta | seqkit sort --by-name - > {output.contigs}
        grep ">" {output.contigs} | cut -c2- > {output.id3}
        
        head -n1 {input.cov1} > {output.cov1}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id3} {input.cov1} | sort >> {output.cov1}

        head -n1 {input.cov2} > {output.cov2}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id3} {input.cov2} | sort >> {output.cov2}
        """

rule Binning_prep_main:
    conda: config["env_2"]
    input:
        lambda wildcards: ''.join([str(wildcards.loc), "/", str(wildcards.sample), "/tmp/binning/round_",str(int(wildcards.round)-1),"/checkm2.tsv"])
    output:
        id="{loc}/{sample}/tmp/binning/round_{round}/contigs.txt",
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
        cov1="{loc}/{sample}/tmp/binning/round_{round}/cov.tsv",
        cov2="{loc}/{sample}/tmp/binning/round_{round}/cov_sub.tsv",
    params:
        round_curr="{round}",
        round_prev=lambda wildcards: str(int(wildcards.round)-1),
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_prep-R{round}.tsv"
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round_curr}" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round_curr}; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round_curr}
        
        comm -23 <(sort {loc}/{sample}/tmp/binning/round_{params.round_prev}/contigs.txt) <(sort {loc}/{sample}/tmp/binning/round_{params.round_prev}/binned.txt) > {output.id}
        seqkit grep -f {output.id} {loc}/{sample}/tmp/binning/round_{params.round_prev}/contigs.fasta | seqkit sort --by-name - > {output.contigs}

        head -n1 {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov.tsv > {output.cov1}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id} {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov.tsv | sort >> {output.cov1}

        head -n1 {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov_sub.tsv > {output.cov2}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id} {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov_sub.tsv | sort >> {output.cov2}
        """

rule Binning_metabat2:
    conda: config["env_2"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/metabat2/bins_metabat2.MemberMatrix.txt"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        mag_len=config["min_mag_len"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_metabat2-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/metabat2" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/metabat2; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/metabat2
        metabat2 -i {input} -a {loc}/{sample}/tmp/binning/round_{params.round}/cov.tsv -o {loc}/{sample}/tmp/binning/round_{params.round}/metabat2/bins_metabat2 -t {threads} -m {params.len} -s {params.mag_len} --saveCls --seed {params.seed}
        if ! grep -q "metabat2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "metabat2 " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi	
        """

rule Binning_vamb:
    conda: config["env_7"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/vamb/vae_clusters_unsplit.tsv"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        mag_len=config["min_mag_len"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_vamb-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/vamb" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/vamb; fi
        vamb bin default --fasta {input} --outdir {loc}/{sample}/tmp/binning/round_{params.round}/vamb -m {params.len} --minfasta {params.mag_len} -p {threads} --abundance_tsv {loc}/{sample}/tmp/binning/round_{params.round}/cov_sub.tsv -o "" --seed {params.seed}
        if ! grep -q "vamb" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "vamb " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi	
        """

rule Binning_semibin2:
    conda: config["env_5"]
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
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_semibin2-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/semibin" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/semibin; fi
        SemiBin2 single_easy_bin -i {input} -b {params.cov} -o {loc}/{sample}/tmp/binning/round_{params.round}/semibin -p {threads} -m {params.len} --minfasta-kbs $(({params.mag_len}/1000)) --self-supervised --sequencing-type long_read --engine cpu --orf-finder {semibin_prot} --compression none --random-seed {params.seed}
        if ! grep -q "semibin" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "semibin " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_binette:
    conda: config["env_8"]
    input:
        asm="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
        vamb="{loc}/{sample}/tmp/binning/round_{round}/vamb/vae_clusters_unsplit.tsv",
        mb2="{loc}/{sample}/tmp/binning/round_{round}/metabat2/bins_metabat2.MemberMatrix.txt",
        sb2="{loc}/{sample}/tmp/binning/round_{round}/semibin/bins_info.tsv"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/binette/final_bins_quality_reports.tsv",
    params:
        round="{round}",
        weight=config["binette_weight"],
        seed=config["seed"],
        compl=lambda wildcards: get_compl(wildcards.round, config["binmode"], config["min_compl_1"], config["min_compl_2"]),
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_binette-R{round}.tsv"
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/binette" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/binette; fi
        export PYTHONHASHSEED={params.seed}
        dir={loc}/{sample}/tmp/binning/round_{params.round}
        binette -d $dir/vamb/bins $dir/metabat2 $dir/semibin/output_bins -c {input.asm} -o $dir/binette -w {params.weight} -m {params.compl} -t {threads}
        if ! grep -q "binette" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "binette " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_comebin:
    conda: config["env_6"]
    input:
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
        bams=expand("{{loc}}/{{sample}}/tmp/binning/round_{{round}}/mapping/{reads}.bam", reads=get_reads(reads_diffcov,fastq,mode))
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/comebin/comebin_res/comebin_res.tsv"
    params:
        round="{round}",
        len=config["min_contig_len"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_comebin-R{round}.tsv"
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
    conda: config["env_8"]
    input:
        lambda wildcards: get_qc_in(wildcards.loc,wildcards.sample,wildcards.round,config["binmode"])
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/checkm2.tsv"
    params:
        round="{round}",
        dir=lambda wildcards: get_qc_dir(wildcards.loc,wildcards.sample,wildcards.round,config["binmode"]),
        compl=lambda wildcards: get_compl(wildcards.round, config["binmode"], config["min_compl_1"], config["min_compl_2"]),
        cont=lambda wildcards: get_cont(wildcards.round, config["binmode"], config["min_cont_1"], config["min_cont_2"]),
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_qc-R{round}.tsv"
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/bins" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/bins; fi

        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/bins
        
        if [ $(basename {input}) == "final_bins_quality_reports.tsv" ];
        then
            rsync {params.dir}/*.fa {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/.
            for bin in {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/*.fa; do mv "$bin" "$(dirname "$bin")/{sample}.bin.{params.round}.$(basename "$bin" | sed 's/^bin_//')"; done
            cut -f1,4,5 {input} | sed 's/^/{sample}.bin.{params.round}./' > {output};
        else
            n=1 && for i in {params.dir}/*.fa; do rsync $i {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/{sample}.bin.{params.round}.${{n}}.fa && n=$(($n+1)); done
            checkm2 predict -x .fa -i {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit -o {loc}/{sample}/tmp/binning/round_{params.round}/checkm2 -t {threads}
            cut -f1,2,3 {loc}/{sample}/tmp/binning/round_{params.round}/checkm2/quality_report.tsv > {output}
        fi
        
        awk -F "\t" '{{ if (($2 >= {params.compl}) && ($3 <= {params.cont})) {{print $1}} }}' {output} > {loc}/{sample}/tmp/binning/round_{params.round}/bins_keep.txt
        if [ $(cat {loc}/{sample}/tmp/binning/round_{params.round}/bins_keep.txt | wc -l) -ge 1 ]; then
        cat {loc}/{sample}/tmp/binning/round_{params.round}/bins_keep.txt | xargs -i --max-procs=1 bash -c 'rsync {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/{{}}.fa {loc}/{sample}/tmp/binning/round_{params.round}/bins/.'
        grep ">" {loc}/{sample}/tmp/binning/round_{params.round}/bins/*.fa | cut -f2 -d">" | sort > {loc}/{sample}/tmp/binning/round_{params.round}/binned.txt; fi
        """

rule Binning_aggregate:
    conda: config["env_2"]
    input:
        expand("{loc}/{sample}/tmp/binning/round_{round}/checkm2.tsv",sample=sample,loc=loc,round=get_rounds(config["binmode"])),
    output:
        os.path.join(loc, sample, "tmp/binning/checkm2.tsv"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_aggregate.tsv")
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/bins_innit; fi
        mkdir {loc}/{sample}/tmp/binning/bins_innit
        rsync {loc}/{sample}/tmp/binning/*/bins/*.fa {loc}/{sample}/tmp/binning/bins_innit/.
        
        head -n1 {loc}/{sample}/tmp/binning/singl/checkm2.tsv > {output}
        for stage in "singl" "round_1" "round_2" "round_3" "round_4" 
        do if [ -d "{loc}/{sample}/tmp/binning/$stage" ]; then if [ $(cat {loc}/{sample}/tmp/binning/$stage/bins_keep.txt | wc -l) -ge 1 ]; then grep -w -f {loc}/{sample}/tmp/binning/$stage/bins_keep.txt {loc}/{sample}/tmp/binning/$stage/checkm2.tsv >> {output}; fi; fi; done
        """

rule Binning_qc2:
    conda: config["env_6"]
    input:
         os.path.join(loc, sample, "tmp/binning/checkm2.tsv"),
    output:
        checkm1=os.path.join(loc, sample, "tmp/binning/checkm1.tsv"),
        bins=os.path.join(loc, sample, "tmp/binning/bins_keep.txt"),
    params:
        compl=config["min_compl_checkm1"],
        cont=config["min_cont_checkm1"],
        score=config["min_score_checkm1"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_qc2.tsv")
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
        os.path.join(loc, sample, "tmp/binning/checkm1.tsv"),
    output:
        cov=os.path.join(loc, sample, "tmp/binning/bin_cov.tsv"),
        abund=os.path.join(loc, sample, "tmp/binning/bin_abund.tsv"),
    params:
        cov_method=config["cov_method"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_summary_cov.tsv")
    threads: proc
    shell:
        """
        coverm genome -t {threads} -b {loc}/{sample}/tmp/binning/mapping/1-*.bam -d {loc}/{sample}/tmp/binning/bins -x fa -o {output.cov} -m {params.cov_method}
        coverm genome -t {threads} -b {loc}/{sample}/tmp/binning/mapping/1-*.bam -d {loc}/{sample}/tmp/binning/bins -x fa -o {output.abund} -m relative_abundance 
        if ! grep -q "coverm" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "coverm " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Summary_stats:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/binning/bin_cov.tsv"),
    output:
        links=os.path.join(loc, sample, "tmp/binning/contig_bin.tsv"),
        stats=os.path.join(loc, sample, "tmp/binning/quast.tsv"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_summary_stats.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/quast" ]; then rm -r {loc}/{sample}/tmp/binning/quast; fi
        grep -H ">" {loc}/{sample}/tmp/binning/bins/*.fa | xargs -n1 basename | sed "s/.fa:>/\t/" | awk -F'\t' '{{print $2 "\t" $1}}' > {output.links}
        quast.py {loc}/{sample}/tmp/binning/bins/*.fa -o {loc}/{sample}/tmp/binning/quast -t {threads} --no-plots --no-html --no-icarus
        cut -f1,14,15,16,17,18,19,20,23 {loc}/{sample}/tmp/binning/quast/transposed_report.tsv > {output.stats}
        if ! grep -q "quast" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "quast " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Summary_usage:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/binning/contig_bin.tsv"),
    output:
        os.path.join(loc, sample, "tmp/logs/summary_mmlong2-lite.tsv")
    threads: proc
    shell:
        """
        head -n1 {loc}/{sample}/tmp/logs/usage_summary_stats.tsv | sed 's/^/stage\tstep\t/' > {output}.tmp
        for log in {loc}/{sample}/tmp/logs/usage_*.tsv; do
            name=$(basename $log)
            IFS=_ read _ stage step <<< "${{name%.tsv}}"
            sed 1d $log | sed "s/^/${{stage}}\t${{step}}\t/" >> {output}.tmp
        done

        awk 'BEGIN {{
            order["stage"]=1; order["assembly"]=2; order["polishing"]=3; order["curation"]=4; order["filtering"]=5; order["singletons"]=6; order["coverage"]=7; order["binning"]=8; order["summary"]=9
        }}{{
            key = ($1 in order) ? order[$1] : 9999
            print key "\t" $0
        }}' {output}.tmp | sort -k1,1n | cut -f2- > {output}

        rm {output}.tmp
        rsync {output} {loc}/{sample}/results/{sample}_usage.tsv
        """

rule Cleanup:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/logs/summary_mmlong2-lite.tsv")
    output:
        os.path.join(loc, sample, "tmp/logs/cleanup.txt")
    params:
        cleanup=config["cleanup_status"],
    threads: proc
    shell:
        """
        if [ "{params.cleanup}" == "TRUE" ]; then

        size_pre=$(du -sb {loc}/{sample}/tmp | awk '{{printf "%.1f", $1/1024/1024/1024}}')
        files_pre=$(find {loc}/{sample}/tmp -type f -printf '.' | wc -c)
        echo "Status before cleanup (MAG production): $files_pre files and $size_pre GB of storage" > {output}

        if [ -d "{loc}/{sample}/tmp/assembly/00-assembly" ]; then rm -r {loc}/{sample}/tmp/assembly/00-assembly; fi
        if [ -d "{loc}/{sample}/tmp/assembly/10-consensus" ]; then rm -r {loc}/{sample}/tmp/assembly/10-consensus; fi
        if [ -d "{loc}/{sample}/tmp/assembly/20-repeat" ]; then rm -r {loc}/{sample}/tmp/assembly/20-repeat; fi
        if [ -d "{loc}/{sample}/tmp/assembly/30-contigger" ]; then rm -r {loc}/{sample}/tmp/assembly/30-contigger; fi
        if [ -d "{loc}/{sample}/tmp/assembly/40-polishing" ]; then rm -r {loc}/{sample}/tmp/assembly/40-polishing; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fasta.map-ont.mmi" ]; then rm {loc}/{sample}/tmp/assembly/assembly.fasta.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fasta" ]; then pigz {loc}/{sample}/tmp/assembly/assembly.fasta; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly_graph.gfa" ]; then pigz {loc}/{sample}/tmp/assembly/assembly_graph.gfa; fi

        if [ -d "{loc}/{sample}/tmp/assembly/0-cleaning_and_unitigs" ]; then rm -r {loc}/{sample}/tmp/assembly/0-cleaning_and_unitigs; fi
        if [ -d "{loc}/{sample}/tmp/assembly/1-light_resolve" ]; then rm -r {loc}/{sample}/tmp/assembly/1-light_resolve; fi
        if [ -d "{loc}/{sample}/tmp/assembly/2-heavy_path_resolve" ]; then rm -r {loc}/{sample}/tmp/assembly/2-heavy_path_resolve; fi
        if [ -d "{loc}/{sample}/tmp/assembly/3-mapping" ]; then rm -r {loc}/{sample}/tmp/assembly/3-mapping; fi
        if [ -d "{loc}/{sample}/tmp/assembly/alternate_assemblies" ]; then rm -r {loc}/{sample}/tmp/assembly/alternate_assemblies; fi
        if [ -d "{loc}/{sample}/tmp/assembly/assembly_graphs" ]; then rm -r {loc}/{sample}/tmp/assembly/assembly_graphs; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fa.map-ont.mmi" ]; then rm {loc}/{sample}/tmp/assembly/assembly.fa.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fa" ]; then pigz {loc}/{sample}/tmp/assembly/assembly.fa; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly_primary.fa" ]; then pigz {loc}/{sample}/tmp/assembly/assembly_primary.fa; fi

        if [ -d "{loc}/{sample}/tmp/assembly/tmp" ]; then rm -r {loc}/{sample}/tmp/assembly/tmp; fi
        if [ -f "{loc}/{sample}/tmp/assembly/contigs.fasta.map-ont.mmi" ]; then rm {loc}/{sample}/tmp/assembly/contigs.fasta.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/tmp/assembly/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/assembly/contigs.fasta; fi

        if [ -f "{loc}/{sample}/tmp/assembly/assembly_custom.fa.map-ont.mmi" ]; then rm {loc}/{sample}/tmp/assembly/assembly_custom.fa.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly_custom.fa" ]; then pigz {loc}/{sample}/tmp/assembly/assembly_custom.fa; fi

        if [ -f "{loc}/{sample}/tmp/polishing/ids.txt" ]; then rm {loc}/{sample}/tmp/polishing/*.txt; fi
        if [ -f "{loc}/{sample}/tmp/polishing/lin_contigs_id_01.hdf" ]; then rm {loc}/{sample}/tmp/polishing/lin_contigs_id_*.hdf; fi
        if [ -f "{loc}/{sample}/tmp/polishing/calls_to_draft.bam" ]; then rm {loc}/{sample}/tmp/polishing/calls_to_draft.bam*; fi
        if [ -f "{loc}/{sample}/tmp/polishing/asm_pol.fasta" ]; then pigz {loc}/{sample}/tmp/polishing/asm_pol.fasta; fi

        if [ -f "{loc}/{sample}/tmp/curation/asm_curated.fasta" ]; then pigz {loc}/{sample}/tmp/curation/asm_curated.fasta; fi
        if [ -f "{loc}/{sample}/tmp/filtering/asm_filt_len.fasta" ]; then pigz {loc}/{sample}/tmp/filtering/asm_filt_len.fasta; fi
        if [ -f "{loc}/{sample}/tmp/filtering/asm_filt_prok.fasta" ]; then pigz {loc}/{sample}/tmp/filtering/asm_filt_prok.fasta; fi
        if [ -f "{loc}/{sample}/tmp/filtering/asm_filt_euk.fasta" ]; then pigz {loc}/{sample}/tmp/filtering/asm_filt_euk.fasta; fi
        if [ -d "{loc}/{sample}/tmp/filtering/whokaryote" ]; then rm -r {loc}/{sample}/tmp/filtering/whokaryote; fi

        if [ -f "{loc}/{sample}/tmp/binning/round_2/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/binning/round_2/contigs.fasta; fi
        if [ -f "{loc}/{sample}/tmp/binning/round_3/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/binning/round_3/contigs.fasta; fi
        if [ -f "{loc}/{sample}/tmp/binning/round_4/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/binning/round_4/contigs.fasta; fi

        if [ -d "{loc}/{sample}/tmp/binning/singl/innit" ]; then rm -r {loc}/{sample}/tmp/binning/singl/innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_1/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_1/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_2/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_2/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_3/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_3/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_4/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_4/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/bins_innit; fi

        if [ -d "{loc}/{sample}/tmp/binning/singl/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/singl/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_1/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/round_1/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_2/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/round_2/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_4/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/round_4/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/checkm1" ]; then rm -r {loc}/{sample}/tmp/binning/checkm1; fi

        if [ -d "{loc}/{sample}/tmp/binning/round_1/binette/temporary_files" ]; then rm -r {loc}/{sample}/tmp/binning/round_1/binette/temporary_files; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_3/binette/temporary_files" ]; then rm -r {loc}/{sample}/tmp/binning/round_3/binette/temporary_files; fi
        if ls {loc}/{sample}/tmp/binning/round_3/binette/final_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/binette/final_bins/*.fa; fi

        if [ -d "{loc}/{sample}/tmp/binning/round_2/comebin/data_augmentation" ]; then rm -r {loc}/{sample}/tmp/binning/round_2/comebin/data_augmentation; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/cluster_res" ]; then rm -r {loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/cluster_res; fi
        if [ -f "{loc}/{sample}/tmp/binning/round_2/contigs.fasta.frag.ffn" ]; then rm {loc}/{sample}/tmp/binning/round_2/contigs.fasta.frag.*; fi
        if ls {loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/comebin_res_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/comebin_res_bins/*.fa; fi

        if ls {loc}/{sample}/tmp/binning/round_1/metabat2/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/metabat2/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_3/metabat2/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/metabat2/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_4/metabat2/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_4/metabat2/*.fa; fi

        if ls {loc}/{sample}/tmp/binning/round_1/vamb/bins/*.fna >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/vamb/bins/*.fna; fi
        if ls {loc}/{sample}/tmp/binning/round_3/vamb/bins/*.fna >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/vamb/bins/*.fna; fi

        if ls {loc}/{sample}/tmp/binning/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_2/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_2/bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_3/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_4/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_4/bins/*.fa; fi

        if ls {loc}/{sample}/tmp/binning/round_3/semibin/output_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/semibin/output_bins/*.fa; fi
        if [ "{binmode}" != "fast" ]; then if ls {loc}/{sample}/tmp/binning/round_1/semibin/output_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/semibin/output_bins/*.fa; fi; fi

        if [ "{binmode}" == "extended" ]; then
        if [ -f "{loc}/{sample}/tmp/binning/round_1/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/binning/round_1/contigs.fasta; fi
        if ls {loc}/{sample}/tmp/binning/singl/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/singl/bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_1/binette/final_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/binette/final_bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_1/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/bins/*.fa; fi; fi

        if [ -f "{loc}/{sample}/tmp/logs/usage_summary_stats.tsv" ]; then
        rsync {loc}/{sample}/tmp/logs/summary_mmlong2-lite.tsv {loc}/{sample}/tmp/logs/summary_mmlong2-lite_$(date +"%Y%m%d-%Hh%Mm%Ss").tsv 
        rm {loc}/{sample}/tmp/logs/usage_*.tsv; fi

        size_post=$(du -sb {loc}/{sample}/tmp | awk '{{printf "%.1f", $1/1024/1024/1024}}')
        files_post=$(find {loc}/{sample}/tmp -type f -printf '.' | wc -c)
        echo "Status after cleanup (MAG production): $files_post files and $size_post GB of storage" >> {output}

        else echo "Cleanup skipped for MAG production" > {output}; fi
        """
