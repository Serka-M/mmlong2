## Frequently asked questions about mmlong2

#### Will my sequenced dataset work with mmlong2?
* It is very situational whether mmlong2 is a good fit for your samples or project, but generally, the workflow is intended for highly complex metagenomes (e.g. soil, sewage sludge, human gut) and is not optimal for samples with very low microbial diversity (e.g. pure cultures, Zymo DNA Standard).
* Please keep in mind that mmlong2 is a long-reads only workflow, designed to work with Nanopore (about 1 % read error rate) or with PacBio HiFi (about 0.1 % read error rate) datasets. Short-read datasets can be used for mapping to improve genome recovery via differential coverage, but the workflow is not designed for short-read metagenomic assembly.
* It also recommended that the input for mmlong2 is at least 1 GB of sequenced data with multiple prokaryotic organisms.

#### Is any data pre-processing required for running mmlong2?
* It is highly recommended to perform read quality filtering (e.g. remove reads with less than Phred Q10 for Nanopore and Phred Q20 for PacBio HiFi as well as short-reads) and adaptor, barcode trimming before performing read assembly with mmlong2.

#### I want to test mmlong2 with my samples, but I don't want to install over 100 Gb of databases
* If you are only interested in getting the MAGs, check out [mmlong2-lite](https://github.com/Serka-M/mmlong2-lite), which is a lightweight version of the pipeline with identical MAG recovery procedure and does not require large database installation.

#### What to do when the mmlong2 workflow crashes?
* If the workflow crashes, it can usually be resumed by re-running the workflow with the exact same code. 
* If you want to resume the workflow from a new installation of mmlong2, it is highly recommended to first run the workflow with the `--touch` option to mark completed files against deletion.

#### What about eukaryotic or viral genomes?
* At the moment, mmlong2 is not designed to perform genome recovery of viruses or eukaryotes. 
* Expansion of the binning features is planned for future releases.

[//]: # (Written by Mantas Sereika)