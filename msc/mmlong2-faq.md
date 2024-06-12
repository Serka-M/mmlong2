## Frequently asked questions about mmlong2

#### Will my sequenced dataset work with mmlong2?
* It is very situational whether mmlong2 is a good fit for your samples or project, but generally, the workflow is intended for highly complex metagenomes (e.g. soil, sewage sludge, human gut) and is not optimal for samples with very low microbial diversity (e.g. pure cultures, Zymo DNA Standard).
* Please also keep in mind that mmlong2 is a long-reads only workflow, designed to work with Nanopore (about 1 % read error rate) or with PacBio HiFi (about 0.1 % read error rate) datasets. Short read datasets can be used for mapping to improve genome recovery via differential coverage, but the workflow is not capable of short-read assembly.
* It also recommended that the input for mmlong2 is at least 1 GB of sequenced data for multiple prokaryotic organisms.

#### Is any preprocessing required for running mmlong2?
* It is highly recommended to perform read quality filtering (e.g. >= Phreq Q10 for Nanopore and >= Phred Q20 for PacBio HiFi) and adaptor, barcode trimming before perming read assembly with mmmlong2

#### I want to test mmlong2 with my samples, but I don't want to install over 100 Gb of databases
* If you are only interested in getting the MAGs, check out [mmlong2-lite](https://github.com/Serka-M/mmlong2-lite), which is a lightweight version of the pipeline with identical MAG recovery procedure and does not require large database installation.

*Comments:*
* The workflow assumes that the input reads have been quality-filtered and adapter/barcode sequences have been trimmed off.
* The workflow is long-read-based and requires either Nanopore or PacBio HiFi reads. It doesn't feature an Illumina-only mode.
* If the workflow crashes, it can be resumed by re-running the same command. Some of the intermediary files might have to be removed for compatibility.

[//]: # (Written by Mantas Sereika)
