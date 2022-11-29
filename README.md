# mmlong2
Automated long-read metagenomics workflow, using either PacBio HiFi or Nanopore sequencing reads as input to generate characterized MAGs.
Mmlong2 is a continuation of [mmlong](https://github.com/SorenKarst/mmlong).

**Note:** At the moment, mmlong2 is an in-house pipeline at Aalborg University Bioservers. A distributable version of the pipeline is scheduled for a future release.

**Comments:**
* The workflow assumes that the input reads have been quality-filtered and adapter/barcode sequences have been trimmed off.
* The workflow is long-read-based and requires either Nanopore or PacBio HiFi reads. It doesn't feature an Illumina-only mode.
* If the workflow crashes, it can be resumed by re-running the same command.
* Medaka parallelisation is capped at 20 instances.
* MetaBinner and Vamb binners require a complex metagenome to run, and thus will exit with an error when reads for a low-complexity metagenome are provided.
* When the workflow is running with a large amount of threads (eg. 100), CheckM2 or GTDB-tk can get stuck. This can be fixed by running the workflow at 40 threads.
