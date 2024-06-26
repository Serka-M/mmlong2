## Manual installation of databases used by mmlong2

Various databases are used by the mmlong2 workflow to characterize the recovered genomes.
Since there are multiple ways of setting up databases, code examples on how to download the relevant databases are presented below.
<br/>

The downloaded databases can be provided to mmlong2 by either using the database settings of the wrapper script (e.g. `mmlong2 --database_gtdb /path/to/gtdb/database`) or by modifying the `mmlong2-proc-config.yaml` configuration file to use the databases by default.
<br/>

#### MiDAS (Microbial Database for Activated Sludge) for 16S rRNA taxonomy
```
wget -O MiDAS-v5.3-sintax.fasta https://www.midasfieldguide.org/files/downloads/taxonomies/SINTAX.fa%20MiDAS%205.3.fa
```

#### SILVA database for 16S rRNA taxonomy
```
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
gzip -d SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
sed 's/^ *>[^ ]* />/' -i SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta
```

#### Kaiju database for contig taxonomy
```
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_nr_2023-05-10.tgz
tar -xvzf kaiju_db_nr_2023-05-10.tgz
```

#### GUNC database for contig taxonomy
```
wget https://swifter.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz
gzip -d gunc_db_progenomes2.1.dmnd.gz
```

#### GTDB (Genome Taxonomy Database) database for microbial genome taxonomy
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
tar -xvzf gtdbtk_r220_data.tar.gz
```

#### Bakta databases for microbial genome annotation (requires [AMRFinderPlus](https://github.com/ncbi/amr/wiki))
```
wget https://zenodo.org/records/10522951/files/db.tar.gz
tar -xvzf db.tar.gz
amrfinder_update --database db/amrfinderplus-db
```

[//]: # (Written by Mantas Sereika)
