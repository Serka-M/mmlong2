## Manual installation of databases used by mmlong2

Multiple databases are used by the mmlong2 workflow to characterize the recovered genomes.
Since there are various ways of setting up databases, code examples on how to download the relevant databases are presented below.
<br/>

The downloaded databases can be provided to mmlong2 by either using the database settings of the wrapper script (e.g. `mmlong2 --database_gtdb /path/to/gtdb/database`) or by modifying the `mmlong2-proc-config.yaml` configuration file to use the databases by default.
<br/>

#### Greengenes2 for 16S rRNA taxonomy
```
wget -O db_rrna_greengenes2-2024.09.udb.gz https://zenodo.org/records/17174373/files/greengenes2_2024.09.udb.gz
gzip -d db_rrna_greengenes2-2024.09.udb.gz
```

#### Metabuli database for contig taxonomy
```
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_nr_2023-05-10.tgz
tar -xvzf kaiju_db_nr_2023-05-10.tgz
```

#### GTDB (Genome Taxonomy Database) database for microbial genome taxonomy
```
wget -O gtdb.tar.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
tar -xvzf gtdb.tar.gz
```

#### GUNC database for genome QC
```
wget -O db_gunc.dmnd.gz https://swifter.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz
gzip -d db_gunc.dmnd.gz
```

#### Bakta databases for microbial genome annotation (requires [AMRFinderPlus](https://github.com/ncbi/amr/wiki))
```
wget https://zenodo.org/records/10522951/files/db.tar.gz
tar -xvzf db.tar.gz
amrfinder_update --database db/amrfinderplus-db
```

[//]: # (Written by Mantas Sereika)
