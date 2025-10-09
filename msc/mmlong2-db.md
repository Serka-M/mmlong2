## Manual installation of databases used by mmlong2

Multiple databases are used by the mmlong2 workflow to characterize the recovered genomes.
Since there are various ways of setting up databases, code examples on how to download the relevant databases are presented below.
<br/>

The downloaded databases can be provided to mmlong2 by either using the database settings of the wrapper script (e.g. `mmlong2 --database_gtdb /path/to/gtdb/database`) or by modifying the `mmlong2-proc-config.yaml` configuration file to use the databases by default.
<br/>

#### Greengenes2 for 16S rRNA taxonomy
```
wget -O db_rrna.udb.gz https://zenodo.org/records/17174373/files/greengenes2_2024.09.udb.gz
gzip -d db_rrna.udb.gz
```

#### Metabuli database for contig taxonomy
```
wget -r -np -nH --cut-dirs=2 -R "index.html*" -P db_metabuli https://hulk.mmseqs.com/jaebeom/gtdb226db/
```

#### GTDB (Genome Taxonomy Database) for microbial genome taxonomy
```
wget -O db_gtdb.tar.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
tar -xvzf db_gtdb.tar.gz
```

#### GUNC database for genome QC
```
wget -O db_gunc.dmnd.gz https://swifter.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz
gzip -d db_gunc.dmnd.gz
```

#### Bakta databases for microbial genome annotation (requires [AMRFinderPlus](https://github.com/ncbi/amr/wiki))
```
wget -O db_bakta.tar.xz https://zenodo.org/records/14916843/files/db.tar.xz
tar -xJf db_bakta.tar.xz -C .
amrfinder_update --database db/amrfinderplus-db
```

[//]: # (Written by Mantas Sereika)
