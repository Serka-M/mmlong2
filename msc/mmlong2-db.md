### Setting up databases used by mmlong2

Various databases are used by the mmlong2 workflow to characterize the recovered MAGs.
Since there's multiple ways of setting up databases, code examples on how to download the relevant databases are presented below.
<br/>

The downloaded databases can be used with mmlong2 by using the database settings of the wrapper script (e.g. `mmlong2 --gtdb /path/to/gtdb/database`) or by modifying the `mmlong2-proc-config.yaml` file to have the databases used by default.
<br/>

**MiDAS (Microbial Database for Activated Sludge) for 16S rRNA taxonomy**
```
wget -O MiDAS-v4.8.1-sintax.fasta https://www.midasfieldguide.org/files/downloads/taxonomies/SINTAX%20fa%20file%20MiDAS%204.8.1.fa
```
<br/>
