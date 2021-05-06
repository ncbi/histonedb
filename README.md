# HistoneDB
A database for all histone proteins in NR organized by their known non-allelic protein isoforms, called variants. This resource can be used to understand how changes in histone variants affect structure, complex formation, and nucleosome function. For more information, please read our [paper](manuscript/paper.md) (citation below).

The database can be accessed at http://www.ncbi.nlm.nih.gov/projects/HistoneDB2.0/

## Requirements ##

- Python 3.8
- Required python packages are specified in requirements/ (use "pip install -r requirements/full.txt"). They include:
- [HMMER 3.1b2](http://hmmer.janelia.org)
- [BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v 2.2.26
- [EMBOSS](http://emboss.sourceforge.net) v6.5.7
- [MUSCLE](http://www.drive5.com/muscle/) v3.8.31
- [ClustalW2](http://www.clustal.org/clustal2/) v2.1
All executables must be present in the bin dir of virual environment.

## Setup ##

If you want to test the server on your own machine, you must make sure have all of the dependencies listed above and follow these steps.

1) Create MySQL database, and store the login information in the file  HistoneDB/NCBI_databse_info.txt, which is formatted in the following way (key = value):
```
name = DB_NAME
user = DB_USER
password = DB_PASS
host = DB_HOST
port = DB_PORT
SECRET_KEY = DJANGO_SECRET_KEY

#Reference another file with same paramters, useful if it needs be hidden
file = /path/to/other/file.txt

#Root of site when accessed from a browser
STATIC_URL = /static/ 
```
If running on the mweb, these values will already be set.
Finally, make sure the database has the correct charset:
```
ALTER DATABASE DB_NAME CHARACTER SET utf8 COLLATE utf8_general_ci
```

2) Migrate Django models into database

```
python manage.py migrate
```

3) Build NCBI Taxonomy with djangophylocore

```
python manage.py buildncbi
python manage.py loadtaxonomy
python manage.py buildtaxonomytoc
```
WARNING: This will download the entire NCBI taxonomy database and load into the database, which can take a long time.

4) Classify sequences in NR

```
python manage.py buildvariants
```
WARNING: This will by default download the entire nr database and classify all sequences in the nr. If you want to build the HistoneDB using a smaller database of proteins in FASTA format using the NR formatted header (">gi|UNIQUE_GI|anything description [TAXONOMY_NAME]"), run the following command:

```
python manage.py buildvariants --db small_database.fasta
```

5) Build trees from seed sequences

```
python manage.py buildtrees
```

6) Build organism sunbursts for each variant

```
python manage.py buildsunburst
```

7) Build Blast database for custom sequence analysis

```
python manage.py buildblastdb
```
8) Build variantinfo
```
python manage.py buildvariantinfo
```

9) Build GFF sequence features for variants

```
python manage.py buildseedinfo
```
## Update ##
If youe need to update or rebuild the database, e.g. if a new variant is discovered, you must rerun steps 4-8 adding the --force parameter after each command to make sure everything gets updated.

## Adding new variants ##
1) Collect representative sequences of the new variant and create seed alignments using any method you wish. Please read our [paper](manuscript/paper.md) for more info on how we collected the sequences and aligned them.

2) Place seed alignments in appropriate static directory:
```
static/browse/seeds/[HISTONE_TYPE]/[VARIANT].fasta
```
3) Follow update instructions

4) If this is a new variant, please let us know by creating a pull request, a new issue (enhancement), or emailing us.

## Run ##

You have several options to run the Django server. The easiest way is to run it through `manage.py`, specifying a port (we use port 8080 in the example):

```
python manage.py runserver 8080
```

For deployment, we use FastCGI on the NCBI webservers. While this will be deprecated in the next version of Django, it is what NCBI allows. For more info, please read https://docs.djangoproject.com/en/1.8/howto/deployment/fastcgi/

## Cite ##

Coming soon.

## Acknowledgements ##

* Eli Draizen
* Alexey K. Shaytan
* Anna Panchenko
* Leonardo Marino-Ramirez
* David Landsman
* Paul Talbert
