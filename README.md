# HistoneDB
A database for all histone proteins in NR organized by their known non-allelic protein isoforms, called variants. This resource can be used to understand how changes in histone variants affect structure, complex formation, and nucleosome function. For more information, please read our [paper](manuscript/paper.md) (citation below).

The database can be accessed at https://histonedb.bioeng.ru

## Requirements ##

- Python 2.7
- Required python packages are specified in ```requirements.txt``` (use ```pip install -r requirements_full.txt```). They include:
  - Django 1.8.19, django-debug_toolbar, django-filter, django-filters 0.2.1
  - django-extensions 1.5.3
  - MySQL-python 1.2.5
  - BioPython 1.74
  - colour 0.1.5
  - pyparsing, matplotlib, seaborn, sklearn, networkx, more_itertools, numpy, pandas, scipy
  - psycopg2 2.7
- [HMMER 3.1b2](http://hmmer.janelia.org)
- [BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v 2.2.26
- [EMBOSS](http://emboss.sourceforge.net) v6.5.7
- [MUSCLE](http://www.drive5.com/muscle/) v3.8.31
- [ClustalW2](http://www.clustal.org/clustal2/) v2.1
All executables must be present in the bin dir of virual environment.

## Setup ##

If you want to test the server on your own machine, you must make sure have all of the dependencies listed above, or if you have Ubuntu 18.04 and greater you can pass the instrustions from [requirements](system_setup/REQUIREMENTS.md).

Then follow these steps:

1) Create MySQL database (see [database configurations](system_setup/DB_CONFIGURATION.md)), and store the login information in the file ```HistoneDB/NCBI_databse_info.txt```, which is formatted in the following way (key = value):
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
WARNING: This will by default download the entire nr database and classify all sequences in the nr. If you want to build the HistoneDB using a smaller database of proteins in FASTA format using the NR formatted header (">UNIQUE_ACCESSION|anything description [TAXONOMY_NAME]"), run the following command:

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

For deployment, you can use any system as you prefer, but we recomend using modwsgi with configuration described in our [configuration instructions](system_setup/CONFIGURATION.md).

## Cite ##

Coming soon.

## Acknowledgements ##

* Eli Draizen
* Alexey K. Shaytan
* Anna Panchenko
* Leonardo Marino-Ramirez
* David Landsman
* Paul Talbert
