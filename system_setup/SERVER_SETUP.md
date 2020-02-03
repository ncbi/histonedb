# HistoneDB setup instructions

HistoneDB can be run simple from file ```run.sh```. This script contains some steps:

1) Install full list of server dependencies by the script ```system_setup/full_setup.sh```. See [requirements](REQUIREMENTS.md).

2) Apache2 and WSGI configuration. More info described in [configuration instructions](CONFIGURATION.md).

3) Database creation (see [database configurations](DB_CONFIGURATION.md)).

4) Project setup.

## Project setup ##

Starting the project include following steps.

1) Storing the login information in the file ```HistoneDB/NCBI_databse_info.txt```, which is formatted in the following way (key = value):
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

2) Migrations of Django models into database:

```
python manage.py migrate
```

3) Building NCBI Taxonomy with djangophylocore:

```
python manage.py buildncbi
python manage.py loadtaxonomy
python manage.py buildtaxonomytoc
```
WARNING: This will download the entire NCBI taxonomy database and load into the database, which can take a long time.

4) Classification of sequences in NR:

```
python manage.py buildvariants
```
WARNING: This will by default download the entire nr database and classify all sequences in the nr. If you want to build the HistoneDB using a smaller database of proteins in FASTA format using the NR formatted header (">UNIQUE_ACCESSION|anything description [TAXONOMY_NAME]"), run the following command:

```
python manage.py buildvariants --db small_database.fasta
```

5) Building trees from seed sequences:

```
python manage.py buildtrees
```

6) Building organism sunbursts for each variant:

```
python manage.py buildsunburst
```

7) Building Blast database for custom sequence analysis:

```
python manage.py buildblastdb
```
8) Building variantinfo:
```
python manage.py buildvariantinfo
```

9) Building GFF sequence features for variants:

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
