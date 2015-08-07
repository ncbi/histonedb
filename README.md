# HistoneDB
A database for all histone proteins in NR organized by their known non-allelic protein isoforms, called variants. This resource can be used to understand how changes in histone variants affect structure, complex formation, and nucleosome function. For more information, please read our [paper](manuscript/paper.md) (citation below).

The database can be accessed at http://www.ncbi.nlm.nih.gov/projects/histonedb/

## Requirements ##

- Python 2.7
- Flup 1.0.2, if using fastcgi
- Django 1.8, django-debug_toolbar, django-filter, django-filters(?)
- django-extensions 1.5.3
- MySQL-python 1.2.5
- BioPython 1.65
- colour 0.1.1
- pyparsing, matplotlib, seaborn, sklearn
- [HMMER 3.1b2](http://hmmer.janelia.org)
- [BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [EMBOSS](http://emboss.sourceforge.net)
- [MUSCLE](http://www.drive5.com/muscle/)
- [ClustalW2](http://www.clustal.org/clustal2/)

## Setup ##

If you want to test the server on your own machine, you must make sure have all of the dependencies listed above and follow these steps.

1) Create MySQL database, and store the login information in the dictionary NCBI_databse_info in HistoneDB/settings.py, which is formatted in the following way:
```
NCBI_database_info = {
    "name": "DB_NAME",
    "user": "DB_USER",
    "password": "DB_PASS",
    "host": "DB_HOST",
    "port": "DB_PORT",
    "SECRET_KEY": "DJANGO_SECRET_KEY",
    "URL":"/" #Root of site when accessed from a browser
}
```
If running on the mweb, these values will already be set as environment variables.

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

4) Classify sequences in NR

```
python manage.py buildvariants
```

5) Build trees from seed sequences

```
python manage.py buildtrees
```

6) Build organism sunbursts for each variant

```
python manage.py buildsunburst
```

7) Build MSA and GFF sequence features for variants

```
python manage.py buildseedinfo
```
8) Build Blast database for custom sequence analysis

```
python manage.py buildblastdb
```
9) Build variantinfo
```
python manage.py buildvariantinfo
```
## Update ##
If youe need to update or rebuild the database, e.g. if a new variant is discovered, you must rerun steps 4-8 adding the --force parameter after each command to make sure everything gets updated.

## Adding new variants ##
1) Collect representative sequences and create seed alignments using any methods you wish. Please read our [paper](manuscript/paper.md) for more info on how we collected the sequences and aligned them.

2) Place seed alignments in appropriate static directory:
```
static/browse/seeds/[CORE_HISTONE]/[VARIANT].fasta
```
3) Follow update instructions

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
* Leonardu Marino-Ramirez
* David Landsman
* Paul Talbert
