# HistoneDB
A database for all histone proteins in NR organized by their known non-allelic protein isoforms, called variants. This resource can be used to understand how changes in histone variants affect structure, complex formation, and nucleosome function. For more information, please read our [paper](manuscript/paper.md) (citaion below).

The database can be accesed at http://www.ncbi.nlm.nih.gov/projects/histonedb/

## Requirements ##
- Python 2.7
- Flup 1.0.2, if using fastcgi
- Django 1.8
- django-extensions 1.5.3
- MySQL-python 1.2.5
- BioPython 1.65
- colour 0.1.1
- HMMER 3.1b2
- BLAST+
- EMBOSS

## Setup ##
If you want to test the server on your own machine, you must make sure have all of the dependies listed above and follow these steps.

1) Create MySQL database and store login information in HistoneDB/NCBI_database_info.py, which is formated in the follwing way:
```
#This file contains the user name and password for the HistoneDB 2.0
#Keep hidden
name = "DB NAME"
user = "DB USER"
password = "DB PASS"
host = "DB URL"
```

2) Build NCBI Taxonomy with djangophylocore
```
python manage.py buildncbi
python manage.py loadtaxonomy
python manage.py buildtaxonomytoc
```
3) Classify sequences in NR
```
python manage.py buildvariants
```
4) Build trees from seed sequences
```
python manage.py buildtrees
```
5) Build organism sunbursts for each variant
```
python manage.py buildsunburst
```
6) Build MSA and GFF seuqence features for variants
```
python manage.py buildseedinfo
```

## Run ##
You have saveral options to run the Django server. The easiest way is to run it through manage.py specifying a port (we use port 8080 in the example):
```
python manage.py runserver 8080
```
For deployment, we use FastCGI on the NCBI webservers. While this will be depracted in the next version Django, it is what NCBI allows. For more info, please read https://docs.djangoproject.com/en/1.8/howto/deployment/fastcgi/

## Cite ##
Coming soon.

## Acknoledgements ##
Eli Draizen
Alexey K. Shaytan
Anna Panchenko
Leonardu Marino-Ramirez
David Landsman
Paul Talbert
