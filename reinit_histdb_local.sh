#!/bin/bash
#This script runs all the commands to initialize hisoneDB from scratch.
echo "Reinitializing HistDB from scratch"

#python manage.py flush << EOF
#yes
#EOF

#python manage.py migrate
#python manage.py sqlclear browse | python manage.py dbshell
#python manage.py sqlclear djangophylocore | python manage.py dbshell
python manage.py migrate


python manage.py buildncbi
python manage.py loadtaxonomy
python manage.py buildtaxonomytoc
python manage.py buildvariants -f --db nr_small
python manage.py buildtrees
python manage.py buildsunburst
#python manage.py buildblastdb
#python manage.py buildvariantinfo
#python manage.py buildseedinfo
