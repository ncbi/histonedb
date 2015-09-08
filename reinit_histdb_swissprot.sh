#!/bin/bash
#This script runs all the commands to initialize hisoneDB from scratch.
echo "Reinitializing HistDB from scratch"
find static/browse/seeds -name "*.gff" | xargs rm
find static/browse/seeds -name "*.pdf" | xargs rm

python manage.py flush << EOF
yes
EOF

python manage.py migrate


python manage.py buildncbi
python manage.py loadtaxonomy
python manage.py buildtaxonomytoc
python manage.py buildvariants -f --db swissprot 
python manage.py buildtrees -f
python manage.py buildsunburst -f
python manage.py buildblastdb -f
python manage.py buildvariantinfo -f
python manage.py buildseedinfo -f
