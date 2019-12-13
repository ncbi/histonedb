#!/bin/bash
#This script runs all the commands to initialize hisoneDB from scratch.
echo "Reinitializing HistDB from scratch"

#python manage.py flush << EOF
#yes
#EOF

#python manage.py migrate
#python manage.py sqlclear browse | python manage.py dbshell
#python manage.py sqlclear djangophylocore | python manage.py dbshell
#python manage.py migrate #закоментировала временно


#python manage.py buildncbi #закоментировала временно
#python manage.py loadtaxonomy #закоментировала временно
#python manage.py buildtaxonomytoc #закоментировала временно
python manage.py buildvariants -f
#python manage.py buildvariants -f --db nr_small
python manage.py buildtrees -f
python manage.py buildsunburst -f
python manage.py buildblastdb -f
python manage.py buildvariantinfo -f
python manage.py buildseedinfo -f

# проверить время, последняя строка в buildvariants 2019-12-11 14:26:35 tools.load_hmmsearch INFO     Loading variant: H2A.W
# в 18:03, или 15:03 по UTC время изменила
