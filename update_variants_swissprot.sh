#!/bin/bash
python manage.py migrate

python manage.py buildvariants -f --db swissprot
python manage.py buildtrees -f
python manage.py buildsunburst -f
python manage.py buildblastdb -f
python manage.py buildvariantinfo -f
python manage.py buildseedinfo -f
