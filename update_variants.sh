#!/bin/bash
python manage.py buildvariants -f --db=/web/public/data/projects/histonedb/nr 
python manage.py buildtrees -f
python manage.py buildsunburst -f
python manage.py buildblastdb -f
python manage.py buildvariantinfo -f
python manage.py buildseedinfo -f
