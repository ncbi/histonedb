#!/bin/bash

mysqldatadir=/var/lib/mysql/
#/usr/bin/mysqld_safe
echo "Running mysql in the background"
mysqld_safe --datadir=$mysqldatadir --port=13306 &


cd /var/www/histonedb
echo 'Starting apache2'
apachectl -D FOREGROUND

#In singulatiry we can simply run 'apachectl start'