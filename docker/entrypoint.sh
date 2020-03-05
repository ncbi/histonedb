#!/bin/bash

cd /var/www/histonedb
echo 'Starting apache2'
apachectl -D FOREGROUND

#In singulatiry we can simply run 'apachectl start'