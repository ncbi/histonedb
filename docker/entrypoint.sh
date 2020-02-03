#!/bin/bash

#db_create=false
#reinit_histdb=false
#
#while [ -n "$1" ]
#do
#case "$1" in
#-db) db_create=true ;;
#-re) reinit_histdb=true ;;
#*) echo "$1 is not an option" ;;
#esac
#shift
#done

if [ -d /var/lib/mysql/histonedb ] ; then
  echo "Running the existing database localed at /var/lib/mysql/histonedb inside container."
else
  cd /var/www/histonedb
  echo 'Database creating ...'
  sudo mysql -u root --skip-password < system_setup/db_setup_query.sql
fi
