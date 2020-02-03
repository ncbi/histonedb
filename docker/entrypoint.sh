#!/bin/bash

db_create=false
#reinit_histdb=false

while [ -n "$1" ]
do
case "$1" in
-db) db_create=true ;;
#-re) reinit_histdb=true ;;
*) echo "$1 is not an option" ;;
esac
shift
done

if ! sudo mysql -u root --skip-password -e 'use histonedb' || $db_create ; then

  if sudo mysql -u root --skip-password -e 'use histonedb' ; then
    echo "Dropping the existing database histonedb ..."
    sudo mysql --skip-password --execute="DROP DATABASE histonedb;"
  fi

  cd /var/www/histonedb
  echo 'Database creating ...'
  sudo mysql -u root --skip-password < system_setup/db_setup_query.sql
  echo "Database histonedb created at /var/lib/mysql/histonedb inside container."

  echo 'Start initialization ...'
  sh reinit_histdb_local.sh
  echo 'Initialization complete. See loginfo in log/ directory.'
  echo 'For reinitialization run file reinit_histdb_local.sh'

else

  if sudo mysql -u root --skip-password -e 'use histonedb' ; then
    echo "Running the existing database histonedb localed at /var/lib/mysql/histonedb inside container."
  fi

fi
