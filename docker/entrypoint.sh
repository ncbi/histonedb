#!/bin/bash

#echo "Errorlog /dev/stderr" >> /etc/apache2/apache2.conf

easy_setup=false
without_setup=false

while [ -n "$1" ]
do
case "$1" in
-easy_setup) easy_setup=true ;;
-without_setup) without_setup=true ;;
*) echo "$1 is not an option" ;;
esac
shift
done

mysqld --initialize-insecure --user=mysql
#/usr/bin/mysqld_safe
mysqld_safe &

if ! $without_setup ; then

  if ! mysql -u root --skip-password -e 'use histonedb' > /dev/null 2>&1 || $easy_setup ; then

    if mysql -u root --skip-password -e 'use histonedb' > /dev/null 2>&1 ; then
      echo "Dropping the existing database histonedb ..."
      mysql --skip-password --execute="DROP DATABASE histonedb;"
    fi

    cd /var/www/histonedb
    echo 'Database creating ...'
    if mysql -u root --skip-password < system_setup/db_setup_query.sql ; then
      echo "Database histonedb created at /var/lib/mysql/histonedb inside container."
    else
      exit
    fi

    echo 'Start initialization ...'
#    sh reinit_histdb_local.sh > reinit.log 2>error.log
    sh reinit_histdb_local.sh
    echo 'Initialization complete. See loginfo in log/ directory.'
    echo 'For reinitialization run file reinit_histdb_local.sh'

  else

    if mysql -u root --skip-password -e 'use histonedb' > /dev/null 2>&1 ; then
      echo "Running the existing database histonedb localed at /var/lib/mysql/histonedb inside container."
    fi

  fi

else

  echo "Project started without initialization."

fi

#apachectl -DFOREGROUND

#a2enconf wsgi
service apache2 restart