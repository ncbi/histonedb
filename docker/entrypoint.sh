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

apt-get install -y mysql-server
#sed -i '$a wait_timeout = 31536000' /etc/mysql/mysql.conf.d/mysqld.cnf && sed -i '$a interactive_timeout = 31536000' /etc/mysql/mysql.conf.d/mysqld.cnf
#/usr/bin/mysqld_safe
service mysql restart

if ! $without_setup ; then

  if ! mysql -u root --skip-password -e 'use histonedb' || $easy_setup ; then

    if mysql -u root --skip-password -e 'use histonedb' ; then
      echo "Dropping the existing database histonedb ..."
      mysql --skip-password --execute="DROP DATABASE histonedb;"
    fi

    cd /var/www/histonedb
    echo 'Database creating ...'
    mysql -u root --skip-password < system_setup/db_setup_query.sql
    echo "Database histonedb created at /var/lib/mysql/histonedb inside container."

    echo 'Start initialization ...'
    sh reinit_histdb_local.sh
    echo 'Initialization complete. See loginfo in log/ directory.'
    echo 'For reinitialization run file reinit_histdb_local.sh'

  else

    if mysql -u root --skip-password -e 'use histonedb' ; then
      echo "Running the existing database histonedb localed at /var/lib/mysql/histonedb inside container."
    fi

  fi

else

  echo "Project started without initialization."

fi

apachectl -DFOREGROUND

#a2enconf wsgi
#service apache2 restart