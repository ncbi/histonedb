#!/bin/bash

#echo "Errorlog /dev/stderr" >> /etc/apache2/apache2.conf
mysqldatadir=/var/lib/mysql/
easy_setup=false
without_setup=false
db_reinit=false

while [ -n "$1" ]
do
case "$1" in
-easy_setup) easy_setup=true ;;
-without_setup) without_setup=true ;;
-db_reinit) db_reinit=true ;;
*) echo "$1 is not an option" ;;
esac
shift
done

if $db_reinit; then
echo "Deleting database files"
rm -rf $mysqldatadir/*
fi

if [ "$(ls -A $mysqldatadir)" ]; then
     echo "$mysqldatadir is not Empty, no need to initialize mysql"
else
    echo "$mysqldatadir is Empty, initializing mysql"
    mysqld --initialize-insecure --datadir=$mysqldatadir --user=mysql
fi

#/usr/bin/mysqld_safe
echo "Running mysql in the background"
mysqld_safe --datadir=$mysqldatadir --port=13306 &

# if ! $without_setup ; then

  # if ! mysql -u root --skip-password -e 'use histonedb' > /dev/null 2>&1 || $easy_setup ; then

  #   if mysql -u root --skip-password -e 'use histonedb' > /dev/null 2>&1 ; then
  #     echo "Dropping the existing database histonedb ..."
  #     mysql --skip-password --execute="DROP DATABASE histonedb;"
  #   fi
if ! $without_setup; then
cd /var/www/histonedb
echo 'Database creating ...'
echo 'Allow 5 secs for database initialization ...'
sleep 5;
mysql -u root --skip-password < system_setup/db_setup_query.sql 
echo "Database  created "

echo 'Start HistoneDB initialization ...'
#    sh reinit_histdb_local.sh > reinit.log 2>error.log
bash reinit_histdb_local.sh
echo 'Initialization complete. See loginfo in log/ directory.'
fi
# echo 'For reinitialization run file reinit_histdb_local.sh'

#   else

#     if mysql -u root --skip-password -e 'use histonedb' > /dev/null 2>&1 ; then
#       echo "Running the existing database histonedb localed at /var/lib/mysql/histonedb inside container."
#     fi

#   fi

# else

#   echo "Project started without initialization."

# fi
echo 'Starting apache2'
apachectl -D FOREGROUND

#a2enconf wsgi
#service apache2 restart