#!/bin/bash

db_create=false
reinit_histdb=false

while [ -n "$1" ]
do
case "$1" in
-db) db_create=true ;;
-re) reinit_histdb=true ;;
*) echo "$1 is not an option" ;;
esac
shift
done

if [ -d /var/lib/mysql/histonedb || $db_create ] ; then
    # Do Stuff ...
fi
