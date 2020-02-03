#!/bin/bash

# for ubuntu:18.04

echo 'System preparing ...'
sh system_setup/full_setup.sh

echo 'Apache2 and WSGI configurating ...'
sh system_setup/apache_setup/apache_setup.sh

echo 'DataBase creating ...'
sudo mysql -u root --skip-password < /system_setup/db_setup_query.sql

echo 'Start initialization ...'
sh reinit_histdb_local.sh
echo 'Initialization complete. See loginfo in log/ directory.'
echo 'For reinitialization run file reinit_histdb_local.sh'
