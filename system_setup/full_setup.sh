sh system_setup/sys_requirements_setup.txt
sh system_setup/py_requirements_setup.txt
sudo mysql -u root -p < system_setup/db_setup_query.sql
