# System Setup
yes Y | sh system_setup/sys_requirements_setup.txt
# Python Requirements Setup
yes Y | sh system_setup/py_requirements_setup.txt
# Create DB
#sudo mysql -u root --skip-password < system_setup/db_setup_query.sql
#sudo sed -i '$a wait_timeout = 31536000' /etc/mysql/mysql.conf.d/mysqld.cnf
#sudo sed -i '$a interactive_timeout = 31536000' /etc/mysql/mysql.conf.d/mysqld.cnf

sudo service apache2 restart