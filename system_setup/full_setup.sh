# System Setup
yes Y | sh system_setup/sys_requirements_setup.txt
# Python Requirements Setup
yes Y | sh system_setup/py_requirements_setup.txt
# Create DB
sudo mysql -u root --skip-password < system_setup/db_setup_query.sql

# To cinfigure apache2 with mod_wsgi (already installed) please, see
# https://medium.com/faun/how-to-set-up-conda-virtual-environments-with-apache-mod-wsgi-flask-c2043711223e
sudo service apache2 restart