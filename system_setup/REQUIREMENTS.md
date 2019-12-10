# Requirements setup for HistoneDB 

To fully and correctly run the project, some system configurations and dependencies are important. To setup Ubuntu machine system, python packages and MySQL database follow steps below:

1) Create a virtualenv (for correct setup, please, do not use anaconda) in your working directory and activate it
```
virtualenv histdb_py27
source histdb_py27/bin/activate
``` 

2) For system setup run script system_setup/full_setup.sh
```
sh system_setup/full_setup.sh
```
This command will install all system dependencies (system_setup/sys_requirements_setup.txt). Then there will be installed all python dependencies in activated virtualenv (system_setup/py_requirements_setup.txt). The step includes installing mysql-client.

3) Create new MySQL user and database. Change file system_setup/db_setup_query.sql by replacing db_name, db_user and password
```
CREATE DATABASE IF NOT EXISTS db_name;
CREATE USER 'db_user'@'localhost' IDENTIFIED BY 'password';
GRANT ALL PRIVILEGES ON histonedb . * TO 'histonedb_user'@'localhost';
ALTER DATABASE histonedb CHARACTER SET utf8 COLLATE utf8_general_ci;
```
And run the script
```
sudo mysql -u root --skip-password < system_setup/db_setup_query.sql
```

Configure timeout for mysql session
```
sudo sed -i '$a wait_timeout = 31536000' /etc/mysql/mysql.conf.d/mysqld.cnf
sudo sed -i '$a interactive_timeout = 31536000' /etc/mysql/mysql.conf.d/mysqld.cnf
```
