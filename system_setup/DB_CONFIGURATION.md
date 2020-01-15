# Create an configure MySQL database for HistoneDB 

Create MySQL database. Change file ```system_setup/db_setup_query.sql``` by replacing db_name, db_user and password
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
