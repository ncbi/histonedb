CREATE DATABASE IF NOT EXISTS histonedb;
CREATE USER 'histonedb_user'@'localhost' IDENTIFIED BY '468725';
GRANT ALL PRIVILEGES ON histonedb . * TO 'histonedb_user'@'localhost';
ALTER DATABASE histonedb CHARACTER SET utf8 COLLATE utf8_general_ci;
