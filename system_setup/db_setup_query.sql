CREATE DATABASE IF NOT EXISTS db;
CREATE USER 'db_user'@'localhost' IDENTIFIED BY 'pwd';
GRANT ALL PRIVILEGES ON histonedb . * TO 'db_user'@'localhost';
ALTER DATABASE db CHARACTER SET utf8 COLLATE utf8_general_ci;
