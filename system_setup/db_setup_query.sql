CREATE DATABASE IF NOT EXISTS db;
DROP USER IF EXISTS 'db_user'@'localhost';
CREATE USER 'db_user'@'localhost' IDENTIFIED BY 'pwd';
GRANT ALL PRIVILEGES ON db . * TO 'db_user'@'localhost';
ALTER DATABASE db CHARACTER SET utf8 COLLATE utf8_general_ci;
