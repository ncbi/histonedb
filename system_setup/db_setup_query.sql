CREATE DATABASE IF NOT EXISTS db_name;
DROP USER IF EXISTS 'db_user'@'localhost';
CREATE USER 'db_user'@'localhost' IDENTIFIED BY 'db_password';
GRANT ALL PRIVILEGES ON db_name.* TO 'db_user'@'localhost';
ALTER DATABASE db_name CHARACTER SET utf8 COLLATE utf8_general_ci;
