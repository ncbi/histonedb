FROM ubuntu:18.04

LABEL "repository"="https://github.com/intbio/histonedb"
LABEL maintainer="l.singh@intbio.org"

WORKDIR /var/www

# System setup
RUN apt-get clean && apt-get update
RUN apt-get install -y git libmysqlclient-dev mysql-client libpq-dev python-dev muscle clustalw emboss hmmer ncbi-blast+ libpng-dev libfreetype6 libfreetype6-dev pkg-config libblas-dev liblapack-dev gfortran
RUN apt-get install -y apache2 apache2-dev mysql-server
RUN sed -i '$a wait_timeout = 31536000' /etc/mysql/mysql.conf.d/mysqld.cnf && sed -i '$a interactive_timeout = 31536000' /etc/mysql/mysql.conf.d/mysqld.cnf
RUN apt-get install -y libapache2-mod-wsgi
#RUN service mysql restart

# Python setup
RUN apt-get install -y python python-pip
RUN yes | pip install virtualenv
RUN virtualenv histdb_py27 && /bin/bash -c "source histdb_py27/bin/activate"
RUN yes | pip install --upgrade setuptools
RUN yes | pip install biopython==1.74 colour decorator dj-database-url Django==1.8.19 django-debug-toolbar==1.3.2 django-extensions django-filter==0.10.0 django-filters django-google-analytics-id gunicorn>=19.5.0 matplotlib==1.4.3 mock more-itertools MySQL-python networkx==1.10 nose numpy pandas pbr psycopg2==2.7.0 pyparsing python-dateutil pytz scikit-learn scipy seaborn==0.6.0 six sklearn sqlparse wsgiref

# DB creation
#RUN mysql -u root --skip-password < /app/db_setup_query.sql
#RUN mysql -u root -e "CREATE DATABASE IF NOT EXISTS histonedb;" && mysql -u root -e "CREATE USER 'histonedb_user'@'localhost' IDENTIFIED BY '468725';" && mysql -u root -e "GRANT ALL PRIVILEGES ON histonedb . * TO 'histonedb_user'@'localhost';" && mysql -u root -e "ALTER DATABASE histonedb CHARACTER SET utf8 COLLATE utf8_general_ci;"

#WORKDIR /var/www

# Apache2 and wsgi configuration
COPY wsgi.conf wsgi.conf
RUN cp -f wsgi.conf /etc/apache2/mods-available && rm wsgi.conf

# histonedb/docker/
#COPY entrypoint.sh entrypoint.sh
#ENTRYPOINT ["bash", "entrypoint.sh"]