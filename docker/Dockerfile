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
RUN apt-get install -y pv

# Python setup
RUN apt-get install -y python python-pip
RUN yes | pip install virtualenv
#RUN virtualenv histdb_py27 && /bin/bash -c "source histdb_py27/bin/activate" # this is not really used
RUN yes | pip install --upgrade setuptools
RUN yes | pip install biopython==1.74 tqdm colour decorator dj-database-url Django==1.8.19 django-debug-toolbar==1.3.2 django-extensions django-filter==0.10.0 django-filters django-google-analytics-id gunicorn>=19.5.0 matplotlib==1.4.3 mock more-itertools MySQL-python networkx==1.10 nose numpy pandas pbr psycopg2==2.7.0 pyparsing python-dateutil pytz scikit-learn scipy seaborn==0.6.0 six sklearn sqlparse wsgiref

# DB tweaks
RUN mkdir /var/run/mysqld
#RUN mkfifo /var/run/mysqld/mysqld.sock
RUN chown -R mysql:mysql /var/run/mysqld
#RUN chown -R mysql:mysql /var/lib/mysql #Not needed for docker
#RUN chmod 755 /var/lib/mysql # Not needed for docker
RUN mkdir /var/www/mysql-data
RUN mkdir /var/www/histonedb

# Apache2 and wsgi configuration
COPY 000-default.conf 000-default.conf
RUN cp -f 000-default.conf /etc/apache2/sites-available && rm 000-default.conf
RUN sed -i 's/Listen/#Listen/' /etc/apache2/ports.conf

COPY entrypoint.sh entrypoint.sh
ENTRYPOINT ["bash", "/var/www/entrypoint.sh"]