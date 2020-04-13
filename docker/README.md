# HistoneDB via virtual machine

Create project directory and directories where to mount project code and database

```
mkdir project_dir
cd project_dir
mkdir histonedb
mkdir db
```
Then you can run HistoneDB via two ways.

## Running via docker

- Run as a service in docker, this will run apache and attempt to start mysqld

```docker run --name histdb -d -p 8080:10080 -v project_dir/histonedb:/var/www/histonedb -v project_dir/db:/var/lib/mysql intbio/histonedb:0.0.1  ```

- Get into container and start db regeneration

```docker exec -it histdb bash```

- Next in reinit_histdb_local.sh adjust the database you would want to build HistoneDB from (swissprot, nr, etc.)

```bash db_gen.sh -mysql_db_reinit -histdb_reinit```

- To stop the container run

```docker stop histdb```

## Run in singularity 

- Build singularity container

```singularity build --sandbox cont docker://intbio/histonedb:0.0.1```

- Run apache on prot 10080 and attempt to start mysqld

```singularity instance start --writable --bind project_dir/histonedb:/var/www/histonedb,project_dir/db:/var/lib/mysql cont histdb```

- Get into container

```singularity shell instance://histdb```

- Regenerate db

```
apachectl start
cd /var/www
```

- Next in ```reinit_histdb_local.sh``` adjust the database you would want to build HistoneDB from (swissprot, nr, etc.)

```bash db_gen.sh -mysql_db_reinit -histdb_reinit```

- To stop the container run

```singularity instance stop histdb```


**IMPORTANT NOTE:
Looks like singularity is not completele isolated from the host machine kernel security modules, e.g. apparmor.
So it means on the host
/etc/apparmor.d/usr.sbin.mysqld should have
```
# Allow data files dir access
  /var/lib/mysql-files/ r,
  /var/lib/mysql-files/** rwk,
```

### For development
- Profiling
```python -m cProfile -s cumtime manage.py buildvariants```