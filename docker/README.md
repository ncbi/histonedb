### This is how we need to do it now
```docker image build -t intbio/histonedb:0.0.1 .```
```docker push intbio/histonedb:0.0.1```

#### Running via docker
- Run as a service in docker, this will run apache and attempt to start mysqld
```docker run --name histdb -d -p 8080:10080 -v /Users/alexsha/work_HD/histonedb:/var/www/histonedb -v /Users/alexsha/junk/db:/var/lib/mysql intbio/histonedb:0.0.1  ```
- Get into container and start db regeneration
```docker exec -it histdb bash```
- Next in reinit_histdb_local.sh adjust the database you would want to build HistoneDB from (swissprot, nr, etc.)

```bash db_gen.sh -mysql_db_reinit -histdb_reinit```

```docker stop histdb```

#### Run in singularity 

```singularity build --sandbox cont docker://intbio/histonedb:0.0.1```

- This will run apache on prot 10080 and attempt to start mysqld
```singularity instance start --writable --bind /mnt/ramdisk/hdb/histonedb:/var/www/histonedb,/mnt/ramdisk/hdb/db:/var/lib/mysql cont histdb```
```singularity instance start --writable --bind /home/_scratch/hdb/histonedb:/var/www/histonedb,/home/_scratch/hdb/db:/var/lib/mysql cont histdb```

- Regenerate db

```singularity shell instance://histdb```

```apachectl start```
```cd /var/www```
- Next in reinit_histdb_local.sh adjust the database you would want to build HistoneDB from (swissprot, nr, etc.)

```bash db_gen.sh -mysql_db_reinit -histdb_reinit```

```singularity instance stop histdb```

### For development
- Profiling
```python -m cProfile -s cumtime manage.py buildvariants```
