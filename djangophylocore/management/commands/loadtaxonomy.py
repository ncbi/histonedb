from django.core.management.base import NoArgsCommand
from optparse import make_option

import os
from django.conf import settings


class Command(NoArgsCommand):
    help = "Load all taxonomy data into database"
    
    requires_system_checks = True
    
    def handle_noargs(self, **options):
        localDir = os.path.dirname(__file__)
        absDir = os.path.join(os.getcwd(), localDir)
        verbose = options.get("verbose", True)
        if verbose:
            print "loading taxonomy, please wait, it can take a while..."
        path_dumps = os.path.join( absDir,'..','..', 'dumps' ) 
        db_name = settings.DATABASES["default"]["NAME"]
        db_engine = settings.DATABASES["default"]["ENGINE"]
        map_dumps = {}
        # print os.listdir( path_dumps )

        for dump in os.listdir( path_dumps ):
        ###########
            if dump == '.svn': continue # XXX
            name = os.path.splitext( dump )[0]
            map_dumps[name] = os.path.join( path_dumps, dump )
            if db_engine.endswith('sqlite3'):
                cmd = "sqlite3 -separator '|' %s '.import %s server_%s'" % (
                  db_name, map_dumps[name],  name) 
            elif db_engine.endswith('mysql'):
                # cmd = """mysql --local-infile -h %s -u %s -p%s %s -e "SET FOREIGN_KEY_CHECKS=0; LOAD DATA LOCAL INFILE '%s' INTO TABLE djangophylocore_%s FIELDS TERMINATED BY '|';" """ % (
                  # settings.DATABASES["default"]["HOST"], settings.DATABASES["default"]["USER"], settings.DATABASES["default"]["PASSWORD"], db_name, map_dumps[name],  name )
                cmd = """cat %s | pv | mysql --local-infile -h %s -u %s -p%s %s -e "SET FOREIGN_KEY_CHECKS=0; LOAD DATA LOCAL INFILE '/dev/stdin' INTO TABLE djangophylocore_%s FIELDS TERMINATED BY '|';" """ % (
                  map_dumps[name], settings.DATABASES["default"]["HOST"], settings.DATABASES["default"]["USER"], settings.DATABASES["default"]["PASSWORD"], db_name, name )
            
            elif db_engine.endswith('psycopg2'):
                cmd = """psql -h %s -U %s -d %s -c "ALTER TABLE djangophylocore_%s DISABLE TRIGGER ALL; COPY djangophylocore_%s FROM '%s' WITH DELIMITER AS '|'; ALTER TABLE djangophylocore_%s ENABLE TRIGGER ALL; " """ % (
                  settings.DATABASES["default"]["HOST"],settings.DATABASES["default"]["USER"], db_name, name, name, map_dumps[name], name)
            if verbose:
                print cmd
            os.system( cmd )
#        if not os.path.exists( 
#          os.path.join( path_dumps, 'parentsrelation.dmp' ) ):
#            if verbose:
#                print "no parent relations found, generating..."
#            os.system( 'python manage.py generateparents' )

