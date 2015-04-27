from django.core.management.base import NoArgsCommand
from optparse import make_option

from django.conf import settings
from djangophylocore.models import *

import os

class Command(NoArgsCommand):
    option_list = NoArgsCommand.option_list + (
        make_option('--verbose', '-v', action='store_true', dest='verbose', 
            help='Verbose operation'),
    )
    help = "Create parents arborescen into the database"
    
    requires_model_validation = True
    
    def generate_parents( self, taxa, parents, file ):
        index = 0
        results = []
        for parent in reversed(parents):
            self.rel_id += 1
            results.append( "%s|%s|%s|%s\n"% ( self.rel_id, taxa.id, parent.id, index ) )
            #ParentsRelation.objects.create( taxa = taxa , parent = parent, index = index )
            index += 1
        file.write( ''.join( results ) )
        parents2 = parents[:]
        parents2.append( taxa )
        for child in taxa.direct_children.all():
            self.generate_parents( child, parents2, file )

    def handle_noargs(self, **options):
        ParentsRelation.objects.all().delete()
        root, created = Taxa.objects.get_or_create( id=1, name = 'root' )
        parents = [root]
        localDir = os.path.dirname(__file__)
        absDir = os.path.join(os.getcwd(), localDir)
        path_dumps = os.path.join( absDir,'..','..', 'dumps' ) 
        dmp_file_path = os.path.join( path_dumps, "parentsrelation.dmp" )
        if not os.path.exists( dmp_file_path ):
            print "Generating parents..."
            open( dmp_file_path , "w" ).write( '' )
            file = open( os.path.join( path_dumps, "parentsrelation.dmp" ), 'a' )
            self.rel_id = 0
            for child in root.direct_children.all():
                self.generate_parents( child, parents, file )
            file.close()
        print "Loading parents into the database..."
        db_name = settings.DATABASE_NAME
        db_engine = settings.DATABASE_ENGINE
        if db_engine == 'sqlite3':
            cmd = "sqlite3 -separator '|' %s '.import %s djangophylocore_parentsrelation'" % ( 
              db_name, dmp_file_path )
        elif db_engine == 'mysql':
            cmd = """mysql -u %s -p%s %s -e "LOAD DATA LOCAL INFILE '%s' INTO TABLE djangophylocore_parentsrelation FIELDS TERMINATED BY '|';" """ % (
              settings.DATABASE_USER, settings.DATABASE_PASSWORD, db_name, dmp_file_path )
        print cmd
        os.system( cmd )

