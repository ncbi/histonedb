from django.core.management.base import NoArgsCommand
from optparse import make_option
from django.db import connection

import os, codecs
from django.conf import settings

if settings.TAXONOMY_ENGINE == "ncbi":
    from phylocore_ncbi import *
elif settings.TAXONOMY_ENGINE == "itis":
    from phylocore_itis import *

class Command(NoArgsCommand):
    option_list = NoArgsCommand.option_list + (
        make_option('--verbose', '-v', action='store_true', dest='verbose', 
            help='Verbose operation'),
    )
    help = "Load treebase collection data into database"
    
    requires_model_validation = True
    
    def handle_noargs(self, **options):
        localDir = os.path.dirname(__file__)
        absDir = os.path.join(os.getcwd(), localDir)
        verbose = options.get("verbose", False)
        if verbose:
            print "loading treebase, please wait, it can take a while..."
        treebase_collection = open( os.path.join( absDir,'..','..', 'treebase.nexus') ).read()
        TreeCollection.objects.create( source = treebase_collection, delimiter="_" ) #VR sept09 add delimiter
        if verbose:
            print "...done"
#        db_name = settings.DATABASE_NAME
#        db_engine = settings.DATABASE_ENGINE
#        map_dumps = {}
#        cursor = connection.cursor()
#        if settings.DATABASE_ENGINE == 'sqlite3':
#            cursor.execute( """ CREATE TABLE "djangophylocore_reltreecoltaxa1" (
#                "id" integer NOT NULL PRIMARY KEY,
#                "collection_id" integer NOT NULL REFERENCES "djangophylocore_treecollection" ("id"),
#                "tree_id" integer NOT NULL REFERENCES "djangophylocore_tree" ("id"),
#                "taxon_id" integer NULL REFERENCES "djangophylocore_taxonomy" ("id"),
#                "user_taxon_name" varchar(200) NULL
#            );""" )
#        elif settings.DATABASE_ENGINE == 'mysql':
#            cursor.execute( """ CREATE TABLE `djangophylocore_reltreecoltaxa1` (
#                `id` integer NOT NULL PRIMARY KEY,
#                `collection_id` integer NOT NULL REFERENCES `djangophylocore_treecollection` (`id`),
#                `tree_id` integer NOT NULL REFERENCES `djangophylocore_tree` (`id`),
#                `taxon_id` integer NULL REFERENCES `djangophylocore_taxonomy` (`id`),
#                `user_taxon_name` varchar(200) NULL
#            );""" )
#        for dump in os.listdir( path_dumps ):
#            if dump == ".svn": continue #XXX
#            name = os.path.splitext( dump )[0]
#            map_dumps[name] = os.path.join( path_dumps, dump )
#            if name in ['tree', 'treecollection']:
#                separator = '!'
#            else:
#                separator = '|'
#            if db_engine == 'sqlite3':
#                cmd = "sqlite3 -separator '%s' %s '.import %s djangophylocore_%s'" % (
#                  separator, db_name, map_dumps[name],  name) 
#            elif db_engine == 'mysql':
#                cmd = """mysql -u %s -p%s %s -e "LOAD DATA LOCAL INFILE '%s' INTO TABLE djangophylocore_%s FIELDS TERMINATED BY '%s';" """ % (
#                  settings.DATABASE_USER, settings.DATABASE_PASSWORD, db_name, map_dumps[name],  name, separator )
#            if verbose:
#                print cmd
#            os.system( cmd )
#            if settings.DATABASE_ENGINE == 'mysql':
#                cursor.execute( """CREATE INDEX djangophylocore_reltreecoltaxa1_taxon_id ON djangophylocore_reltreecoltaxa1 (taxon_id);""" )
#            cursor.close()
