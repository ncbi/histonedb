from django.core.management.base import NoArgsCommand
from optparse import make_option

import os
from django.conf import settings
from djangophylocore.models import *
try:
    import cPickle as pickle
except:
    import pickle

class Command(NoArgsCommand):
    option_list = NoArgsCommand.option_list + (
#        make_option('--verbose', '-v', action='store_true', dest='verbose', 
#            help='Verbose operation'),
    )
    help = "Build taxonomy toc"
    
    requires_model_validation = True
    
    def handle_noargs(self, **options):
        d = {}
        taxas = Taxonomy.objects.all()
        for taxa in taxas.iterator():
            d[taxa.name] = taxa.id
        localDir = os.path.dirname(__file__)
        absDir = os.path.join(os.getcwd(), localDir)
        path = os.path.join( absDir,'..','..' ) 
        try:
            os.remove( os.path.join( path, 'taxonomy_toc_%s' % settings.TAXONOMY_ENGINE ) )
        except:
            pass
        pickle.dump( d, open( os.path.join( path, 'taxonomy_toc_%s' % settings.TAXONOMY_ENGINE ), 'w' ) )

