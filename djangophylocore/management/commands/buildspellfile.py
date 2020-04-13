from django.core.management.base import NoArgsCommand
from optparse import make_option

import os
from django.conf import settings
from server.phylocore_models import *

class Command(NoArgsCommand):
    option_list = NoArgsCommand.option_list + (
#        make_option('--verbose', '-v', action='store_true', dest='verbose', 
#            help='Verbose operation'),
   )
    help = "Build spell file for suggestions"
    
    requires_system_checks = True
    
    def handle_noargs(self, **options):
        d = {}
        taxas = Taxonomy.objects.all().iterator
        result = []
        for taxa in taxas():
            result.append(taxa.name.encode("iso-8859-15"))
        localDir = os.path.dirname(__file__)
        absDir = os.path.join(os.getcwd(), localDir)
        path = os.path.join(absDir, "..", "..", "lib")
        open(os.path.join(path, "spell_file_ncbi.txt"), "w").write("\n".join(result))

