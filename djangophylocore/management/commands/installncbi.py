from django.core.management.base import NoArgsCommand
from optparse import make_option

import os
from django.conf import settings


class Command(NoArgsCommand):
    option_list = NoArgsCommand.option_list + (
        make_option('--verbose', '-v', action='store_true', dest='verbose', 
            help='Verbose operation'),
   )
    help = "Load all taxonomy data into database"
    
    requires_model_validation = True
    
    def handle_noargs(self, **options):
        verbose = options.get("verbose", False)
        verbose_string = ''
        if verbose:
            verbose_string = '-v'
        if verbose:
            print("building ncbi dumps")
        os.system('python manage.py buildncbi %s' % verbose_string)
        if verbose:
            print("creating database")
        os.system('python manage.py reset_db --noinput')
        os.system('python manage.py syncdb')
        if verbose:
            print("loading taxonomy")
        os.system('python manage.py loadtaxonomy %s' % verbose_string)
        if verbose:
            print("loading treebase informations")
        os.system('python manage.py loadtreebase %s' % verbose_string)
        if verbose:
            print("building taxonomy toc")
        os.system('python manage.py buildtaxonomytoc')

