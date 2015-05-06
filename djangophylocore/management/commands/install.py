from django.core.management.base import NoArgsCommand
from optparse import make_option

import os
from django.conf import settings


class Command(NoArgsCommand):
    option_list = NoArgsCommand.option_list + (
        make_option('--verbose', '-v', action='store_true', dest='verbose', 
            help='Verbose behavior'),
        #make_option('--taxonomy', '-t', dest='taxonomy', default = 'ncbi',
        #    help='Taxonomy to use: itis or ncbi (default ncbi)'),
    )
    help = "Install base information into the database"
    
    requires_model_validation = True
    
    def handle_noargs(self, **options):
        verbose = options.get("verbose", False)
        #taxonomy = options.get("taxonomy", 'ncbi')
        taxonomy = settings.TAXONOMY_ENGINE
        assert taxonomy in ['ncbi', 'itis'], "taxonomy supported : itis or ncbi"
        verbose_string = ''
        if verbose:
            verbose_string = '-v'
        if verbose:
            print "building %s dumps" % taxonomy
        os.system( 'python manage.py build%s %s' % (taxonomy, verbose_string) )
        if verbose:
            print "creating database"
        os.system( 'python manage.py reset_db --noinput' )
        os.system( 'python manage.py syncdb' )
        if verbose:
            print "loading taxonomy"
        os.system( 'python manage.py loadtaxonomy %s' % verbose_string )
        if verbose:
            print "building taxonomy toc"
        os.system( 'python manage.py buildtaxonomytoc' )
        if verbose:
            print "building suggestion file"
        os.system( 'python manage.py buildspellfile' )
        if verbose:
            print "loading treebase informations"
        os.system( 'python manage.py loadtreebase %s' % verbose_string )

