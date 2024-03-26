from django.core.management.base import BaseCommand

import os
from django.conf import settings
from djangophylocore.models import *

try:
    import pickle as pickle
except:
    import pickle


class Command(BaseCommand):
    help = "Build taxonomy toc"

    requires_system_checks = [True]

    def handle(self, **options):
        d = {}
        taxas = Taxonomy.objects.all()
        for taxa in taxas.iterator():
            d[taxa.name] = taxa.id
        localDir = os.path.dirname(__file__)
        absDir = os.path.join(os.getcwd(), localDir)
        path = os.path.join(absDir, '..', '..')
        try:
            os.remove(os.path.join(path, 'taxonomy_toc_ncbi'))
        except:
            pass
        pickle.dump(d, open(os.path.join(path, 'taxonomy_toc_ncbi'), 'wb'))
