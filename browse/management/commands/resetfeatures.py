from django.core.management.base import BaseCommand, CommandError
from browse.models import Sequence
from tools.load_hmmsearch import update_features

class Command(BaseCommand):
    help = 'Reset sequence features'

    def add_arguments(self, parser):
        pass
        
    def handle(self, *args, **options):
        for s in Sequence.objects.all():
            update_features(s, variant_model=s.variant)
