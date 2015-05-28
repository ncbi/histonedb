from django.core.management.base import BaseCommand, CommandError
from browse.models import *
from tools import load_hmmsearch
from subprocess import Popen
import os

class Command(BaseCommand):
    help = 'Build the HistoneDB by training HMMs with seed sequences found in seeds directory in the top directory of thi project and using those HMMs to search the NR database.'

    def add_arguments(self, parser):
        parser.add_argument('-v', '--verbose' help="Print info")

    def handle(self, *args, **options):
        """
        1. Build Variants
        2. Build Trees
        3. Build Sunbursts
        """