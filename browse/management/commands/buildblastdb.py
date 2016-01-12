from django.core.management.base import BaseCommand, CommandError
from browse.models import *
from django.conf import settings

from Bio import SeqIO

import os
import sys
import subprocess

class Command(BaseCommand):
    help = 'Build the blast database used to analyze custom sequences.'

    def add_arguments(self, parser):
        parser.add_argument(
            "-f", 
            "--force", 
            default=False, 
            action="store_true", 
            help="Force the recreation of the blast db.")

    def handle(self, *args, **options):
        """Create new BLAST databse with seuqences in the HistoneDB."""
        force = options["force"]
    
        seqs_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "HistoneDB_sequences.fa")
        if not os.path.isfile(seqs_file) or force:
            with open(seqs_file, "w") as seqs:
                for s in Sequence.objects.filter(reviewed=True): #here we restrict the blast DB to reviewed seqs
                    SeqIO.write(s.to_biopython(ungap=True), seqs, "fasta")

        makeblastdb = os.path.join(os.path.dirname(sys.executable), "makeblastdb")
        subprocess.call(["makeblastdb", "-in", seqs_file, "-dbtype", "prot","-title", "HistoneDB"])
