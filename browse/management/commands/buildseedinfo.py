import os
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from tools.L_shade_hist_aln import write_alignments
from tools.hist_ss import get_gff_from_align
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

class Command(BaseCommand):
    help = 'Reset sequence features'
    seed_directory = os.path.join(settings.STATIC_ROOT, "browse", "seeds")

    def add_arguments(self, parser):
        pass
        
    def handle(self, *args, **options):
        for core_type, seed in self.get_seeds():
            write_alignments([seed], seed[:-6], save_dir=os.path.dirname(seed))

            with open("{}.gff".format(seed[:-6]), "w") as gff:
                msa = MultipleSeqAlignment(list(SeqIO.parse(seed, "fasta")))
                get_gff_from_align(msa, gff, save_dir=os.path.dirname(seed))

    def get_seeds(self):
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            for seed in files: 
                if not seed.endswith(".fasta"): continue
                core_type = seed.split(".")[0]
                yield core_type, os.path.join(root, seed)
