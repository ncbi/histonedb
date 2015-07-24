import os
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from tools.L_shade_hist_aln import write_alignments
from tools.hist_ss import get_gff_from_align
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from browse.models import Sequence

class Command(BaseCommand):
    help = 'Reset sequence features'
    seed_directory = os.path.join(settings.STATIC_ROOT, "browse", "seeds")

    def add_arguments(self, parser):
         parser.add_argument(
            "-f", 
            "--force", 
            default=False, 
            action="store_true", 
            help="Force the creation of PDFs, GFFs even if the files exist")
        
    def handle(self, *args, **options):
        for variant, seed in self.get_seeds():
            if not os.path.exists("{}.pdf".format(seed[:-6])) or options["force"]:
                #Write PDF
                write_alignments([seed], seed[:-6], save_dir=os.path.dirname(seed))

            if not os.path.exists("{}.gff".format(seed[:-6])) or options["force"]:
                #Write GFF
                with open("{}.gff".format(seed[:-6]), "w") as gff:
                    msa = MultipleSeqAlignment(list(SeqIO.parse(seed, "fasta")))
                    get_gff_from_align(msa, gff, save_dir=os.path.dirname(seed))

            #Set reviewed to True
            not_found = {}
            for num_seq, s in enumerate(SeqIO.parse(seed, "fasta")):
                fields = s.id.split("|")
                try:
                    #Core seed gi is first index, but we want to ignore them
                    id = int(fields[0])
                    continue
                except ValueError:
                    try:
                        #Variant seed is second index
                        id = int(fields[1])
                    except ValueError:
                        try:
                            not_found[seed[:-6]].append(s.id)
                        except KeyError:
                            not_found[seed[:-6]] = [s.id]
                try:
                    s = Sequence.objects.get(id=str(id), variant__id=variant)
                    s.reviewed = True
                    s.save()
                except Sequence.DoesNotExist:
                    try:
                        not_found[variant].append(id)
                    except KeyError:
                        not_found[variant] = [id]
            print not_found

    def get_seeds(self):
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            for seed in files: 
                if not seed.endswith(".fasta"): continue
                variant = os.path.basename(seed)[:-6]
                if i == 0:
                    variant = "canonical{}".format(variant)
                yield variant, os.path.join(root, seed)