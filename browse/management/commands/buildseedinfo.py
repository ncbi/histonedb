import os
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from tools.L_shade_hist_aln import write_alignments
from tools.hist_ss import get_gff_from_align, get_variant_features
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from browse.models import Sequence
from Bio.Align.AlignInfo import SummaryInfo


class Command(BaseCommand):
    help = 'Reset sequence features'
    seed_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "seeds")

    def add_arguments(self, parser):
         parser.add_argument(
            "-f", 
            "--force", 
            default=False, 
            action="store_true", 
            help="Force the creation of PDFs, GFFs even if the files exist")
        
    def handle(self, *args, **options):
        save_dir = os.path.join(os.path.sep, "tmp", "HistoneDB")
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        for variant, seed in self.get_seeds():
            #PDF currently contaminates the dir with many files.
            #if not os.path.exists("{}.pdf".format(seed[:-6])) or options["force"]:
                #Write PDF
                #write_alignments([seed], seed[:-6], save_dir=os.path.dirname(seed))

            if not os.path.exists("{}.gff".format(seed[:-6])) or options["force"]:
                #Write GFF
                with open("{}.gff".format(seed[:-6]), "w") as gff:
                    msa = MultipleSeqAlignment(list(SeqIO.parse(seed, "fasta")))
                    a = SummaryInfo(msa)
                    cons = Sequence(id="Consensus", variant_id=variant, taxonomy_id=1, sequence=a.dumb_consensus(threshold=0.1, ambiguous='X').tostring())
                    features = get_variant_features(cons, save_dir=save_dir)
                    print >> gff, features
                    # get_gff_from_align(msa, gff, save_dir=os.path.dirname(seed))

            #Set reviewed to True
            not_found = {}
            for num_seq, s in enumerate(SeqIO.parse(seed, "fasta")):
                fields = s.id.split("|")
                try:
                    #Core seed gi is first index, but we want to ignore them
                    #This part is tricky, in seed alignmnets of variants fist argument is taxonomy name, second gi.
                    #but in cores (canonicals???) - gi is first.
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
