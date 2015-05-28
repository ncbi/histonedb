from django.core.management.base import BaseCommand, CommandError
from browse.models import *
from tools.load_hmmsearch import parseHmmer
from tools.test_model import test_model
import subprocess
import os, sys

from Bio import SeqIO

class Command(BaseCommand):
    help = 'Build the HistoneDB by training HMMs with seed sequences found in seeds directory in the top directory of thi project and using those HMMs to search the NR database.'
    seed_directory = os.path.join("static", "browse", "seeds")
    hmm_directory = os.path.join("static", "browse", "hmms")
    combined_varaints_file = os.path.join(hmm_directory, "combined_variants.hmm")
    pressed_combined_varaints_file = os.path.join(hmm_directory, "combined_variants.h3f")
    results_file = "{}.out".format(combined_varaints_file)
    nr_file = "nr"
    model_evalution = os.path.join(hmm_directory, "model_evalution")

    def add_arguments(self, parser):
        parser.add_argument(
            "-f", 
            "--force", 
            default=False, 
            action="store_true", 
            help="Force the recreation of the HistoneDB. If True and an hmmsearch file is not present, the program will redownload the nr db if not present, build and press the HMM, and search the nr databse using the HMMs.")

        parser.add_argument(
            "--nr",
            default=False,
            action="store_true",
            help="Redownload NR database. If force is also True, this option is redundant")
    def handle(self, *args, **options):
        self.get_nr(force=options["nr"])

        if not options["force"] and os.path.isfile(self.results_file):
            self.load()

        elif not options["force"] and os.path.isfile(self.pressed_combined_varaints_file):
            self.test()
            self.search()
            self.load()

        else:
            #Force to rebuild everything
            self.build()
            self.test()
            self.press()
            self.search()
            self.load()


    def get_nr(self, force=False):
        """Download nr if not present"""
        if not os.path.isfile(self.nr_file) or force:
            print >> self.stdout, "Downloading nr..."
            subprocess.call(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz", ">", "nr.gz"])
            subprocess.call(["tar", "-xf", "nr.gz"])
            os.remove("nr.gz")

    def build(self):
        """Build HMMs"""
        print >> self.stdout, "Building HMMs..."
        
        with open(self.combined_varaints_file, "w") as combined_variants:
            for core_type, seed in self.get_seeds():
                hmm_dir = os.path.join(self.hmm_directory, core_type)
                if not os.path.exists(hmm_dir):
                    os.makedirs(hmm_dir)
                hmm_file = os.path.join(hmm_dir, "{}.hmm".format(seed[:-6]))
                subprocess.call(["hmmbuild", "-n", seed[:-6], hmm_file, os.path.join(self.seed_directory, core_type, seed)])
                with open(hmm_file) as hmm:
                    print >> combined_variants, hmm.read().rstrip()

    def press(self):
        """Press the HMMs into a simgle HMM file"""
        print >> self.stdout, "Pressing HMMs..."
        subprocess.call(["hmmpress", self.combined_varaints_file])

    def search(self, db=None, sequences=None, out=None, E=10):
        """Use HMMs to search the nr database"""
        print >> self.stdout, "Searching HMMs..."

        if db is None:
            db = self.combined_varaints_file

        if sequences is None:
            sequences = self.nr_file

        if out is None:
            out = self.results_file

        subprocess.call(["hmmsearch", "-o", out, "-E", str(E), "--cpu", "4", "--notextw", db, sequences])

    def load(self):
        """Load data into the histone database"""
        print >> self.stdout, "Loading data into HistoneDB..."
        parseHmmer(self.results_file, self.nr_file)

    def get_seeds(self):
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            core_type = os.path.basename(root)
            if i==0:
                #Skip parent directory, only allow variant hmms to be built/searched
                continue
            for seed in files: 
                if not seed.endswith(".fasta"): continue
                yield core_type, seed

    def test(self, specificity=0.95):
        for core_type, seed1 in self.get_seeds():
            #Seed1 is positive examples
            variant = seed1[:-6]
            hmm_file = os.path.join(self.hmm_directory, core_type, "{}.hmm".format(variant))

            output_dir = os.path.join(self.model_evalution, core_type)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            positive_examples_file = os.path.join(output_dir, "{}_postive_examples.fasta".format(variant))
            positive_examples = os.path.join(output_dir, "{}_postive_examples.out".format(variant))
            with open(positive_examples_file, "w") as positive_file:
                for s in SeqIO.parse(os.path.join(self.seed_directory, core_type, seed1), "fasta"):
                    s.seq = s.seq.ungap("-")
                    SeqIO.write(s, positive_file, "fasta")
            
            self.search(db=hmm_file, sequences=positive_examples_file, out=positive_examples, E=500)
            
            #Build negative examples from all other varaints
            negative_examples_file = os.path.join(output_dir, "{}_negative_examples.fasta".format(variant))
            negative_examples = os.path.join(output_dir, "{}_negative_examples.out".format(variant))
            with open(negative_examples_file, "w") as negative_file:
                for core_type2, seed2 in self.get_seeds():
                    if seed1 == seed2: continue
                    for s in SeqIO.parse(os.path.join(self.seed_directory, core_type2, seed2), "fasta"):
                        s.seq = s.seq.ungap("-")
                        SeqIO.write(s, negative_file, "fasta")

            self.search(db=hmm_file, sequences=negative_examples_file, out=negative_examples, E=500)

            parameters = test_model(variant, output_dir, positive_examples, negative_examples)

            try:
                variant_model = Variant.objects.get(id=variant)
            except:
                if "H2A" in variant:
                    core_histone = Histone.objects.get(id="H2A")
                elif "H2B" in variant:
                    core_histone = Histone.objects.get(id="H2B")
                elif "H3" in variant:
                    core_histone = Histone.objects.get(id="H3")
                elif "H4" in variant:
                    core_histone = Histone.objects.get(id="H4")
                elif "H1" in variant:
                    core_histone = Histone.objects.get(id="H1")
                else:
                    continue
                variant_model = Variant(id=variant, core_type=core_histone)
                variant_model.save()
            variant_model.hmmthreshold = parameters["threshold"]
            variant_model.aucroc = parameters["roc_auc"]
            variant_model.save()


