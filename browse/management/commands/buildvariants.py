from django.core.management.base import BaseCommand, CommandError
from browse.models import *
from tools.load_hmmsearch import parseHmmer
import subprocess
import os, sys

class Command(BaseCommand):
    help = 'Build the HistoneDB by training HMMs with seed sequences found in seeds directory in the top directory of thi project and using those HMMs to search the NR database.'
    seed_directory = os.path.join("static", "browse", "seeds")
    hmm_directory = os.path.join("static", "browse", "hmms")
    combined_varaints_file = os.path.join(hmm_directory, "combined_variants.hmm")
    pressed_combined_varaints_file = os.path.join(hmm_directory, "combined_variants.h3f")
    results_file = "{}.out".format(combined_varaints_file)
    nr_file = "nr"

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
        self.get_nr(force=options["force"] or options["nr"])

        if not options["force"] and os.path.isfile(self.results_file):
            self.load()

        elif not options["force"] and os.path.isfile(self.pressed_combined_varaints_file):
            self.search()
            self.load()
        else:
            #Force to rebuild everything
            self.build()
            self.press()
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
            for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
                core_type = os.path.basename(root)
                print core_type
                if i==0:
                    #Skip parent directory, only allow variant hmms to be built/searched
                    continue
                for seed in files:
                    if not seed.endswith(".fasta"): continue
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

    def search(self):
        """Use HMMs to search the nr database"""
        print >> self.stdout, "Searching HMMs..."
        subprocess.call(["hmmsearch", "-o", self.results_file, "--cpu", "4", "--notextw", "--incdomE", "0.1", self.combined_varaints_file, self.nr_file])

    def load(self):
        """Load data into the histone database"""
        print >> self.stdout, "Loading data into HistoneDB..."
        parseHmmer(self.results_file, self.nr_file)
