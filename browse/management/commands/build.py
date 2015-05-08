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
        #Download nr
        nr_file = "nr"
        if not os.isfile(nr_file):
            process = Popen(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz", ">", "nr.gz"])
            process.communicate()
            process = Popen(["tar", "-xf", "nr.gz"])
            os.remove("nr.gz")

        #Build and search HMMs
        seed_directory = os.path.join("static", "browse", "seeds")
        hmm_directory = os.path.join("static", "browse", "hmms", options["core_type"])
        combined_varaints_file = os.path.join(hmm_directory, "combined_variants.hmm")
        with open(combined_varaints_file) as combined_variants:
            for root, _, files in os.walk(seed_directory):
                core_type = os.path.basename(root)
                if not core_type:
                    #Skip parent directory, only allow variant hmms to be built/searched
                    continue
                for seed in files:
                    if not seed.endswith(".fasta"): continue
                    hmm_file = os.path.join(hmm_directory, "{}.hmm".format(seed.split("_")[0]))
                    process = Popen(["hmmbuild", hmm, os.path.join("seeds", seed)])
                    process.communicate()
                    with open(hmm_file) as hmm:
                        print >> combined_variants, hmm.read().rstrip()

        #Press the HMMs
        process = Popen(["hmmpress", combined_varaints_file])
        process.communicate()

        #Use HMMs to search the nr database
        results_file = os.path(hmm_directory, "{}.out".format(options["core_type"]))
        process = Popen(["hmmsearch", "-o", results_file, "--cpu", "4", "--notextw", "--incdomE", "0.1", combined_varaints_file, nr_file])
        process.communicate()

        #Load data into the histone database
        load_hmmsearch(results_file, nr_file)