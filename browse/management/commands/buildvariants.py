from django.core.management.base import BaseCommand, CommandError
from browse.models import *
from tools.load_hmmsearch import parseHmmer
from tools.test_model import test_model
import subprocess
import os, sys

from Bio import SeqIO

class Command(BaseCommand):
    help = 'Build the HistoneDB by training HMMs with seed sequences found in seeds directory in the top directory of this project and using those HMMs to search the NR database. To add or update a single variant, you must rebuild everything using the --force option'
    seed_directory = os.path.join("static", "browse", "seeds")
    hmm_directory = os.path.join("static", "browse", "hmms")
    combined_varaints_file = os.path.join(hmm_directory, "combined_variants.hmm")
    combined_core_file = os.path.join(hmm_directory, "combined_cores.hmm")
    pressed_combined_varaints_file = os.path.join(hmm_directory, "combined_variants.h3f")
    pressed_combined_cores_file = os.path.join(hmm_directory, "combined_cores.h3f")
    results_file = "{}.out".format(combined_varaints_file)
    core_results_file = "{}.out".format(combined_core_file)
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

        parser.add_argument(
            "--only_variants",
            default=False,
            action="store_true",
            help="Build and evaluate only variant HMMs. Default False, evalute both core and variant HMMs")

        parser.add_argument(
            "--only_cores",
            default=False,
            action="store_true",
            help="Build and evaluate only core HMMs. Default False, evalute both core and variant HMMs")



    def handle(self, *args, **options):
        self.get_nr(force=options["nr"])

        if not options["force"] and os.path.isfile(self.results_file):
            self.load()

        elif not options["force"] and 
          ((not options["only_cores"] and os.path.isfile(self.pressed_combined_varaints_file)) or
            (not options["only_varaints"] and os.path.isfile(self.pressed_combined_combined_file)):
            self.test(only_cores=options["only_cores"], only_variants=options["only_variants"])
            if not options["only_cores"]: 
                self.search_varaint()
                self.load_variants()
            if not options["only_variants"]: 
                self.search_core()
                self.load_cores()

        else:
            #Force to rebuild everything
            self.build(only_cores=options["only_cores"], only_variants=options["only_variants"])
            self.test(only_cores=options["only_cores"], only_variants=options["only_variants"])
            if not options["only_cores"]:
                self.press_variants()
                self.search_varaints()
                self.load_variants()
            if not options["only_variants"]: 
                self.press_cores()
                self.search_cores()
                self.load_cores()


    def get_nr(self, force=False):
        """Download nr if not present"""
        if not os.path.isfile(self.nr_file) or force:
            print >> self.stdout, "Downloading nr..."
            subprocess.call(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz", ">", "nr.gz"])
            subprocess.call(["tar", "-xf", "nr.gz"])
            os.remove("nr.gz")

    def build(self, only_cores=False, only_variants=False):
        """Build HMMs"""
        print >> self.stdout, "Building HMMs..."
        
        with open(self.combined_varaints_file, "w"), open(self.combined_core_file, "w") as combined_variants, combine_core:
            for core_type, seed in self.get_seeds(core=True):
                if seed is None and not only_variants:
                    #Build Core HMM
                    core_hmm_file = os.path.join(self.hmm_directory, "{}.hmm".format(core_type))
                    self.build_hmm(core_type, core_hmm_file, "{}.fasta".format(core_hmm_file[:-6]))
                    with open(core_hmm_file) as core_hmm:
                        print >> combine_core, core_hmm.read().rstrip()
                    continue

                if only_cores:
                    continue

                #Build Variant HMM
                hmm_dir = os.path.join(self.hmm_directory, core_type)
                if not os.path.exists(hmm_dir):
                    os.makedirs(hmm_dir)
                hmm_file = os.path.join(hmm_dir, "{}.hmm".format(seed[:-6]))
                self.build_hmm(seed[:-6], hmm_file, os.path.join(self.seed_directory, core_type, seed))
                with open(hmm_file) as hmm:
                    print >> combined_variants, hmm.read().rstrip()

    def build_hmm(self, name, db, seqs):
        subprocess.call(["hmmbuild", "-n", name, db, seqs])

    def press_variants(self):
        return self.press(self.combined_varaints_file)

    def press_cores(self):
        return self.press(self.combined_core_file)

    def press(self, combined_hmm):
        """Press the HMMs into a single HMM file, overwriting if present"""
        print >> self.stdout, "Pressing HMMs..."
        subprocess.call(["hmmpress", "-f", combined_hmm])

    def search_variants(self):
        return self.search(db=self.combined_varaints_file, out=self.results_file)

    def search_core(self):
        return self.search(db=self.combined_core_file, out=self.core_results_file)

    def search(self, db, out, sequences=None, E=10):
        """Use HMMs to search the nr database"""
        print >> self.stdout, "Searching HMMs..."

        if sequences is None:
            sequences = self.nr_file

        subprocess.call(["hmmsearch", "-o", out, "-E", str(E), "--cpu", "4", "--notextw", db, sequences])

    def load(self):
        """Load data into the histone database"""
        print >> self.stdout, "Loading data into HistoneDB..."
        load_variants(self.results_file, self.nr_file)
        load_cores(self.core_results_file)

    def get_seeds(self, core=False):
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            core_type = os.path.basename(root)
            if i==0:
                if core:
                    yield core_type, None
                #Skip parent directory, only allow variant hmms to be built/searched
                continue
            for seed in files: 
                if not seed.endswith(".fasta"): continue
                yield core_type, seed

    def test(self, only_cores=False, only_variants=False, specificity=0.95):
        for core_type, seed1 in self.get_seeds(core=True):
            #Seed1 is positive examples
            if seed1 == None:
                #Test Core HMM
                variant = "canonical{}".format(core_type)
                if only_variants: continue
            else:
                #Test Varaint HMM
                variant = seed1[:-6]
                if only_cores: continue
            
            hmm_file = os.path.join(self.hmm_directory, core_type, "{}.hmm".format(variant))

            output_dir = os.path.join(self.model_evalution, core_type)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            positive_examples_file = os.path.join(output_dir, "{}_postive_examples.fasta".format(variant))
            positive_examples = os.path.join(output_dir, "{}_postive_examples.out".format(variant))
            negative_examples_file = os.path.join(output_dir, "{}_negative_examples.fasta".format(variant))
            negative_examples = os.path.join(output_dir, "{}_negative_examples.out".format(variant))

            with open(positive_examples_file, "w") as positive_file:
                for s in SeqIO.parse(os.path.join(self.seed_directory, core_type, seed1), "fasta"):
                    s.seq = s.seq.ungap("-")
                    SeqIO.write(s, positive_file, "fasta")
            
            self.search(db=hmm_file, out=positive_examples, sequences=positive_examples_file, E=500)
            
            #Build negative examples from all other varaints
            
            with open(negative_examples_file, "w") as negative_file:
                for core_type2, seed2 in self.get_seeds(core=True):
                    if seed1 is None and seed2 is None:
                        #Compare Cores
                        if core_type == core_type2:
                            #Don't compare cores if they are the same
                            continue
                        os.path.join(self.seed_directory, "{}.fasta".format(core_type2))
                    elif (seed1 is not None and seed2 is None) or (seed1 is None and seed2 is not None):
                        #Do not compare core to varaints
                        continue
                    if not seed1 == seed2:
                        sequences = os.path.join(self.seed_directory, core_type2, seed2)
                    else:
                        #Do not compare files if they are the same
                        continue

                    for s in SeqIO.parse(sequences, "fasta"):
                        s.seq = s.seq.ungap("-")
                        SeqIO.write(s, negative_file, "fasta")

            self.search(db=hmm_file, out=negative_examples, sequences=negative_examples_file, E=500)

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


