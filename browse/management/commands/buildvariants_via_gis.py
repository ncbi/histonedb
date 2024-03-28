from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from browse.models import Histone, Variant, Sequence, Score, Feature
from tools.load_hmmsearch_with_gis import load_hmm_results, add_score, get_many_prot_seqrec_by_gis
from tools.test_model import test_model
import subprocess
import os, sys
import re
import io
from tools.taxonomy_from_gis import taxonomy_from_header, easytaxonomy_from_header, taxids_from_gis, update_taxonomy_for_gis

from Bio import SearchIO
from Bio import SeqIO

#This command is the main one in creating the histone database system from seed alignments
#and by using HMMs constructed based on these alignment to classify the bigger database.
#see handle() for the workflow description.
#INPUT NEEDED: static/browse/seeds_gis
#db_file (nr) in the main directory - the database of raw sequences.

class Command(BaseCommand):
    help = 'Build HistoneDB by first loading the seed sequences and then parsing the database file'
    seed_directory = os.path.join(settings.STATIC_ROOT, "browse", "seeds_gis")
    hmm_directory = os.path.join(settings.STATIC_ROOT, "browse", "hmms")
    combined_hmm_file = os.path.join(hmm_directory, "combined_hmm.hmm")
    pressed_combined_hmm_file = os.path.join(hmm_directory, "combined_hmm.h3f")
    db_search_results_file = os.path.join(hmm_directory, "db_search_results.out")
    curated_all_fasta=os.path.join(hmm_directory, "curated_all.fasta")
    curated_search_results_file = os.path.join(hmm_directory, "curated_all_search_results.out")
    model_evaluation = os.path.join(hmm_directory, "model_evaluation")
    ids_file = os.path.join(settings.STATIC_ROOT, "browse", "blast", "HistoneDB_sequences.ids")
    full_length_seqs_file = os.path.join(settings.STATIC_ROOT, "browse", "blast", "HistoneDB_sequences.fasta")

    def add_arguments(self, parser):
        parser.add_argument(
            "-f", 
            "--force", 
            default=False,
            action="store_true", 
            help="Force the regeneration of HMM from seeds, HUMMER search in db_file, Test models and loading of results to database")

        parser.add_argument(
            "--db",
            dest="db_file",
            default="nr",
            help="Specify the database file, by default will use or download nr")

        parser.add_argument(
            "--adjust_thrshehold",
            default=False,
            action="store_true",
            help="Adjust threshold")


    def handle(self, *args, **options):
        ##If no nr file is present in the main dir, will download nr from the NCBI ftp.
        self.db_file=options['db_file']
        if self.db_file == "nr":
            self.get_nr()
        if self.db_file == "swissprot":
            self.get_swissprot()

        if options["force"]:
            #Clean the DB, removing all sequence/variants/etc
            Sequence.objects.all().delete()
            Score.objects.all().delete()
            Variant.objects.all().delete()
            Histone.objects.all().delete()

        #Populate our Histone types table add descriptions
        if options["force"]:
            self.create_histone_types()

        if options["force"] or not os.path.isfile(self.combined_hmm_file):
            #Create HMMs from seeds and compress them to one HMM file tp faster search with hmmpress.
            self.build_hmms_from_seeds()
            self.press_hmms()

            #Determine HMMER thresholds params used to decide if sequence is a variant ot not. If sequence is above 
            self.estimate_thresholds()

            #Load our curated sets taken from seed alignments into the database an run classification algorithm
            self.load_curated()
            self.get_scores_for_curated_via_hmm()

        if options["force"] or not os.path.isfile(self.db_search_results_file):
            #Search inputted seuqence database using our variantt models
            self.search_in_db()

        #Load the sequences and classify them based on thresholds
        self.load_from_db()
        self.extract_full_sequences_from_ncbi()

        # self.extract_full_sequences()
        self.canonical2H2AX()

    def canonical2H2AX(self):
        """Fix an issue where the canonical variant takes over sequence from H2A.X. 
        The H2A.X motif SQ[ED][YFL]$ is not strong enough, but is the correct variant.
        """
        for s in Sequence.objects.filter(variant="canonical_H2A",reviewed=False, sequence__regex="SQ[ED][YFLIA]$"):
            old_score = s.all_model_scores.get(used_for_classification=True)
            old_score.used_for_classification = False
            old_score.save()
            new_score, created = Score.objects.get_or_create(variant__id="H2A.X",sequence=s)
            new_score.used_for_classification = True
            new_score.regex = True
            s.variant_id="H2A.X"
            new_score.save()
            s.save()

    def create_histone_types(self):
        """Create basic histone types
        Check that these names are consistent with the folder names in static/browse/seeds_gis/
        """
        for i in ['H3','H4','H2A','H2B']:
            obj,created = Histone.objects.get_or_create(id=i,taxonomic_span="Eukaryotes",\
                      description="Core histone")
        if created:
            print("Histone ", obj," type was created in database.")

        obj,created = Histone.objects.get_or_create(id="H1",taxonomic_span="Eukaryotes",\
                      description="Linker histone")
        if created:
            print("Histone ", obj," type was created in database.")

    def get_nr(self):
        """Download nr if not present"""
        if not os.path.isfile(self.db_file):
            print("Downloading nr...", file=self.stdout)
            with open("nr.gz", "w") as nrgz:
                subprocess.call(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"], stdout=nrgz)
            subprocess.call(["gunzip", "nr.gz"])

    def get_swissprot(self):
        """Download nr if not present"""
        if not os.path.isfile(self.db_file):
            print("Downloading swissprot...", file=self.stdout)
            with open("swissprot.gz", "w") as swissprotgz:
                subprocess.call(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz"], stdout=swissprotgz)
            subprocess.call(["gunzip", "swissprot.gz"])


    def build_hmms_from_seeds(self):
        """Build HMMs from seed histone sequences,
        and outputing them to
        static/browse/hmms/
        to individual dirs as well as combining to pne file combined_hmm_file
        """
        print("Building HMMs...", file=self.stdout)
        
        with open(self.combined_hmm_file, "w") as combined_hmm:
            for hist_type, seed in self.get_seeds():
                #Build HMMs
                hmm_dir = os.path.join(self.hmm_directory, hist_type)
                if not os.path.exists(hmm_dir):
                    os.makedirs(hmm_dir)
                hmm_file = os.path.join(hmm_dir, "{}.hmm".format(seed[:-6]))
                self.build_hmm(seed[:-6], hmm_file, os.path.join(self.seed_directory, hist_type, seed))
                with open(hmm_file) as hmm:
                    print(hmm.read().rstrip(), file=combined_hmm)

    def build_hmm(self, name, db, seqs):
        print(["hmmbuild", "-n", name, db, seqs])
        subprocess.call(["hmmbuild", "-n", name, db, seqs])

    def press_hmms(self):
        return self.press(self.combined_hmm_file)

    def press(self, combined_hmm):
        """Press the HMMs into a single HMM file, overwriting if present"""
        print("Pressing HMMs...", file=self.stdout)
        subprocess.call(["hmmpress", "-f", combined_hmm])

    def search_in_db(self):
        return self.search(hmms_db=self.combined_hmm_file, out=self.db_search_results_file, sequences=self.db_file)

    def search(self, hmms_db, out, sequences=None, E=10):
        """Use HMMs to search the nr database"""
        print("Searching HMMs...", file=self.stdout)

        if sequences is None:
            sequences = self.db_file
        print(" ".join(["hmmsearch", "-o", out, "-E", str(E), "--cpu", "4", "--notextw", hmms_db, sequences]))
        subprocess.call(["hmmsearch", "-o", out, "-E", str(E), "--cpu", "4", "--notextw", hmms_db, sequences])

    def extract_full_sequences(self, sequences=None):
        """Create database to extract full length sequences"""

        if sequences is None:
            sequences = self.db_file

        #1) Create and index of sequence file
        print("Indexing sequence database {}...".format(sequences))
        subprocess.call(["esl-sfetch", "--index", sequences])

        #2) Extract all ids 
        print("Extracting full length sequences...")
        subprocess.call(["esl-sfetch", "-o", self.full_length_seqs_file, "-f", sequences, self.ids_file])

        #3) Update sequences with full length NR sequences -- is there a faster way?
        print("Updating records with full length sequences...")
        for record in SeqIO.parse(self.full_length_seqs_file, "fasta"):
            headers = record.description.split(" >")
            for header in headers:
                gi = header.split("|")[1]
                try:
                    seq = Sequence.objects.get(id=gi)
                    print("Updating sequence:", seq.description)
                    seq.sequence = str(record.seq)
                    seq.save()
                except Sequence.DoesNotExist:
                    pass

    def get_scores_for_curated_via_hmm(self):
        """
        For every curated variant we want to generate a set of scores against HMMs.
        This is needed to supply the same type of information for curated as well as for automatic seqs.
        """
        #Construct the one big file from all cureated seqs.
        with open(self.curated_all_fasta, "w") as f:
            for hist_type, seed in self.get_seeds():
                seed_aln_file = os.path.join(self.seed_directory, hist_type, seed)
                for s in SeqIO.parse(seed_aln_file, "fasta"):
                    s.seq = s.seq.ungap("-")
                    SeqIO.write(s, f, "fasta")
        #Search it by our HMMs
        self.search(hmms_db=self.combined_hmm_file, out=self.curated_search_results_file,sequences=self.curated_all_fasta)
        ##We need to parse this results file;
        ##we take here a snippet from load_hmmsearch.py, and tune it to work for our curated seq header format
        for variant_query in SearchIO.parse(self.curated_search_results_file, "hmmer3-text"):
            print("Loading hmmsearch for variant:", variant_query.id)
            variant_model=Variant.objects.get(id=variant_query.id)
            for hit in variant_query:
                gi = hit.id.split("|")[1]
                seq = Sequence.objects.get(id=gi)
                # print hit
                try: #sometimes we get this:    [No individual domains that satisfy reporting thresholds (although complete target did)]
                    best_hsp = max(hit, key=lambda hsp: hsp.bitscore)
                    add_score(seq, variant_model, best_hsp, seq.variant==variant_model)
                except:
                    pass

    def load_from_db(self,reset=True):
        """Load data into the histone database"""
        print("Loading data into HistoneDB...", file=self.stdout)
        load_hmm_results(self.db_search_results_file, self.ids_file)

    def load_curated(self):
        """
        Extracts sequences from seed alignments in static/browse/seeds_gis
        Loads them into the database with flag reviewed=True (which means curated)
        An important fact:
        the seqs in seeds, should have a special header currently:
        >Ixodes|241122402|macroH2A Ixodes_macroH2A
        we accept only this patterns to extract GIs
        """
        gis=[]
        for hist_type, seed in self.get_seeds():
            variant_name = seed[:-6]
            print(variant_name,"===========")
            seed_aln_file = os.path.join(self.seed_directory, hist_type, seed)
            for s in SeqIO.parse(seed_aln_file, "fasta"):
                s.seq = s.seq.ungap("-")
                gi = s.id.split("|")[1]
                if gi.startswith("NOGI"):
                    print("NO GI detected ", s.id)
                    taxid= easytaxonomy_from_header(s.id).id
                else:
                    #trick to make taxid retrieval faster
                    # taxonomy = taxonomy_from_header("", gi=gi)
                    taxid=1
                    gis.append(gi)
                print("Loading ", s.id)
                
                seq = Sequence(
                    id       = gi,
                    variant_id  = variant_name,
                    gene     = None,
                    splice   = None,
                    taxonomy_id = taxid,
                    header   = "CURATED SEQUENCE: {}".format(s.description),
                    sequence = s.seq,
                    reviewed = True,
               )
                seq.save()

        # Now let's lookup taxid for those having GIs via NCBI.
        update_taxonomy_for_gis(gis)

    def get_seeds(self):
        """
        Goes through static/browse/seeds_gis directories and returns histone type names and fasta file name of variant (without path).
        """
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            hist_type = os.path.basename(root)
            if hist_type=="seeds_gis": #means we are in top dir, we skip,
            # combinded alignmnents for hist types are their, but we do not use them in database constuction,
            #only in visualization on website
                continue
            for seed in files: 
                if not seed.endswith(".fasta"): continue
                yield hist_type, seed

    def estimate_thresholds(self, specificity=0.95):
        """
        Estimate HMM threshold that we will use for variant classification.
        Construct two sets for every variant:
            negative: The seed alignmnents from every other variant
            positive: the current seed alignment for the variant
        And estimate params from ROC-curves.
        """
        for hist_type_pos, seed_pos in self.get_seeds():
            variant_name = seed_pos[:-6]

            #Getting all paths right
            positive_seed_aln_file = os.path.join(self.seed_directory, hist_type_pos, seed_pos)
            hmm_file = os.path.join(self.hmm_directory, hist_type_pos, "{}.hmm".format(variant_name))

            output_dir = os.path.join(self.model_evaluation, hist_type_pos)

            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            positive_examples_file = os.path.join(output_dir, "{}_postive_examples.fasta".format(variant_name))
            positive_examples_out = os.path.join(output_dir, "{}_postive_examples.out".format(variant_name))
            negative_examples_file = os.path.join(output_dir, "{}_negative_examples.fasta".format(variant_name))
            negative_examples_out = os.path.join(output_dir, "{}_negative_examples.out".format(variant_name))

            #Unagapping all sequence from seed aln - this will be the positive example
            with open(positive_examples_file, "w") as pf:
                for s in SeqIO.parse(positive_seed_aln_file, "fasta"):
                    s.seq = s.seq.ungap("-")
                    SeqIO.write(s, pf, "fasta")

            #Searching the positive examples set
            self.search(hmms_db=hmm_file, out=positive_examples_out, sequences=positive_examples_file, E=500)
            
            #Build negative examples from all other varaints
            with open(negative_examples_file, "w") as nf:
                for hist_type_neg, seed_neg in self.get_seeds():
                    if((hist_type_pos == hist_type_neg) and (seed_neg == seed_pos)):
                        continue
                    else:
                        sequences = os.path.join(self.seed_directory, hist_type_neg, seed_neg)
                    
                    for s in SeqIO.parse(sequences, "fasta"):
                        s.seq = s.seq.ungap("-")
                        SeqIO.write(s, nf, "fasta")

            #Searching through negative example set
            self.search(hmms_db=hmm_file, out=negative_examples_out, sequences=negative_examples_file, E=500)

            #Here we are doing ROC curve analysis and returning parameters
            specificity = 0.8 if ("canonical" in variant_name) else 0.9 #Hack to make canoical have a lower threshold and ther variants higher threshold
            specificity = 0.05 if ("generic" in variant_name)  else specificity #Hack to make canoical have a lower threshold and ther variants higher threshold

            parameters = test_model(variant_name, output_dir, positive_examples_out, negative_examples_out, measure_threshold=specificity)

            #Let's put the parameter data to the database,
            #We can set hist_type directly by ID, which is hist_type_pos in this case - because it is the primary key in Histone class.
            variant_model, create = Variant.objects.get_or_create(id=variant_name,hist_type_id=hist_type_pos) #,hist_type_id=hist_type_pos)
            if create:
                print("Created ",variant_model," variant model in database")
            print("Updating thresholds for ", variant_model)
            variant_model.hmmthreshold = parameters["threshold"]
            variant_model.aucroc = parameters["roc_auc"]
            variant_model.save()

    def extract_full_sequences_from_ncbi(self):
        """Exract full seq by direct call to NCBI servers"""
        print("Getting full sequences of automatically annotated proteins from NCBI====")
        gis=Sequence.objects.filter(reviewed=False).values_list('id', flat=True)
        fasta_dict=get_many_prot_seqrec_by_gis(gis)

        #3) Update sequences with full length NR sequences -- is there a faster way?
        for gi,record in fasta_dict.items():
            headers = record.description.split(" >")
            for header in headers:
                gi = header.split("|")[1]
                # print gi
                try:
                    seq = Sequence.objects.get(id=gi)
                    seq.sequence = str(record.seq)
                    seq.save()
                except Sequence.DoesNotExist:
                    pass

