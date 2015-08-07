from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from browse.models import Histone, Variant, Sequence, Score, Features 
from tools.load_hmmsearch import load_hmm_results, add_score, easytaxonomy_from_header
from tools.test_model import test_model
import subprocess
import os, sys
import re
from tools.taxonomy_from_gis import taxids_from_gis

from Bio import SearchIO
from Bio import SeqIO

#This command is the main one in creating the histone database system from seed alignments
#and by using HMMs constructed based on these alignment to classify the bigger database.
#see handle() for the workflow description.
#INPUT NEEDED: static/browse/seeds
#db_file (nr) in the main directory - the database of raw sequences.

class Command(BaseCommand):
    help = 'Build HistoneDB by first loading the seed sequences and then parsing the database file'
    seed_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "seeds")
    hmm_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "hmms")
    combined_hmm_file = os.path.join(hmm_directory, "combined_hmm.hmm")
    pressed_combined_hmm_file = os.path.join(hmm_directory, "combined_hmm.h3f")
    dbauto_search_results_file = os.path.join(hmm_directory, "dbauto_search_results.out")
    curated_all_fasta=os.path.join(hmm_directory, "curated_all.fasta")
    curated_search_results_file = os.path.join(hmm_directory, "curated_all_search_results.out")
    model_evalution = os.path.join(hmm_directory, "model_evalution")

    def add_arguments(self, parser):
        parser.add_argument(
            "-f", 
            "--force", 
            default=False,
            action="store_true", 
            help="Force the regeneration of HMM from seeds, HUMMER search in db_file, Test models and loading of results to database")

        parser.add_argument(
            "--db",
            nargs=1,
            dest="dbauto_file",
            default="nr",
            help="Specify the database file, by default will use or download nr")


    def handle(self, *args, **options):
        ##If no nr file is present in the main dir, will download nr from the NCBI ftp.
        self.dbauto_file=options['dbauto_file']
        if(self.dbauto_file=="nr"):
            self.get_nr()
        ##If force=True, we do the default procedure:
        if options["force"]:
            ####Clean the DB first of all
            Sequence.objects.all().delete()
            Features.objects.all().delete()
            Score.objects.all().delete()
            Variant.objects.all().delete()
            Histone.objects.all().delete()
            ####Let's start by populating our Histone types table, add some descriptions to them.
            self.create_histone_types()
            ####We start next by creating HMMs from seeds and compressing them to one HMM file.
            self.build_hmms_from_seeds()
            ####HMMpress the combined list of HMMs to use it in fast search
            self.press_hmms()
            ####We need to determing HUMMER thresholds params,
            ####that we would in HMM search resutls analysis to decide if seq fits to be a certain variant
            ####AND create variant recoreds in database and load variant params into database.
            self.estimate_thresholds()
            ####Now we need to load our curated sets (which is indeed taken from seed alignments)
            ####into the database
            ####but we will subject them to the same classification as automatic seqs.
            self.load_curated()
            self.get_scores_for_curated_via_hmm()
            ####Now we will do the atomated seach in the dbauto_file (nr or other) database
            #### and load all automatically classified sequnces, as well as add scores to them.
            self.search_in_dbauto()
            self.load_from_dbauto()
            exit()
            ####Now we need to annotate sequneces with features: secondary structure, strucutral features, functional features.
            exit()

        #esle we only try to do things that are lacking
        #do every step only if the corresponding files are missing
        else:
            #TODO: Need to modify this section!!!
            print "Not currently implemented!!! Use -f"
            pass
############An ad-hoc snippet to load the seeds sequence into the database.###########
        try:
            os.remove('/tmp/all.fasta')
        except:
            pass
        for core_type, seed1 in self.get_seeds(core=True):
            print core_type, seed1
            #Seed1 is positive examples
            if seed1 == None:
                #Test Core HMM
                variant = "canonical{}".format(core_type)
                positive_path = os.path.join(self.seed_directory, "{}.fasta".format(core_type))
            else:
                #Test Varaint HMM
                variant = seed1[:-6]
                positive_path = os.path.join(self.seed_directory, core_type, seed1)
            print "Saving", positive_path, "into", "/tmp/all.fasta"
            with open("/tmp/all.fasta", "a") as positive_file:
                for s in SeqIO.parse(positive_path, "fasta"):
                    s.seq = s.seq.ungap("-")
                    SeqIO.write(s, positive_file, "fasta")
#Search variants
        self.search(self.combined_varaints_file, out=self.results_file,sequences="/tmp/all.fasta")
#Search cores
        self.search(db=self.combined_core_file, out=self.core_results_file,sequences="/tmp/all.fasta")
        self.load_variants(reset=False)
        self.load_cores(reset=False)

    #####################

    def create_histone_types(self):
        """Create basic histone types
        Check that these names are consistent with the folder names in static/browse/seeds/
        """
        for i in ['H3','H4','H2A','H2B']:
            obj,created = Histone.objects.get_or_create(id=i,taxonomic_span="Eukaryotes",\
                      description="Core histone")
        if created:
            print "Histone ", obj," type was created in database."

        obj,created = Histone.objects.get_or_create(id="H1",taxonomic_span="Eukaryotes",\
                      description="Linker histone")
        if created:
            print "Histone ", obj," type was created in database."

    def get_nr(self):
        """Download nr if not present"""
        if not os.path.isfile(self.dbauto_file):
            print >> self.stdout, "Downloading nr..."
            with open("nr.gz", "w") as nrgz:
                subprocess.call(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"], stdout=nrgz)
            subprocess.call(["gunzip", "nr.gz"])

    def build_hmms_from_seeds(self):
        """Build HMMs from seed histone sequences,
        and outputing them to
        static/browse/hmms/
        to individual dirs as well as combining to pne file combined_hmm_file
        """
        print >> self.stdout, "Building HMMs..."
        
        with open(self.combined_hmm_file, "w") as combined_hmm:
            for hist_type, seed in self.get_seeds():
                #Build HMMs
                hmm_dir = os.path.join(self.hmm_directory, hist_type)
                if not os.path.exists(hmm_dir):
                    os.makedirs(hmm_dir)
                hmm_file = os.path.join(hmm_dir, "{}.hmm".format(seed[:-6]))
                self.build_hmm(seed[:-6], hmm_file, os.path.join(self.seed_directory, hist_type, seed))
                with open(hmm_file) as hmm:
                    print >> combined_hmm, hmm.read().rstrip()

    def build_hmm(self, name, db, seqs):
        print ["hmmbuild", "-n", name, db, seqs]
        subprocess.call(["hmmbuild", "-n", name, db, seqs])

    def press_hmms(self):
        return self.press(self.combined_hmm_file)

    def press(self, combined_hmm):
        """Press the HMMs into a single HMM file, overwriting if present"""
        print >> self.stdout, "Pressing HMMs..."
        subprocess.call(["hmmpress", "-f", combined_hmm])

    def search_in_dbauto(self):
        return self.search(hmms_db=self.combined_hmm_file, out=self.dbauto_search_results_file,sequences=self.dbauto_file)

    def search(self, hmms_db, out, sequences=None, E=10):
        """Use HMMs to search the nr database"""
        print >> self.stdout, "Searching HMMs..."

        if sequences is None:
            sequences = self.dbauto_file

        subprocess.call(["hmmsearch", "-o", out, "-E", str(E), "--cpu", "4", "--notextw", hmms_db, sequences])


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
            print "Loading hmmsearch for variant:", variant_query.id
            variant_model=Variant.objects.get(id=variant_query.id)
            for hit in variant_query:
                try:
                    gi=re.search('(\d+)\|',hit.id).group(1)
                except:
                    try:
                        gi=re.search('(NOGI_\d+)\|',hit.id).group(1)
                    except:
                        print "NOT Loading ", hit.id
                        continue
                seq = Sequence.objects.get(id=gi)
                for i, hsp in enumerate(hit):
                    add_score(seq, variant_model, hsp)



    def load_from_dbauto(self,reset=True):
        """Load data into the histone database"""
        print >> self.stdout, "Loading data into HistoneDB..."
        load_hmm_results(self.dbauto_search_results_file, reset=reset)


    def load_curated(self):
        """
        Extracts sequences from seed alignments in static/browse/seeds
        Loads them into the database with flag reviewed=True (which means curated)
        An important fact:
        the seqs in seeds, should have a special header currently:
        >Ixodes|241122402|macroH2A Ixodes_macroH2A
        we accept only this patterns to extract GIs
        """
        no_gi_index=1
        for hist_type, seed in self.get_seeds():
            variant_name = seed[:-6]
            print variant_name,"==========="
            seed_aln_file = os.path.join(self.seed_directory, hist_type, seed)
            gis=[]
            for s in SeqIO.parse(seed_aln_file, "fasta"):
                # print s.id, "!!!!!!"
                # print s.description
                s.seq = s.seq.ungap("-")
                try:
                    gi=re.search('(\d+)\|',s.id).group(1)
                    temp_tax=1
                    gis.append(gi)
                except:
                    gi=re.search('(NOGI_\d+)\|',s.id).group(1)
                    print "NO GI detected ", s.id
                    print "Assigning NOGI_",no_gi_index
                    temp_tax=easytaxonomy_from_header(s.id).id
                    # continue # we will now try to put them in database
                print "Loading ", s.id
                #here is a trick, to assign dummy taxid at first and then batch update
                seq = Sequence(
                id       = gi,
                variant_id  = variant_name,
                gene     = None,
                splice   = None,
                taxonomy_id = temp_tax,
                #header   = s.id+" "+s.description, #desc includes id
                header   = "CURATED SEQUENCE: "+s.description,
                sequence = s.seq,
                reviewed = True,
                )
                seq.save()
            #Now let's lookup taxid for those having GIs via NCBI.
            for taxid,gi in zip(taxids_from_gis(gis),gis):
                seq=Sequence.objects.get(pk=gi)
                seq.taxonomy_id=taxid
                seq.save()






    def get_seeds(self):
        """
        Goes through static/browse/seeds directories and returns histone type names and fasta file name of variant (without path).
        """
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            hist_type = os.path.basename(root)
            for seed in files: 
                if not seed.endswith(".fasta"): continue
                yield hist_type, seed

    def estimate_thresholds(self, specificity=0.95):
        """
        Estimate HMM threshold that we will use for variant classification.
        Construct two sets for every variant - negative and positive.
        And estimate params from ROC-curves.
        """
        for hist_type_pos, seed_pos in self.get_seeds():
            variant_name = seed_pos[:-6]

            #Getting all paths right
            positive_seed_aln_file = os.path.join(self.seed_directory, hist_type_pos, seed_pos)
            hmm_file = os.path.join(self.hmm_directory, hist_type_pos, "{}.hmm".format(variant_name))

            output_dir = os.path.join(self.model_evalution, hist_type_pos)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            positive_examples_file = os.path.join(output_dir, "{}_postive_examples.fasta".format(variant_name))
            positive_examples_out = os.path.join(output_dir, "{}_postive_examples.out".format(variant_name))
            negative_examples_file = os.path.join(output_dir, "{}_negative_examples.fasta".format(variant_name))
            negative_examples_out = os.path.join(output_dir, "{}_negative_examples.out".format(variant_name))

            #Unagapping all sequence from seed aln - this will be the positive example
            # print "Saving", positive_seed_aln_file, "into", positive_examples_file
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
                    # print "Saving negatives from", sequences, "into", negative_examples_file
                    for s in SeqIO.parse(sequences, "fasta"):
                        s.seq = s.seq.ungap("-")
                        SeqIO.write(s, nf, "fasta")
            #Searching through negative example set
            self.search(hmms_db=hmm_file, out=negative_examples_out, sequences=negative_examples_file, E=500)

            #Here we are doing ROC curve analysis and returning parameters
            #TODO still needs to be refurbished, breaks on H1.8
            try:
                parameters = test_model(variant_name, output_dir, positive_examples_out, negative_examples_out)
            except:
                print "Test model, failed for ", variant_name
                parameters["threshold"] = 100.

            # #Some addhoc, because the test_model breaks on it, redid it as seen above
            # if variant_name == "H1.8":
            #     parameters["threshold"] = 100.

            #Let's put the parameter data to the database,
            #We can set hist_type directly by ID, which is hist_type_pos in this case - because it is the primary key in Histone class.
            variant_model, create = Variant.objects.get_or_create(id=variant_name,hist_type_id=hist_type_pos)
            if create:
                print "Created ",variant_model," variant model in database"
            print "Updating thresholds for ", variant_model
            variant_model.hmmthreshold = parameters["threshold"]
            variant_model.aucroc = parameters["roc_auc"]
            variant_model.save()


