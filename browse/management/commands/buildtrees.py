from django.core.management.base import BaseCommand, CommandError
from browse.models import *
from tools.load_hmmsearch import parseHmmer
from tools.test_model import test_model
import subprocess
import os, sys
from Bio import Phylo
from Bio.Phylo import PhyloXMLIO

class Command(BaseCommand):
    help = 'Build the HistoneDB by training HMMs with seed sequences found in seeds directory in the top directory of thi project and using those HMMs to search the NR database.'
    seed_directory = os.path.join("static", "browse", "seeds")
    trees_path = os.path.join("static", "browse", "trees")
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
        self.make_trees()
        self.add_features()

    def get_variants(self, core_type=None):
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            if core_type and not os.path.basename(root) == core_type: continue
            if i==0:
                #Skip parent directory, only allow variant hmms to be built/searched
                continue
            for seed in files: 
                if not seed.endswith(".fasta"): continue
                if core_type is None:
                    yield core_type, seed
                else:
                    yield seed[:-6]

    def make_trees(self):
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            core_type = os.path.basename(root)
            if i==0:
                #Skip parent directory, only allow variant hmms to be built/searched
                continue
            combined_seed_file = os.path.join(self.trees_path, "{}.fasta".format(core_histone))
            combined_seed_aligned = os.path.join(self.trees_path, "{}_aligned.fasta".format(core_histone))
            with open(combined_seed_file, "w") as combined_seed:
                for seed in files: 
                    if not seed.endswith(".fasta"): continue
                    for s in SeqIO.parse(os.path.join(self.seed_directory, core_type, seed), "fasta"):
                        s.seq = s.seq.ungap("-")
                        SeqIO.write(s, combined_seed, "fasta")

            tree = os.path.join(self.trees_path, "{}.nwk".format(core_histone))
            subprocess.call(["muscle", "-in", combined_seed, '-out', combined_seed_aligned])
            subprocess.call(["clustalw2", "-infile", combined_seed_aligned, '-tree', '-out', tree])
            Phylo.convert(tree, 'newick',
                          os.path.join(self.trees_path, "{}.xml".format(core_histone)), 'phyloxml')
    
    def add_features(self):
        for core_histone in ["H2A", "H2B", "H3", "H4", "H1"]:
            phx = PhyloXMLIO.read(os.path.join(self.trees_path, "{}.xml".format(core_histone)), 'phyloxml')
            print phx
            print dir(phx)
            """for phylogeny in tree.phylogenies:
                for variant in get_variants(core_histone):
                    '<{0} fill="#000" stroke="#000" opacity="0.7" label="{0}" labelStyle="sectorHighlightText" />'.format(variant.id)
                for clade in phylogeny:
                    genus, gi, variant = clade.name.split("|")
                    clade.name = "{}_{}".format(genus, gi)
                    "<chart><group>{}</group></chart>".format(variant)
            PhyloXMLIO.write(phx, 'ex_no_other.xml')"""



