from django.core.management.base import BaseCommand, CommandError
from browse.models import *
from tools.load_hmmsearch import parseHmmer
from tools.test_model import test_model
import subprocess
import os, sys
from itertools import cycle
import StringIO

from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo import PhyloXML
from Bio.Phylo import PhyloXMLIO

import xml.etree.ElementTree as ET

import seaborn as sns
colors = cycle(sns.color_palette("Set2", 7))

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child] = clade
    return parents

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
            core_histone = os.path.basename(root)
            print "Creating tree for", core_histone
            if i==0:
                #Skip parent directory, only allow variant hmms to be built/searched
                continue

            #Combine all varaints for a core histone type into one unaligned fasta file
            combined_seed_file = os.path.join(self.trees_path, "{}.fasta".format(core_histone))
            combined_seed_aligned = os.path.join(self.trees_path, "{}_aligned.fasta".format(core_histone))
            with open(combined_seed_file, "w") as combined_seed:
                for seed in files: 
                    if not seed.endswith(".fasta"): continue
                    for s in SeqIO.parse(os.path.join(self.seed_directory, core_histone, seed), "fasta"):
                        s.seq = s.seq.ungap("-")
                        SeqIO.write(s, combined_seed, "fasta")

            tree = os.path.join(self.trees_path, "{}_aligned.ph".format(core_histone))
            subprocess.call(["muscle", "-in", combined_seed_file, '-out', combined_seed_aligned])
            subprocess.call(["clustalw2", "-infile={}".format(combined_seed_aligned), '-tree'])
            Phylo.convert(tree, 'newick',
                          os.path.join(self.trees_path, "{}_no_features.xml".format(core_histone)), 'phyloxml')
    
    def add_features(self):
        for core_histone in ["H2A", "H2B", "H3", "H4", "H1"]:
            tree_path = os.path.join(self.trees_path, "{}_no_features.xml".format(core_histone))
            phx = PhyloXMLIO.read(tree_path)

            tree = ET.parse(tree_path)

            render = PhyloXML.Other("render")
            
            parameters = PhyloXML.Other("parameters", namespace="")
            circular = PhyloXML.Other("circular", namespace="")
            circular.children.append(PhyloXML.Other("bufferRadius", value="0.5", namespace=""))
            parameters.children.append(circular)
            render.children.append(parameters)
            
            charts = PhyloXML.Other("charts", namespace="")
            charts.children.append(PhyloXML.Other("group", attributes={"type":"integratedBinary", "thickness":"20"}, namespace=""))
            render.children.append(parameters).append(charts)

            styles = PhyloXML.Other("styles", namespace="")
            for variant in self.get_variants(core_histone):
                color = colors.next()
                styles.children.append(PhyloXML.Other("{}".format(variant.replace("_","")), attributes={"fill":color, "stroke":color})
                styles.children.append(PhyloXML.Other("markdown{}".format(variant.replace("_","")), attributes={"fill":"#000", "stroke":"#000", "opacity":"0.7", "label":variant, "labelStyle":"sectorHighlightText"}, namespace=""))
            styles.children.append(PhyloXML.Other("sectorHighlightText", attributes={"font-family":"Verdana", "font-size":"14", "font-weight":"bold", "fill":"#FFFFFF", "rotate":"90"}, namespace=""))
            render.children.append(styles)

            phx.other.append(render)

            for phylogeny in phx.phylogenies:
                parents = all_parents(phylogeny)
                for clade in phylogeny.find_clades():
                    if clade.name:
                        try:
                            genus, gi, variant = clade.name.split("|")
                        except ValueError:
                            parent = parents[clade]
                            child = parent.clades.index(clade)
                            del parent.clades[child]
                            continue
                        clade.name = genus
                        chart = PhyloXML.Other("chart", namespace="")
                        chart.children.append(PhyloXML.Other("group", value=variant, namespace=""))
                        clade.other.append(chart)

                        annotation = PhyloXML.Other("annotation", namespace="")
                        annotation.children.append(PhyloXML.Other("desc", value="Variant {} from {} ({})".format(variant, genus, gi), namespace=""))
                        annotation.children.append(PhyloXML.Other("uri", value="variant/{}/".format(variant), namespace=""))
                        clade.other.append(annotation)

            with open(os.path.join(self.trees_path, "{}.xml".format(core_histone)), "w") as outfile:
                tree = StringIO.StringIO()
                PhyloXMLIO.write(phx, tree)
                tree = tree.getvalue().replace("ns0:", "")
                outfile.write(tree)

            



