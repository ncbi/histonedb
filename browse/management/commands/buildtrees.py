from django.core.management.base import BaseCommand, CommandError
from browse.models import *

import os
from itertools import cycle
import StringIO
import subprocess

from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo import PhyloXML
from Bio.Phylo import PhyloXMLIO

import xml.etree.ElementTree as ET
ET.register_namespace("", "http://www.phyloxml.org/1.10/phyloxml.xsd")

colors = [
    "#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
    "#a6d854",
    "#ffd92f",
    "#e5c494"]

class Command(BaseCommand):
    help = 'Building data for variant trees using ClustalW2'
    seed_directory = os.path.join("static", "browse", "seeds")
    trees_path = os.path.join("static", "browse", "trees")
    def add_arguments(self, parser):
        parser.add_argument(
            "-f", 
            "--force", 
            default=False, 
            action="store_true", 
            help="Force the recreation of the varaint trees. If True and an phyloxml file is not present, the program will re-build each tree and add jsPhyloSVG features")

    def handle(self, *args, **options):
        self.make_trees(force=options["force"])
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

    def make_trees(self, force=False):
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            if i==0: #skip base path
                continue
            hist_type = os.path.basename(root)
            print "Creating tree for", hist_type

            final_tree_name = os.path.join(self.trees_path, "{}_no_features.xml".format(hist_type))
            if not force and os.path.isfile(final_tree_name):
                continue

            if not os.path.exists(self.trees_path):
                os.makedirs(self.trees_path)
            #Combine all varaints for a core histone type into one unaligned fasta file
            combined_seed_file = os.path.join(self.trees_path, "{}.fasta".format(hist_type))
            combined_seed_aligned = os.path.join(self.trees_path, "{}_aligned.fasta".format(hist_type))
            with open(combined_seed_file, "w") as combined_seed:
                for seed in files: 
                    if not seed.endswith(".fasta"): continue
                    for s in SeqIO.parse(os.path.join(self.seed_directory, hist_type, seed), "fasta"):
                        s.seq = s.seq.ungap("-")
                        SeqIO.write(s, combined_seed, "fasta")

            #Create trees and convert them to phyloxml
            tree = os.path.join(self.trees_path, "{}_aligned.ph".format(hist_type))
            subprocess.call(["muscle", "-in", combined_seed_file, '-out', combined_seed_aligned])
            subprocess.call(["clustalw2", "-infile={}".format(combined_seed_aligned), '-tree'])
            Phylo.convert(tree, 'newick', final_tree_name, 'phyloxml')
    
    def add_features(self):
        for hist_type in ["H2A", "H2B", "H3", "H1", "H4"]:
            print hist_type
            tree_path = os.path.join(self.trees_path, "{}_no_features.xml".format(hist_type))
            tree = ET.parse(tree_path)
            parent_map = {c: p for p in tree.getiterator() for c in p}

            for phylogeny in tree.iter("{http://www.phyloxml.org}phylogeny"):
                render = ET.Element("render")
                
                parameters = ET.Element("parameters")
                circular = ET.Element("circular")
                radius = ET.Element("bufferRadius")
                radius.text = "0.5"
                circular.append(radius)
                parameters.append(circular)
                render.append(parameters)
                
                charts = ET.Element("charts")
                group = ET.Element("group", attrib={"type":"integratedBinary", "thickness":"20"})
                charts.append(group)
                render.append(charts)

                styles = ET.Element("styles")
                variants = list(Variant.objects.filter(hist_type__id=hist_type).values_list("id", flat=True))
                for i, variant in enumerate(variants):
                    color = colors[i]
                    background = ET.Element("{}".format(variant.replace(".","")), attrib={"fill":color, "stroke":color})
                    if not variant.startswith(hist_type):
                        #Remove descriptor
                        start, name = variant[:variant.find(hist_type)], variant[variant.find(hist_type):]
                        if len(start) > 3 and start != "canonical":
                            start = start[0]
                        name = start+name
                    else:
                        name = variant
                    print "Adding", name
                    label = ET.Element("markgroup{}".format(variant.replace(".","")), attrib={"fill":"#000", "stroke":"#000", "opacity":"0.7", "label":name, "labelStyle":"sectorHighlightText"})
                    styles.append(background)
                    styles.append(label)
                label_sector = ET.Element("sectorHighlightText", attrib={"font-family":"Verdana", "font-size":"14", "font-weight":"bold", "fill":"#FFFFFF", "rotate":"90"})
                styles.append(label_sector)
                render.append(styles)

                phylogeny.insert(0, render)

                remove = []
                for clade in phylogeny.iter("{http://www.phyloxml.org}clade"):
                    name = clade.find("{http://www.phyloxml.org}name")
                    print "CLADE", clade
                    try:
                        print name.text
                    except:
                        pass
                    if name is not None:
                        print name.text
                    
                        try:
                            genus, gi, partial_variant = name.text.split("|")
                        except ValueError:
                            remove.append(clade)
                            continue

                        for variant in variants:
                            if variant in partial_variant:
                                break
                        else:
                            continue

                        name.attrib = {"bgStyle": variant.replace(".", "")}

                        if "canonical" in variant:
                            name.text = "."
                        else:
                            name.text = genus
                        
                        chart = ET.Element("chart")
                        group = ET.Element("group")
                        group.text = "markgroup{}".format(variant.replace(".", ""))
                        chart.append(group)
                        clade.append(chart)

                        annotation = ET.Element("annotation")
                        desc = ET.Element("desc")
                        desc.text = "Variant {} from {} ({})".format(variant, genus, gi)
                        uri = ET.Element("uri")
                        uri.text = "variant/{}/".format(variant)
                        annotation.append(desc)
                        annotation.append(uri)
                        clade.append(annotation)
                if remove:
                    for clade in remove:
                        parent = parent_map[clade]
                        child = parent.remove(clade)

            with open(os.path.join(self.trees_path, "{}.xml".format(hist_type)), "w") as outfile:
                treestr = StringIO.StringIO()
                tree.write(treestr)
                treestr = treestr.getvalue().replace("phy:", "")
                header, treestr = treestr.split("\n", 1)
                treestr = '<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">\n'+treestr
                outfile.write(treestr)
