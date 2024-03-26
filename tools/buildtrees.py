from itertools import cycle
import io
import subprocess

from browse.models import *

from Bio import SeqIO
from Bio import Phylo

import xml.etree.ElementTree as ET

ET.register_namespace("", "http://www.phyloxml.org/1.10/phyloxml.xsd")

colors = cycle([
    "#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
    "#a6d854",
    "#ffd92f",
    "#e5c494"])


class BuildTrees(object):
    help = 'Build the HistoneDB by training HMMs with seed sequences found in seeds directory in the top directory of thi project and using those HMMs to search the NR database.'
    seed_directory = os.path.join("static", "browse", "seeds")
    trees_path = os.path.join("static", "browse", "trees")

    def __init__(self, *args, **options):
        # self.make_trees()
        self.add_features()

    def get_variants(self, core_type=None):
        for i, (root, _, files) in enumerate(os.walk(self.seed_directory)):
            if core_type and not os.path.basename(root) == core_type: continue
            if i == 0:
                # Skip parent directory, only allow variant hmms to be built/searched
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
            print("Creating tree for", core_histone)
            if i == 0:
                # Skip parent directory, only allow variant hmms to be built/searched
                continue

            # Combine all varaints for a core histone type into one unaligned fasta file
            combined_seed_file = os.path.join(self.trees_path, "{}.fasta".format(core_histone))
            combined_seed_aligned = os.path.join(self.trees_path, "{}_aligned.fasta".format(core_histone))
            with open(combined_seed_file, "w") as combined_seed:
                for seed in files:
                    if not seed.endswith(".fasta"): continue
                    for s in SeqIO.parse(os.path.join(self.seed_directory, core_histone, seed), "fasta"):
                        s.seq = s.seq.replace("-", "")
                        SeqIO.write(s, combined_seed, "fasta")

            tree = os.path.join(self.trees_path, "{}_aligned.ph".format(core_histone))
            subprocess.call(["muscle", "-in", combined_seed_file, '-out', combined_seed_aligned])
            subprocess.call(["clustalw2", "-infile={}".format(combined_seed_aligned), '-tree'])
            Phylo.convert(tree, 'newick',
                          os.path.join(self.trees_path, "{}_no_features.xml".format(core_histone)), 'phyloxml')

    def add_features(self):
        for core_histone in ["H2A", "H2B", "H3", "H1", "H4"]:
            print(core_histone)
            tree_path = os.path.join(self.trees_path, "{}_no_features.xml".format(core_histone))
            print(tree_path)
            tree = ET.parse(tree_path)
            print(tree)
            parent_map = {c: p for p in tree.getiterator() for c in p}

            for phylogeny in tree.iter("{http://www.phyloxml.org}phylogeny"):
                print(phylogeny)
                render = ET.Element("render")

                parameters = ET.Element("parameters")
                circular = ET.Element("circular")
                radius = ET.Element("bufferRadius")
                radius.text = "0.5"
                circular.append(radius)
                parameters.append(circular)
                render.append(parameters)

                charts = ET.Element("charts")
                group = ET.Element("group", attrib={"type": "integratedBinary", "thickness": "20"})
                charts.append(group)
                render.append(charts)

                styles = ET.Element("styles")
                for variant in Variant.objects.filter(core_type=core_histone):
                    color = next(colors)
                    background = ET.Element("{}".format(variant.replace(".", "")),
                                            attrib={"fill": color, "stroke": color})

                    if not variant.id.startswith(core_type.id):
                        # Remove descriptor
                        start, name = variant.id[:variant.id.find(core_type.id)], variant.id[
                                                                                  variant.id.find(core_type.id):]
                        if len(start) > 3:
                            start = start[0]
                        name = start + name

                    label = ET.Element("markgroup{}".format(variant.replace(".", "")),
                                       attrib={"fill": "#000", "stroke": "#000", "opacity": "0.7", "label": name,
                                               "labelStyle": "sectorHighlightText"})
                    styles.append(background)
                    styles.append(label)
                label_sector = ET.Element("sectorHighlightText",
                                          attrib={"font-family": "Verdana", "font-size": "14", "font-weight": "bold",
                                                  "fill": "#FFFFFF", "rotate": "90"})
                styles.append(label_sector)
                render.append(styles)
                phylogeny.insert(0, render)

                for clade in phylogeny.iter("{http://www.phyloxml.org}clade"):
                    name = clade.find("{http://www.phyloxml.org}name")
                    try:
                        print(name.text)
                    except:
                        pass
                    if name is not None:
                        print(name.text)

                        try:
                            genus, gi, variant = name.text.split("|")
                            name.attrib = {"bgStyle": variant.replace(".", "")}
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
                        except ValueError:
                            parent = parent_map[clade]
                            child = parent.remove(clade)

            with open(os.path.join(self.trees_path, "{}.xml".format(core_histone)), "w") as outfile:
                treestr = io.BytesIO()
                tree.write(treestr)
                treestr = treestr.getvalue().decode("utf-8").replace("phy:", "")
                header, treestr = treestr.split("\n", 1)
                treestr = '<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">\n' + treestr
                outfile.write(treestr)


BuildTrees()
