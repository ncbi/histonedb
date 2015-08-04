from django.core.management.base import BaseCommand, CommandError
from browse.models import *
from djangophylocore.models import Rank
import os
from itertools import chain
import pprint as pp
import json
from colour import Color
from django.db.models import Max, Min, Count, Avg
from math import floor

class Command(BaseCommand):
    help = 'Build the sunburst json files for each core histone and its variants'

    def add_arguments(self, parser):
        parser.add_argument(
            "-f", 
            "--force", 
            default=False, 
            action="store_true", 
            help="Force the recreation of sunburst. This will overwrite any existing sunburst json files.")
        parser.add_argument(
            "--all_taxonomy",
            default=False,
            action="store_true",
            help="Build a sunburst of the enitre NCBI Taxonomy databse.")
        
    def handle(self, *args, **options):
        path = os.path.join("static", "browse", "sunbursts")
        if options["all_taxonomy"]:
            sb = self.build_sunburst(all_taxonomy=True)
            with open(os.path.join(path, "all_taxa.json"), "w") as all_taxa:
                all_taxa.write(sb)

        for core_histone in Histone.objects.all():
            print "Saving", core_histone.id
            #sb = self.build_sunburst(variant__core_type=core_histone)
            #with open(os.path.join(path, "{}.json".format(core_histone.id)), "w") as core_burst:
            #    core_burst.write(sb)

            vpath = os.path.join(path, core_histone.id)
            if not os.path.exists(vpath):
                os.makedirs(vpath)

            for variant in Variant.objects.filter(core_type=core_histone):
                print "  Saving", variant.id
                sb = self.build_sunburst(variant=variant)
                with open(os.path.join(vpath, "{}.json".format(variant.id)), "w") as variant_burst:
                    variant_burst.write(sb)

    def build_sunburst(self, **filter):
        """Build the sunburst
        """
        sequences = Sequence.objects.filter(**filter).filter(all_model_scores__used_for_classifiation=True).annotate(score=Avg("all_model_scores__score"))
        return json.dumps(build_sunburst(sequences))

def build_sunburst(sequences):
    scores = sequences.values_list("score", flat=True)
    scores_min = min(scores) if scores else 0
    scores_max = max(scores) if scores else 100

    green = Color("#66c2a5")
    red = Color("#fc8d62")
    color_range = list(red.range_to(green, 100))

    def get_color_for_taxa(taxon): 
        #ids = set()
        #ids.add(taxon.id)
        #children = set(taxon.children.values_list("id", flat=True))
        #ids |= children
        avg_score = taxon.children.filter(sequence__all_model_scores__used_for_classifiation=True).aggregate(score=Avg("sequence__all_model_scores__score"))["score"]
        avg_score = avg_score if avg_score else scores_min
        scaled = int(floor((float(avg_score-scores_min)/float(scores_max-scores_min))*100))
        color_index = scaled if scaled <= 99 else 99
        color_index = color_index if color_index >= 0 else 0
        return str(color_range[color_index])

    taxa = sequences.values_list("taxonomy__parent__parent__parent", flat=True).distinct()

    from djangophylocore.models import TaxonomyReference
    reference = TaxonomyReference()
    import networkx as nx
    allow_ranks = set(Rank.objects.filter(name__in=["kingdom", "phylum", "order", "tribe", "forma"]))
    tree = reference.get_filtered_reference_graph(taxa, allow_ranks=allow_ranks)
    from networkx.readwrite import json_graph
    import networkx as NX
    NX.set_node_attributes(tree, "colour", {n:get_color_for_taxa(Taxonomy.objects.get(name=n)) for n,d in tree.out_degree_iter() if d==0})
    
    return json_graph.tree_data(tree, "root", attrs={'children': 'children', 'id': 'name'})


    sunburst = {"name":"root", "children":[]}
    paths = {}
    curr_path = paths
    for taxa in taxas:
        taxa = Taxonomy.objects.get(pk=taxa)
        try:
            taxa = taxa.get_scientific_names()[0]
        except IndexError:
            pass
        path = chain(taxa.parents.reverse().all()[1:], [taxa])
        root = sunburst
        stopForOrder = False
        
        for i, curr_taxa in enumerate(path):
            if stopForOrder:
                break
            if curr_taxa.rank.name in ["no rank", "class", "sub", "super"] or "sub" in curr_taxa.rank.name or  "super" in curr_taxa.rank.name:
                continue
            if "order" in curr_taxa.rank.name or "family" in curr_taxa.rank.name:
                stopForOrder = True
            try:
                curr_path = curr_path[curr_taxa.name]
            except KeyError:
                curr_path[curr_taxa.name] = {}
                curr_path = curr_path[curr_taxa.name]

        continue
        for i, curr_taxa in enumerate(path):
            if stopForOrder:
                break
            if curr_taxa.rank.name in ["no rank", "class"] or "sub" in curr_taxa.rank.name or  "super" in curr_taxa.rank.name:
                continue
            if "order" in curr_taxa.rank.name or "family" in curr_taxa.rank.name:
                stopForOrder = True
            #print "-"*(i+1), curr_taxa.name, curr_taxa.rank.name
            for node in root["children"]:
                if node["name"] == curr_taxa.name:
                    break
                else:
                    pass
            else:
                #print "(Adding)"
                node = {"name":curr_taxa.name, "children":[]}
                root["children"].append(node)
            root = node
        try:
            root["colour"] = get_color_for_taxa(taxa)
        except:
            root["color"] = str(red)

    assert 0, paths
    return sunburst
