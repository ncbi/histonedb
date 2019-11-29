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
import logging

class Command(BaseCommand):
    help = 'Build the sunburst json files for each core histone and its variants'

    # Logging info
    logging.basicConfig(filename='log_buildsunburst.log',
                        format='%(asctime)s %(name)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    log = logging.getLogger(__name__)

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
        self.log.info('=======================================================')
        self.log.info('===               buildsunburst START               ===')
        self.log.info('=======================================================')
        path = os.path.join("static", "browse", "sunbursts")
        if options["all_taxonomy"]:
            sb = self.build_sunburst(all_taxonomy=True)
            with open(os.path.join(path, "all_taxa.json"), "w") as all_taxa:
                all_taxa.write(sb)

        for hist_type in Histone.objects.exclude(id="Unknown"):
            self.log.info("Saving {}".format(hist_type.id))

            vpath = os.path.join(path, hist_type.id)
            if not os.path.exists(vpath):
                os.makedirs(vpath)

            for variant in Variant.objects.filter(hist_type=hist_type):
                self.log.info("Saving {}".format(variant.id))
                sb = self.build_sunburst(variant=variant)
                with open(os.path.join(vpath, "{}.json".format(variant.id)), "w") as variant_burst:
                    variant_burst.write(sb)

        self.log.info('=======================================================')
        self.log.info('===       buildsunburst SUCCESSFULLY finished       ===')
        self.log.info('=======================================================')

    def build_sunburst(self, **filter):
        """Build the sunburst
        """
        sequences = Sequence.objects.filter(**filter).filter(all_model_scores__used_for_classification=True).annotate(score=Avg("all_model_scores__score"))
        return json.dumps(build_sunburst(sequences))

def build_sunburst(sequences):
    from djangophylocore.models import TaxonomyReference
    import networkx as nx
    from networkx.readwrite import json_graph

    scores = sequences.values_list("score", flat=True)
    scores_min = min(scores) if scores else 0
    scores_max = max(scores) if scores else 100

    green = Color("#66c2a5")
    red = Color("#fc8d62")
    color_range = list(red.range_to(green, 100))

    def get_color_for_taxa(taxon): 
        avg_score = taxon.children.filter(sequence__all_model_scores__used_for_classification=True).aggregate(score=Avg("sequence__all_model_scores__score"))["score"]
        avg_score = avg_score if avg_score else scores_min
        scaled = int(floor((float(avg_score-scores_min)/float(scores_max-scores_min))*100))
        color_index = scaled if scaled <= 99 else 99
        color_index = color_index if color_index >= 0 else 0
        return str(color_range[color_index])

    taxa = list(sequences.values_list("taxonomy__parent__parent__parent", flat=True).distinct())
    allow_ranks = ["kingdom", "phylum", "order"]
    tree = TaxonomyReference().get_filtered_reference_graph(taxa, allow_ranks=allow_ranks)
    nx.set_node_attributes(tree, "colour", {n:get_color_for_taxa(Taxonomy.objects.get(name=n)) for n,d in tree.out_degree_iter() if d==0})
    
    return json_graph.tree_data(tree, "root", attrs={'children': 'children', 'id': 'name'})
