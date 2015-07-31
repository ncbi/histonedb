from django.core.management.base import BaseCommand, CommandError
from browse.models import *
import os
from itertools import chain
import pprint as pp
import json
from colour import Color
from django.db.models import Max, Min, Count
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
        green = Color("#66c2a5")
        red = Color("#fc8d62")
        color_range = list(red.range_to(green, 100))

        taxas = Sequence.objects.filter(**filter).filter(all_model_scores__used_for_classifiation=True).annotate(score=Max("all_model_scores__score"))
        scores_min = min(taxas.values_list("score", flat=True) or [0])
        scores_max = max(taxas.values_list("score", flat=True) or [0])

        taxas = taxas.values_list("taxonomy", flat=True).distinct()

        def get_color_for_taxa(taxon):
            ids = set()
            ids.add(taxon.id)
            children = set(taxon.children.values_list("id", flat=True))
            ids |= children
            scores = list(Sequence.objects.filter(taxonomy__in=list(ids)).filter(all_model_scores__used_for_classifiation=True).annotate(score=Max("all_model_scores__score")).values_list("score", flat=True))
            avg = float(sum(scores))/len(scores) if len(scores) > 0 else 0.
            scaled = int(floor((float(avg-scores_min)/float(scores_max-scores_min))*100))
            color_index = scaled if scaled <= 99 else 99
            color_index = color_index if color_index >= 0 else 0
            return str(color_range[color_index])

        sunburst = {"name":"root", "children":[]}
        colors = {"eukaryota":"#6600CC", "prokaryota":"#00FF00", "archea":"#FF6600"}
        for taxa in taxas:
            if filter.get("all_taxonomy"):
                path = chain([taxa], taxa.children.all())
            else:
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
                try:
                    if curr_taxa.rank.name in ["no rank", "class"] or "sub" in curr_taxa.rank.name or  "super" in curr_taxa.rank.name:
                        continue
                except:
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

            #pp.pprint(sunburst)

        return json.dumps(sunburst)
