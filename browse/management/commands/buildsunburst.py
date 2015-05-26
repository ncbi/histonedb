from django.core.management.base import BaseCommand, CommandError
from browse.models import *
import os
from itertools import chain
import pprint as pp
import json

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
            sb = self.build_sunburst(variant__core_type=core_histone)
            with open(os.path.join(path, "{}.json".format(core_histone.id)), "w") as core_burst:
                core_burst.write(sb)

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
        if filter.get("all_taxonomy"):
            taxas = Taxonomy.objects.filter(name="root", type_name="scientific name")
        else:
            taxas = Sequence.objects.filter(**filter).values_list("taxonomy", flat=True).distinct()

        sunburst = [{"name":"root", "children":[]}]
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
            root = sunburst[0]
            for i, curr_taxa in enumerate(path):
                #print "-"*i, curr_taxa.name, 
                for node in root["children"]:
                    if node["name"] == curr_taxa.name:
                        #print "(Exists)"
                        break
                    else:
                        #print "(DNE)"
                        pass
                else:
                    #print "(Adding)"
                    node = {"name":curr_taxa.name, "children":[]}
                    root["children"].append(node)
                root = node

            #pp.pprint(sunburst)

        return json.dumps(sunburst)
