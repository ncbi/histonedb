import os
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from browse.models import Variant, OldStyleVariant
from djangophylocore.models import Taxonomy
import json

class Command(BaseCommand):
    help = 'Reset sequence features'
    info_directory = os.path.join(settings.STATIC_ROOT, "browse", "info")

    def add_arguments(self, parser):
         parser.add_argument(
            "-f", 
            "--force", 
            default=False, 
            action="store_true", 
            help="Force the creation of PDFs, GFFs even if the files exist")
        
    def handle(self, *args, **options):
        variant_info_path = os.path.join(self.info_directory, "variants.json")

        with open(variant_info_path) as variant_info_file:  
            variant_info = json.load(variant_info_file)

        for core_histone_name, variants in variant_info.iteritems():
            for variant_name, info in variants.iteritems():
                variant = Variant.objects.get(id=variant_name)
                variant.description = info["description"]
                variant.taxonmic_span = info["taxonmic_span"]
                variant.save()

                for alternate_name in info["alternate_names"]:
                    tax_name = alternate_name.get("taxonomy", "eukaryotes")
                    print tax_name
                    alt_variant = OldStyleVariant(
                        updated_variant = variant,
                        name            = alternate_name["name"],
                        gene            = alternate_name.get("gene"),
                        splice          = alternate_name.get("splice"),
                        taxonomy        = Taxonomy.objects.get(name=tax_name)
                    )
                    alt_variant.save()
