import os
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from browse.models import Variant, OldStyleVariant, Histone, Publication
from djangophylocore.models import Taxonomy
import json

class Command(BaseCommand):
    help = 'Reset sequence features'
    info_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "info")

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

        for hist_type_name, variants in variant_info.iteritems():
            for variant_name, info in variants.iteritems():
                variant = Variant.objects.get(id=variant_name)
                variant.description = info["description"]
                variant.taxonomic_span = info["taxonomic_span"]
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
                    try:
                        alt_variant.save()
                    except:
                        from django.db import connection
                        cursor = connection.cursor()
                        cursor.execute("ALTER DATABASE histdb CHARACTER SET utf8 COLLATE utf8_general_ci")
                        alt_variant.save()

                for publication_id in info["publications"]:
                    publication, created = Publication.objects.get_or_create(id=publication_id, cited=False)
                    publication.variants.add(variant)

        type_info_path = os.path.join(self.info_directory, "types.json")

        with open(type_info_path) as type_info_file:  
            type_info = json.load(type_info_file)


        for type_name, info in type_info.iteritems():
            type = Histone.objects.get(id=type_name)
            type.description = info["description"]
            type.save()


        feature_info_path = os.path.join(self.info_directory, "features.json")
        with open(feature_info_path) as feature_info_file:  
            feature_info = json.load(feature_info_file)

        
        #Features not working yet
        for type_name, variants in feature_info.iteritems():
            type = Histone.objects.get(id=type_name)
            for variant_name, info in variants.iteritems():
                if variant_name == "General":
                    variant = type
                else:
                    variant = Variant.objects.get(id=variant_name)

                try:
                    lines = [l.rstrip() for l in info["format"].split("\n") if l.rstrip()]
                    sequence = lines[0]
                    position_lines = lines[1:]
                except IndexError:
                    raise RuntimeError("Invalid feature format")

                assert False not in [len(position)==len(sequence) for position in position_lines]

                try:
                    taxonomy = Taxononomy.objects.get(name=info.get("taxonomy", "unidentified"))
                except Taxononomy.DoesNotExist:
                    taxonomy = Taxononomy.objects.get(name="unidentified")

                template, created = TemplateSequence.objects.get_or_create(taxonomy=taxonomy, variant=variant.id)
                if created:
                    SeqIO.write(
                        SeqRecord(Sequence(sequence), id="{}_{}".format(type_name, variant)),
                        template.path()
                    )

                for positions in position_lines:
                    for feature_name, group in groupby(enumerate(positions), key=lambda x:x[1]):
                        if not feature_name in [" ", "="]:
                            Feature(
                                template    = template,
                                start       = group[0][0],
                                end         = group[-1][0],
                                name        = info["features"][feature_name]["name"],
                                description = info["features"][feature_name]["description"],
                            )
                            Feature.save()
