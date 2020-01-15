import os
import logging
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from browse.models import Variant, OldStyleVariant, Histone, Publication, TemplateSequence, Feature
from djangophylocore.models import Taxonomy
import json
from Bio import SeqIO
from itertools import groupby

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Command(BaseCommand):
    help = 'Reset sequence features'
    info_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "info")

    # Logging info
    logging.basicConfig(filename='log/buildvariantinfo.log',
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
            help="Force the creation of PDFs, GFFs even if the files exist")
        
    def handle(self, *args, **options):
        self.log.info('=======================================================')
        self.log.info('===             buildvariantinfo START              ===')
        self.log.info('=======================================================')
        Publication.objects.all().delete()
        variant_info_path = os.path.join(self.info_directory, "variants.json")

        with open(variant_info_path) as variant_info_file:  
            variant_info = json.load(variant_info_file)

        for hist_type_name, variants in variant_info.iteritems():
            for variant_name, info in variants.iteritems():
                self.log.info(variant_name)
                variant = Variant.objects.get(id=variant_name)
                variant.description = info["description"]
                variant.taxonomic_span = info["taxonomic_span"]
                variant.save()

                for alternate_name in info["alternate_names"]:
                    tax_name = alternate_name.get("taxonomy", "eukaryota")
                    self.log.info(tax_name)
                    # tax_eu = Taxonomy.objects.get(name='eucaryotes')
                    # print('===================== {}'.format(tax_eu))
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
                        cursor.execute("ALTER TABLE browse_feature CONVERT TO CHARACTER SET utf8 COLLATE     utf8_general_ci;")
                        alt_variant.save()

                for publication_id in info["publications"]:
                    publication, created = Publication.objects.get_or_create(id=publication_id, cited=False)
                    publication.variants.add(variant)

        type_info_path = os.path.join(self.info_directory, "types.json")

        with open(type_info_path) as type_info_file:  
            type_info = json.load(type_info_file)


        for type_name, info in type_info.iteritems():
            hist_type = Histone.objects.get(id=type_name)
            hist_type.description = info["description"]
            hist_type.save()

        feature_info_path = os.path.join(self.info_directory, "features.json")
        with open(feature_info_path) as feature_info_file:  
            feature_info = json.load(feature_info_file)
        
        Feature.objects.all().delete()
        for type_name, variants in feature_info.iteritems():
            self.log.info("Making features for {}".format(type_name))
            hist_type = Histone.objects.get(id=type_name)
            for variant, info in variants.iteritems():
                self.log.info("Making features for {}".format(variant_name))
                if variant_name.startswith("General"):
                    variant = "General{}".format(hist_type)

                sequence = str(info["sequence"])
                position_lines = [str(position) for key, position in info.iteritems() if key.startswith("feature") and not key.endswith("info")]
                assert False not in [len(position)==len(sequence) for position in position_lines], "Sequence and feaures must have the same number of characters!\n{}\n{}".format(sequence, "\n".join(position_lines))

                try:
                    taxonomy = Taxonomy.objects.get(name=info.get("taxonomy", "root"))
                except Taxonomy.DoesNotExist:
                    taxonomy = Taxonomy.objects.get(name="root")

                template, created = TemplateSequence.objects.get_or_create(taxonomy=taxonomy, variant=variant)
                # if not os.path.isfile(template.path()): #we need to rewrite it!!!
                SeqIO.write(
                        SeqRecord(Seq(sequence), id=str(template)),
                        template.path(),
                        "fasta"
                )
                used_features = {}
                for positions in position_lines:
                    for feature_name, group in groupby(enumerate(positions), key=lambda x:x[1]):
                        group = list(group)
                        if not feature_name in [" ", "="]:
                            name = info["feature_info"][feature_name]["name"]
                            feature = Feature(
                                id          = "{}_{}{}".format(template, name, " "+str(used_features.get(name, ""))),
                                template    = template,
                                start       = int(group[0][0]),
                                end         = int(group[-1][0]),
                                name        = name,
                                description = info["feature_info"][feature_name]["description"],
                                color       = info["feature_info"][feature_name]["color"],
                            )
                            feature.save()
                            try:
                                used_features[name] += 1
                            except KeyError:
                                used_features[name] = 1

        self.log.info('=======================================================')
        self.log.info('===     buildvariantinfo SUCCESSFULLY finished      ===')
        self.log.info('=======================================================')
