from django.core.management.base import BaseCommand

import os
from typing import Set
from django.conf import settings
from tqdm import tqdm
from djangophylocore.models import *


class Command(BaseCommand):
    help = "Load all taxonomy data into database using Django ORM"
    inserted = set()
    objects = []

    requires_system_checks = [True]

    def inserting_rank(self, path_dumps):
        print("Inserting Taxonomy Rank")
        self.objects = []
        with open(os.path.join(path_dumps, "rank.dmp")) as fin:
            for line in fin.readlines():
                fields = line.strip().split("|")
                self.objects.append(Rank(id=int(fields[0]), name=fields[1]))
        Rank.objects.bulk_create(self.objects)
        self.objects = []

    def insert_taxonomy(self, taxonomies: dict, tax_id: int):
        taxonomy = taxonomies[tax_id]
        if taxonomy.id != taxonomy.parent_id:
            if int(taxonomy.parent_id) not in self.inserted:
                self.insert_taxonomy(taxonomies, taxonomy.parent_id)
        if tax_id not in self.inserted:
            self.inserted.add(tax_id)
            self.objects.append(taxonomy)
            if len(self.objects) == 500:
                Taxonomy.objects.bulk_create(self.objects)
                self.objects = []

    def inserting_taxonomy(self, path_dumps):
        print("Loading Taxonomies")
        with open(os.path.join(path_dumps, "taxonomy.dmp")) as fin:
            taxonomies = {}
            for line in tqdm(fin.readlines()):
                fields = line.strip().split("|")
                taxonomies[int(fields[0])] = Taxonomy(id=int(fields[0]),
                                                      name=fields[1],
                                                      type_name=fields[2],
                                                      parent_id=int(fields[3]),
                                                      rank_id=int(fields[4]))
            print("\nInserting Taxonomies")
            for taxid, taxonomy in tqdm(taxonomies.items()):
                self.insert_taxonomy(taxonomies, taxid)
            if self.objects:
                Taxonomy.objects.bulk_create(self.objects)
            self.objects = []

    def insert_parent_relation(self, path_dumps):
        print("Loading parent relationship")
        self.objects = []
        with open(os.path.join(path_dumps, "parentsrelation.dmp")) as fin:
            for line in tqdm(fin.readlines()):
                fields = line.strip().split("|")
                self.objects.append(ParentsRelation(id=int(fields[0]),
                                                    index=int(fields[1]),
                                                    taxon_id=int(fields[3]),
                                                    parent_id=int(fields[2])))
                if len(self.objects) == 10000:
                    ParentsRelation.objects.bulk_create(self.objects)
                    self.objects = []
        if self.objects:
            ParentsRelation.objects.bulk_create(self.objects)
            self.objects = []

    def insert_common_taxa(self, path_dumps):
        print("Loading common taxa")
        self.objects = []
        with open(os.path.join(path_dumps, "relcommontaxa.dmp")) as fin:
            for line in tqdm(fin.readlines()):
                fields = line.strip().split("|")
                self.objects.append(RelCommonTaxa(id=int(fields[0]),
                                                  language=fields[1],
                                                  taxon_id=int(fields[2]),
                                                  common_id=int(fields[3])))
                if len(self.objects) == 10000:
                    RelCommonTaxa.objects.bulk_create(self.objects)
                    self.objects = []
        if self.objects:
            RelCommonTaxa.objects.bulk_create(self.objects)
            self.objects = []

    def insert_relsynonymtaxa(self, path_dumps):
        print("Loading relation synonym taxa")
        self.objects = []
        with open(os.path.join(path_dumps, "relsynonymtaxa.dmp")) as fin:
            for line in tqdm(fin.readlines()):
                fields = line.strip().split("|")
                self.objects.append(RelSynonymTaxa(id=int(fields[0]),
                                                   taxon_id=int(fields[1]),
                                                   synonym_id=int(fields[2])))
                if len(self.objects) == 10000:
                    RelSynonymTaxa.objects.bulk_create(self.objects)
                    self.objects = []
        if self.objects:
            RelSynonymTaxa.objects.bulk_create(self.objects)
            self.objects = []

    def insert_relhomonymtaxa(self, path_dumps):
        print("Loading relation synonym taxa")
        self.objects = []
        with open(os.path.join(path_dumps, "relhomonymtaxa.dmp")) as fin:
            for line in tqdm(fin.readlines()):
                fields = line.strip().split("|")
                self.objects.append(RelHomonymTaxa(id=int(fields[0]),
                                                   taxon_id=int(fields[1]),
                                                   homonym_id=int(fields[2])))
                if len(self.objects) == 10000:
                    RelHomonymTaxa.objects.bulk_create(self.objects)
                    self.objects = []
        if self.objects:
            RelHomonymTaxa.objects.bulk_create(self.objects)
            self.objects = []

    def handle(self, **options):
        localDir = os.path.dirname(__file__)
        absDir = os.path.join(os.getcwd(), localDir)
        verbose = options.get("verbose", True)
        if verbose:
            print("loading taxonomy, please wait, it can take a while...")
        path_dumps = os.path.join(absDir, '..', '..', 'dumps')
        self.inserting_rank(path_dumps)
        self.inserting_taxonomy(path_dumps)
        self.insert_parent_relation(path_dumps)
        self.insert_common_taxa(path_dumps)
        self.insert_relsynonymtaxa(path_dumps)
        self.insert_relhomonymtaxa(path_dumps)
