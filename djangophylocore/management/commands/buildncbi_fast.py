# UNDER DEVELOPMENT
from django.core.management.base import BaseCommand

import os
import sys
from django.conf import settings
from tqdm import tqdm
from joblib import Parallel, delayed

localDir = os.path.dirname(__file__)
absDir = os.path.join(os.getcwd(), localDir)
DUMP_PATH = os.path.join(absDir, '..', '..', 'dumps')

global TBI
global TBN

TBI = {}
TBN = {}


def loop1(line):
    id = line.split("|")[0].strip()
    name = line.split("|")[1].strip().lower()

    name = clean_name(name)  # VR sept 09 clean name
    homonym = line.split("|")[2].strip().lower()
    if homonym:
        homonym = clean_homonym(name, homonym)  # VR sept 09 clean name
    # bad way to handle homonym not declare by NCBI e.g. Influenza A virus (A/turkey/Beit_Herut/1265/03(H9N2))
    if not homonym and name in TBN:
        name = name + " ncbiid " + id
    type_name = line.split("|")[3].strip()
    synonym = "synonym" in type_name
    common = "common name" in type_name
    # self.index = int(id)
    if type_name == "scientific name":
        # Creating TAXONOMY_BY_ID
        TBI[id] = {}
        if homonym:
            TBI[id]["name"] = homonym
        else:
            TBI[id]["name"] = name
        TBI[id]["common"] = []
        TBI[id]["homonym"] = []
        TBI[id]["synonym"] = []
        TBI[id]["parent"] = []
        TBI[id]["parents"] = []
        TBI[id]["type_name"] = 'scientific name'
        # Creating TAXONOMY_BY_NAME
        if homonym:
            TBN[homonym] = {}
            TBN[homonym]["id"] = id
            TBN[homonym]["homonym"] = name
        else:
            TBN[name] = {}
            TBN[name]["id"] = id
    return (int(id))


class Command(BaseCommand):
    help = "Download and build the ncbi database"
    requires_system_checks = []

    TEST = False
    # TBI = {}
    # TBN = {}
    NAMES = "names.dmp"
    NODES = "nodes.dmp"

    def handle(self, **options):
        global DUMP_PATH
        self.list_id = []
        verbose = options.get("verbose", True)
        if verbose:
            print("loading taxonomy, please wait, it can take a while...")
        if not os.path.exists(DUMP_PATH):
            os.system('mkdir %s' % DUMP_PATH)
        else:
            os.system('rm %s/*' % DUMP_PATH)
        # VR 15 aout 2009 otherwise the file is always loaded from the web
        os.system('rm taxdump.tar.gz')
        self.download_ncbi(verbose)
        self.generate_structure(verbose)
        if verbose:
            print("making rank.dmp")
        self.make_rank()
        if verbose:
            print("making taxonomy.dmp")
        self.make_taxa()
        self.make_taxonomy_plus(verbose)
        if verbose:
            print("making parents")
        self.make_parents()

    # VRaou09 debug
    #        os.system('rm nodes.dmp names.dmp')
    #        os.system('rm taxdump.tar.gz')

    def download_ncbi(self, verbose):
        if not os.path.exists('./taxdump.tar.gz'):
            if verbose:
                print("Downloading NCBI database on the web")
            os.system("curl -# ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz > taxdump.tar.gz ")
            # os.system("wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
        if verbose:
            print("Extracting database... please wait")
        os.system("tar xf taxdump.tar.gz names.dmp nodes.dmp")

    def generate_structure(self, verbose):
        global DUMP_PATH

        def getParents(id, TBI):
            lp = []
            while id != "1":
                id_parent = TBI[id]["parent"]
                lp.append(id_parent)
                id = id_parent
            return lp

        if verbose:
            print("Generating structure...")
        # Retrieving all scientific names
        index = 0
        # Parallel(n_jobs=2)(delayed(loop1)(line) for line in tqdm(open(self.NAMES).readlines()))
        for line in tqdm(open(self.NAMES).readlines()):
            self.index = loop1(line)
        if verbose:
            print("Extracting parents...")
        for node in tqdm(open(self.NODES).readlines()):
            id = node.split("|")[0].strip()
            parent = node.split("|")[1].strip()
            rank = node.split('|')[2].strip()
            name = TBI[id]["name"]
            name = TBI[id]['rank'] = rank
            if id not in TBI:
                TBI[id] = {}
            TBI[id]["parent"] = parent
            if name not in TBN:
                TBN[name] = {}
            TBN[name]["parent"] = parent
        if verbose:
            print("Filling parents...")
        for node in tqdm(open(self.NODES).readlines()):
            id = node.split("|")[0].strip()
            TBI[id]["parents"] = getParents(id, TBI)

    def make_taxonomy_plus_old(self, verbose):
        if verbose:
            print("Adding synonyms, homonyms and common names...")
        # Adding synonyms, homonyms and common names
        index = self.index
        list_synonym = []
        list_common = []
        list_homonym = []
        list_relsynonymtaxa = []
        list_relcommontaxa = []
        list_relhomonymtaxa = []
        index_relsynonym = 0
        index_relcommon = 0
        index_relhomonym = 0
        homonym_toc = {}
        common_toc = {}
        synonym_toc = {}
        for line in tqdm(open(self.NAMES).readlines()):
            type_name = line.split("|")[3].strip()
            synonym = "synonym" in type_name
            common = "common name" in type_name
            homonym = line.split("|")[2].strip().lower()
            id = line.split("|")[0].strip()
            if self.TEST:
                if id not in self.list_id:
                    continue
            name = line.split("|")[1].strip().lower()
            name = name.replace(")", " ").replace("(", " ").replace(",", " ").replace(":", " ").replace(";",
                                                                                                        " ").replace(
                "'", " ");
            if synonym or common:
                if homonym:  # We do not want synonym wich have homonym
                    continue
                base_name = TBI[id]["name"]
                if synonym:
                    if name not in synonym_toc:
                        index += 1
                        list_synonym.append("%s|%s|synonym|2|\n" % (index, name))
                        synonym_toc[name] = index
                    index_relsynonym += 1
                    list_relsynonymtaxa.append(
                        "%s|%s|%s\n" % (index_relsynonym, index, id))
                if common:
                    if name not in common_toc:
                        index += 1
                        list_common.append("%s|%s|common|2|\n" % (index, name))
                        common_toc[name] = index
                    index_relcommon += 1
                    list_relcommontaxa.append(
                        "%s|%s|%s|english\n" % (index_relcommon, index, id))
            if type_name == "scientific name" and homonym:
                if name not in homonym_toc:
                    index += 1
                    list_homonym.append('%s|%s|homonym|2|\n' % (index, name))
                    homonym_toc[name] = index
                index_relhomonym += 1
                list_relhomonymtaxa.append(
                    "%s|%s|%s\n" % (index_relhomonym, homonym_toc[name], id))
        taxonomy_file = open(os.path.join(DUMP_PATH, 'taxonomy.dmp'), 'a')
        taxonomy_file.write(''.join(list_synonym))
        taxonomy_file.write(''.join(list_common))
        taxonomy_file.write(''.join(list_homonym))
        taxonomy_file.close()
        open(os.path.join(DUMP_PATH, 'relsynonymtaxa.dmp'), 'w').write(
            ''.join(list_relsynonymtaxa))
        open(os.path.join(DUMP_PATH, 'relcommontaxa.dmp'), 'w').write(
            ''.join(list_relcommontaxa))
        open(os.path.join(DUMP_PATH, 'relhomonymtaxa.dmp'), 'w').write(
            ''.join(list_relhomonymtaxa))

    def make_taxonomy_plus(self, verbose):
        if verbose:
            print("Adding synonyms, homonyms and common names...")
        # Adding synonyms, homonyms and common names
        index = self.index
        list_synonym = []
        list_common = []
        list_homonym = []
        list_relsynonymtaxa = []
        list_relcommontaxa = []
        list_relhomonymtaxa = []
        index_relsynonym = 0
        index_relcommon = 0
        index_relhomonym = 0
        homonym_toc = {}
        common_toc = {}
        synonym_toc = {}

        for line in open(self.NAMES).readlines():
            type_name = line.split("|")[3].strip()
            synonym = "synonym" in type_name
            common = "common name" in type_name
            homonym = line.split("|")[2].strip().lower()
            id = line.split("|")[0].strip()
            if self.TEST:
                if id not in self.list_id:
                    continue
            name = line.split("|")[1].strip().lower()
            # name = name.replace(")", " ").replace("(", " ").replace(",", " ").replace(":", " ").replace(";", " ").replace("'", " ")
            name = clean_name(name)
            if synonym:
                if homonym:  # We do not want synonym wich have homonym
                    continue
                base_name = TBI[id]["name"]
                # we don't want synonyms similar to scientific names
                if name in TBN:
                    continue
                if synonym:
                    if name not in synonym_toc:
                        index += 1
                        list_synonym.append("%s|%s|synonym|%s|2\n" % (index, name, index))
                        synonym_toc[name] = index
                        # VR sept09 bug correction next line were out of the if when the name was encounter a second time the same line was added
                        index_relsynonym += 1
                        list_relsynonymtaxa.append(
                            "%s|%s|%s\n" % (index_relsynonym, id, index))
            if type_name == "scientific name" and homonym:
                if name not in homonym_toc:
                    index += 1
                    list_homonym.append('%s|%s|homonym|%s|2\n' % (index, name, index))
                    homonym_toc[name] = index
                index_relhomonym += 1
                list_relhomonymtaxa.append(
                    "%s|%s|%s\n" % (index_relhomonym, homonym_toc[name], id))

        # commons
        for line in open(self.NAMES).readlines():
            type_name = line.split("|")[3].strip()
            synonym = "synonym" in type_name
            common = "common name" in type_name
            homonym = line.split("|")[2].strip().lower()
            id = line.split("|")[0].strip()
            if self.TEST:
                if id not in self.list_id:
                    continue
            name = line.split("|")[1].strip().lower()
            # name= name.replace(")", " ").replace("(", " ").replace(",", " ").replace(":", " ").replace(";", " ").replace("'", " ");
            name = clean_name(name)
            if common:
                if homonym:  # We do not want synonym wich have homonym
                    continue
                base_name = TBI[id]["name"]
                if common:
                    if name not in common_toc and name not in list(
                            homonym_toc.keys()) and name not in synonym_toc and name not in TBN:
                        index += 1
                        list_common.append("%s|%s|common|%s|2\n" % (index, name, index))
                        common_toc[name] = index
                        # VR sept 09 the two following line should be in the if
                        index_relcommon += 1
                        list_relcommontaxa.append(
                            "%s|english|%s|%s\n" % (index_relcommon, id, index))

        taxonomy_file = open(os.path.join(DUMP_PATH, 'taxonomy.dmp'), 'a')
        taxonomy_file.write(''.join(list_synonym))
        taxonomy_file.write(''.join(list_common))
        taxonomy_file.write(''.join(list_homonym))
        taxonomy_file.close()
        open(os.path.join(DUMP_PATH, 'relsynonymtaxa.dmp'), 'w').write(
            ''.join(list_relsynonymtaxa))
        open(os.path.join(DUMP_PATH, 'relcommontaxa.dmp'), 'w').write(
            ''.join(list_relcommontaxa))
        open(os.path.join(DUMP_PATH, 'relhomonymtaxa.dmp'), 'w').write(
            ''.join(list_relhomonymtaxa))

    ################################################
    #               Generating dumps               #
    ################################################
    RANK = {}

    def make_rank(self):
        global DUMP_PATH
        l_rank = []
        list_line = []
        index = 0
        for species in tqdm(sorted(list(TBI.keys()), key=lambda x: int(x))):
            rank = TBI[species]['rank']
            if rank not in l_rank:
                index += 1
                line = '%s|%s' % (index, rank)
                list_line.append(line)
                l_rank.append(rank)
                self.RANK[rank] = index
        open(os.path.join(DUMP_PATH, 'rank.dmp'), 'w').write('\n'.join(list_line))

    def make_taxa(self):
        global DUMP_PATH
        # Taxa.dmp
        list_line = []
        for species in tqdm(sorted(list(TBI.keys()), key=lambda x: int(x))):
            if self.TEST:
                if species not in self.list_id:
                    continue
            line = '%s|%s|%s|%s|%s\n' % (
                species,
                # self.TBI[species]['name'].replace(")", " ").replace("(", " ").replace(",", " ").replace(":", " ").replace(";", " ").replace("'", " "),
                # self.clean_name(self.TBI[species]['name']),
                TBI[species]['name'],
                # VR sept09 .replace(")", " ").replace("(", " ").replace(",", " ").replace(":", " ").replace(";", " ").replace("'", " "),
                TBI[species]['type_name'],
                TBI[species]['parent'],
                self.RANK[TBI[species]['rank']]
            )
            list_line.append(line)
        open(os.path.join(DUMP_PATH, 'taxonomy.dmp'), 'w').write(''.join(list_line))

    def make_parents(self):
        global DUMP_PATH
        id_rel = 0
        list_parents = []
        for species in tqdm(sorted(list(TBI.keys()), key=lambda x: int(x))):
            index = 0
            if self.TEST:
                if species not in self.list_id:
                    continue
            for parent_id in TBI[species]["parents"]:
                id_rel += 1
                line = '%s|%s|%s|%s\n' % (id_rel, index, parent_id, species)
                list_parents.append(line)
                index += 1
        open(os.path.join(DUMP_PATH, 'parentsrelation.dmp'), 'w').write(''.join(list_parents))


def clean_homonym(name, homonym):
    cleanName = clean_name(name)
    cleanHomonyn = clean_name(homonym)
    homonymBegin = homonym.find("<")
    homonymEnd = homonym.find(">")
    if (homonymBegin > 0 and homonymEnd > 0):
        cleanHomonyn = cleanHomonyn[(homonymBegin + 1):homonymEnd]
    return cleanName + " <" + cleanHomonyn + ">"


def clean_name(name):
    cleanName = "";
    allowed = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890+-*<>";
    spaceOpen = False
    for l in name:
        if l in allowed:
            #            if l in "'":
            #                if spaceOpen:
            #                    l="pr "
            #                else:
            #                    l=" pr "
            #                spaceOpen=True;
            cleanName = cleanName + l;
            spaceOpen = False
        else:
            if not spaceOpen:
                spaceOpen = True
                cleanName = cleanName + " "
    cleanName = cleanName.strip()
    # cleanName = cleanName.replace(" ", "_")
    return cleanName
    # cleanNameNwk = name.replace(")", "_").replace("(", "_").replace(",", " ").replace(":", " ").replace(";", " ")
    # cleanName = cleanNameNwk.replace("'", " ").replace("`"," ").replace('"',' ').strip()
