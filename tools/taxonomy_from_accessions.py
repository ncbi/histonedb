import re
from Bio import Entrez, SeqIO

# *Always* tell NCBI who you are
Entrez.email = "l.singh@intbio.org"
from djangophylocore.models import Taxonomy
from browse.models import Sequence
from itertools import ifilter

import sys
import logging


def fetch_seq(accessions):
    if len(accessions) == 0:
        data = []
    else:
        for i in range(10):
            try:
                # print gis
                post_results = Entrez.read(Entrez.epost("protein", id=",".join(accessions)))
                webenv = post_results["WebEnv"]
                query_key = post_results["QueryKey"]
                handle = Entrez.efetch(db="protein", rettype="gb", retmode="text", webenv=webenv, query_key=query_key)
                data = list(SeqIO.parse(handle, "gb"))
                # print(data)
                if (len(accessions) == len(data)):
                    break
            except:
                logging.error("Unexpected error: {}".format(sys.exc_info()[0]))
                # continue
    for s in data:
        yield s


def taxonomy_from_accessions(accessions):
    """
    """
    for s in fetch_seq(accessions):
        logging.info(s.annotations["organism"])
        yield s.annotations["organism"]


def fetch_taxids(accessions):
    """
    """
    for s in fetch_seq(accessions):
        try:
            for a in s.features[0].qualifiers['db_xref']:
                text = re.search('(\S+):(\S+)', a).group(1)
                id = re.search('(\S+):(\S+)', a).group(2)
                if (text == "taxon"):
                    logging.info("Fetched taxid from NCBI {}".format(id))
                    yield id
                else:
                    continue
        except:
            logging.error("!!!!!!Unable to get TAXID for \n {} setting it to 1".format(s))
            yield 1  # unable to identify


from Bio import Entrez, SeqIO

already_exists = []


def taxonomy_from_header(header_init, accession=None, species_re=None):
    header = header_init.replace("_", " ")
    if species_re is None:
        species_re = re.compile(r'\[(.*?)\]')
    match = species_re.findall(header)
    if match:
        organism = match[-1]
    elif accession:
        logging.info("No taxonomy match for {}: {}, get it from NCBI".format(accession, header))
        for i in xrange(10):
            try:
                organism = taxonomy_from_accessions([accession]).next()
                break
            except StopIteration:
                pass
        else:
            return Taxonomy.objects.get(name="unidentified")
    else:
        return Taxonomy.objects.get(name="unidentified")
    organism = organism.replace(":", " ")
    organism = re.sub(r'([a-zA-Z0-9]+)[\./#_:]([a-zA-Z0-9]+)', r'\1 \2', organism)
    organism = re.sub(r'([a-zA-Z0-9]+)\. ', r'\1 ', organism)
    organism = re.sub(r"['\(\)\.]", r'', organism)
    organism = organism.lower()
    try:
        return Taxonomy.objects.get(name=organism.lower())
    except:
        try:
            genus = Taxonomy.objects.get(name=organism.split(" ")[0].lower())
            if genus.type_name != "scientific name":
                genus = genus.get_scientific_names()[0]
        except:
            logging.info(header)
            return Taxonomy.objects.get(name="unidentified")

        # Maybe the var is wrong
        if "var" in organism:
            try:
                org_end = organism.split("var ")[1]
            except IndexError:
                return Taxonomy.objects.get(name="unidentified")
            best_orgs = genus.children.filter(name__endswith=org_end).all()
            if len(best_orgs):
                return best_orgs[0]
            else:
                return Taxonomy.objects.get(name="unidentified")
        else:
            return Taxonomy.objects.get(name="unidentified")


easyspecies_re = re.compile(r'^(.*?)\|')


def easytaxonomy_from_header(header):
    """
  same as previous but different regex rule and no gi
    """
    return taxonomy_from_header(header, species_re=easyspecies_re)


def update_taxonomy(accessions):
    for taxid, accession in zip(fetch_taxids(accessions), accessions):
        seq = Sequence.objects.get(pk=accession)
        seq.taxonomy_id = taxid
        seq.save()


def get_tax_for_gi_taxdump(gis):
    with open("gi_taxid_prot.dmp") as gi2taxid:
        for line in ifilter(lambda l: l.split()[0] in gis, gi2taxid):
            yield list(reversed(line.strip().split()))  # Line is (gi, taxid)

