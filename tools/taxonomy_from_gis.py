import re
from Bio import Entrez, SeqIO
# *Always* tell NCBI who you are
Entrez.email = "eli.draizen@nih.gov"
from djangophylocore.models import Taxonomy

def seq_from_gi(gis):
    post_results = Entrez.read(Entrez.epost("protein", id=",".join(gis)))
    webenv = post_results["WebEnv"]
    query_key = post_results["QueryKey"]
    handle = Entrez.efetch(db="protein", rettype="gb",retmode="text", webenv=webenv, query_key=query_key)
    for s in SeqIO.parse(handle, "gb"):
        yield s

def taxonomy_from_gis(gis):
    """
    """
    for s in seq_from_gi(gis):
        print s.annotations["organism"]
        yield s.annotations["organism"]

def taxids_from_gis(gis):
    """
    """
    for s in seq_from_gi(gis):
        a=s.features[0].qualifiers['db_xref'][0]
        id=re.search(':(\d+)',a).group(1)
        print "Taxid ",id
        yield id

already_exists = []

def taxonomy_from_header(header, gi=None, species_re=None):
  if species_re is None:
    species_re = re.compile(r'\[(.*?)\]')
  match = species_re.findall(header)
  if match:
    organism = match[-1]
  elif gi:
    print "No match for {}: {}, searh NCBI".format(gi, header)
    for i in xrange(10):
      try:
        organism = taxonomy_from_gis([gi]).next()
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
      print header
      return Taxonomy.objects.get(name="unidentified")

    #Maybe the var is wrong
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

