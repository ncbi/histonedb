import re
from Bio import Entrez, SeqIO
# *Always* tell NCBI who you are
Entrez.email = "eli.draizen@nih.gov"

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

