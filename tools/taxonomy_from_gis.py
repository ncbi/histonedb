from Bio import Entrez, SeqIO
# *Always* tell NCBI who you are
Entrez.email = "eli.draizen@nih.gov"

def taxonomy_from_gis(gis):
    """
    """
    post_results = Entrez.read(Entrez.epost("protein", id=",".join(gis)))
    webenv = post_results["WebEnv"]
    query_key = post_results["QueryKey"]
    handle = Entrez.efetch(db="protein", rettype="gb",retmode="text", webenv=webenv, query_key=query_key)
    for s in SeqIO.parse(handle, "gb"):
        yield s.annotations["organism"]

species_re = re.compile(r'\[(.+)\]$')
def taxonomy_from_header(header, gi):
	match = species_re.find(header)
	if match:
		return match.groups(0)
	else:
		return taxonomy_from_gis([gi])

