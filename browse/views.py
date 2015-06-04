from django.shortcuts import render
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.http import JsonResponse
from django.shortcuts import redirect
from django.shortcuts import get_list_or_404

from browse.forms import SearchForm, FilterForm, UploadFileForm
from browse.search import HistoneSearch
from browse.process_upload import process_upload

#Django libraires
from browse.models import *
from djangophylocore.models import *

#BioPython
from Bio import SeqIO

from django.db.models import Min, Count

#Set2 Brewer, used in variant colors
colors = [
	"#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
    "#a6d854",
    "#ffd92f",
    "#e5c494"
    ]

def treemap(request):
	return render(request, 'circle.html', {"starburst_url": "data/type/json/{}/species/".format("H2A")})

def help(request):
	return render(request, 'help.html', {})

def browse_types(request):
	"""Home"""
	return render(request, 'browse_types.html', {})

def browse_variants(request, histone_type):
	try:
		core_histone = Histone.objects.get(id=histone_type)
	except:
		return "404"

	#Store sequences in session, accesed in get_sequence_table_data
	HistoneSearch(request, {"core_type":histone_type}, reset=True)

	variants = core_histone.variants.annotate(num_sequences=Count('sequences')).all().values_list("id", "num_sequences")
	variants = [(id, num, color) for (id, num), color in zip(variants, colors)]


	data = {
		"histone_type": histone_type,
		"histone_description": core_histone.description,
		"browse_section": "type",
		"name": histone_type,
		"variants": variants,
		"tree_url": "browse/trees/{}.xml".format(core_histone.id),
		"seed_file":"browse/seeds/{}.fasta".format(core_histone.id),
		"filter_form": FilterForm(),
	}

	return render(request, 'browse_variants.html', data)

def browse_variant(request, histone_type, variant):
	try:
		variant = Variant.objects.get(id=variant)
	except:
		return "404"

	HistoneSearch(request, {"variant":variant.id}, reset=True)

	data = {
		"core_type": variant.core_type.id,
		"variant": variant.id,
		"name": variant.id,
		"sunburst_url": "browse/sunbursts/{}/{}.json".format(variant.core_type.id, variant.id),
		"seed_file":"browse/seeds/{}/{}.fasta".format(variant.core_type.id, variant.id),
		"browse_section": "variant",
		"description": variant.description,
	}
	return render(request, 'browse_variant.html', data)

def search(request):
	data = {}
	if request.method == "POST":
		result = HistoneSearch(request, request.POST, reset=True, navbar=True)
	else:
		#Show all sequence
		result = HistoneSearch.all(request)
		
	if result.redirect: 
		return result.redirect
		
	return render(request, 'search.html', {"result":result})

def upload(request):
	if request.method == "POST":
		type = request.POST.get("type")
		if request.POST.get("sequences"):
			format = "text"
			sequences = request.POST["sequences"]
		elif request.POST.get("file"):
			format="file"
			sequences = request.POST["file"]
		data = process_upload(type, sequences, format)
	else:
		data = {'form':UploadFileForm()}

	return render(request, 'upload.html', data)

def get_sequence_table_data(request):
	"""Downloads the previos search and converts into json required by Bootstrap table
	"""

	if request.method == "GET":
		parameters =  request.GET.dict()
	else:
		#Returning 'false' stops Bootstrap table
		return "false"

	#Continues to filter previous search, unless paramters contains key 'reset'
	results = HistoneSearch(request, parameters)

	if len(results.errors) > 0:
		#Returning 'false' stops Bootstrap table
		return "false"
	
	return JsonResponse(results.get_dict())

def get_all_scores(request, ids=None):
	if ids is None and request.method == "GET" and "id" in request.GET:
		ids = request.GET["id"]
	else:
		#Returning 'false' stops Bootstrap table
		return "false"
	variants = Variant.objects.all().values_list("id", flat=True)
	rows = [["none"]*len(ids) for _ in xrange(len(variants))]
	for i, id in enuemrate(ids):
		try:
			sequence = Sequence.objects.get(id=id)
		except:
			return "404"
		scores = sequence.scores.values("variant").annotate(highest=Max("scores__score"))
		for variant in enumerate(scores):
			rows[i][variants.index(score.variant.id)] = score.score

	return JsonResponse({"total":len(rows), "rows":rows})


def get_starburst_json(request, browse_type, search, debug=False):
	"""
	"""

	

def create_sunburst(root, sunburst=None, level=0):
	if sunburst is None:
		sunburst = {"name":root.name, "children":[]}

	for curr_taxa in root.direct_children.filter(type_name="scientific name").all():
		#print "{}{}".format(" "*level, curr_taxa.name),
		child = {"name":curr_taxa.name, "children":[]}
		#print child
		sunburst["children"].append(create_sunburst(curr_taxa, child, level=level+1))
	return sunburst

def get_phyloxml_of_type(request):
	return None
