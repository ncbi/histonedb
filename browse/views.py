from django.shortcuts import render
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.http import JsonResponse
from django.shortcuts import get_list_or_404

from browse.forms import SearchForm, FilterForm
from browse.search import HistoneSearch

#Django libraires
from browse.models import *
from djangophylocore.models import *

#BioPython
from Bio import SeqIO

from django.db.models import Min

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

	HistoneSearch.reset()

	data = {
		"histone_type": histone_type,
		"histone_description": core_histone.description,
		"browse_section": "type",
		"name": histone_type,
		"tree_url": "browse/trees/{}.xml".format(core_histone.id),
		"seed_file":"browse/seeds/{}.fasta".format(core_histone.id),
		"filter_form": FilterForm(),
	}

	return render(request, 'browse_variants.html', data)

def browse_variant(request, variant):
	try:
		variant = Variant.objects.get(id=variant)
	except:
		return "404"

	HistoneSearch.reset()

	data = {
		"core_type": variant.core_type.id,
		"variant": variant.id,
		"name": variant.id
		"seed_file":"browse/seeds/{}/{}.fasta".format(variant.core_type.id, variant.id),
		"browse_section": "variant",
		"description": variant.description,
	}
	return render(request, 'browse_variant.html', data)

def search(request):
	data = {}
	if request.method == "POST":
		result = HistoneSearch(request, request.POST, reset=True, navbar=True)
		
		if type(result) == type(redirect): 
			return result.redirect

		if len(result.errors) > 0:
			data['errors'] = result.errors
	else:
		
	return render(request, 'search.html', data)

def upload(request):
	if request.method == "POST":
		data = {'result':None}
	else:
		data = {}

	return render(request, 'upload.html', data)

def get_sequence_table_data(request, browse_type, search):
	"""result = Sequence.objects.filter(variant__core_type="H2A", taxonomy__name="homo sapiens").annotate(evalue=Min("scores__evalue")).order_by("evalue")[:10]
	result = [{"gi":r.id, "variant":r.variant_id, "gene":r.gene, "splice":r.splice, "species":r.taxonomy.name, "evalue":r.evalue, "header":r.header} for r in result]
	return JsonResponse({"count": len(result), "rows":result})

	json = {"count":1, "rows":[{"gi": 1, "variant":"H2A.Z", "gene":1, "splice":1, "species":"Homo sapien", "header":"H2A.Z.1.s1 [Homo sapien]", "evalue":1e-100, "program_train":"hmmer3.1b2", "program_test":"hmmer3.1b2"}]}
	#json = [[1, "H2A.Z", 1, 1, "Homo sapien", "H2A.Z.1.s1 [Homo sapien]", 600.2, 1e-100, "hmmer3.1b2"]]
	return JsonResponse(json)"""

	if request.method == "GET":
		parameters =  request.GET.dict()
	else:
		#Returning 'false' stops Bootstrap table
		return "false"

	if browse_type == "type":
		parameters["core_type"] = search
		parameters["core_type_search_type"] = "is"
	elif browse_type == "variant":
		parameters["variant"] = search
		parameters["variant_search_type"] = "is"

	#Continues to filter previous search, unless paramters contains key 'reset'
	results = HistoneSearch(request, parameters)

	if len(results.errors) > 0:
		#Returning 'false' stops Bootstrap table
		return "false"
	
	return JsonResponse(results.get_dict())

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
