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

def browse_types(request):
	"""Home"""
	return render(request, 'browse_types.html', {})

def browse_variants(request, histone_type):
	data = {
		"histone_type": histone_type,
		"histone_description": "NOPE", #histone_description,
		"browse_section": "type",
		"name": histone_type,
		"filter_form": FilterForm(),
	}
	print data
        print "WHY ISN't it loading?"
	return render(request, 'browse_variants.html', data)

def browse_variant(request, histone_type):
	data = {
		"histone_type": histone_type,
		"histone_description": "NOPE", #histone_description,
	}
	return render(request, 'browse_varaint.html', {})

def search(request):
	data = {}
	if request.method == "POST":
		print "POST", request.POST
		HistoneSearch.current_search = None
		result = HistoneSearch.search(request.POST, navbar=True)

		if type(result) == HttpResponseRedirect: 
			return result

		if len(result.errors) > 0:
			HistoneSearch.current_search = None
			data['errors'] = result.errors
	else:
		#All sequences
		HistoneSearch.current_search = None
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
		parameters = request.GET
	else:
		return "false"

	if browse_type == "type":
		parameters["core_type"]=search
	elif browse_type == "variant":
		parameters["variant"]=search

	#Continues to filter previous search, unless paramters contains key 'reset'
	results = HistoneSearch.search(parameters)

	if len(results.errors) > 0:
		return "false"
	
	return JsonResponse(results.get_dict())

def get_starburst_json(request, browse_type, search, debug=True):
	"""
	"""
	if debug:
		taxas = Taxonomy.objects.filter(name="root", type_name="scientific name")
		print len(taxas)
	elif browse_type == "type":
		taxas = Sequence.objects.filter(variant__core_type=search).values_list("taxonomy", flat=True).distinct()
	elif browse_type == "variant":
		taxas = Sequence.objects.filter(variant=search).values_list("taxonomy", flat=True).distinct()
	else:
		raise Http404("Must only search for core histone types 'type' or 'variants')")

	sunburst = []
	colors = {"eukaryota":"#6600CC", "prokaryota":"#00FF00", "archea":"#FF6600"}
	for taxa in taxas:
		print type(taxa)
		print taxa.children
		if debug:
			path = [taxa] + taxa.children.all()
			print path
		else:
			path = reversed(taxa.parents.all())+[taxa]
		root = sunburst
		print path
		return
		for i, curr_taxa in enumerate(path):
			print curr_taxa
			continue
			for node in root:
				if node["name"] == curr_taxa.name:
					break
			else:
				node = {"name":curr_taxa, "children":[]}
				root.append()
			root = node["children"]
		try:
			root["size"] += 1
		except KeyError:
			root["size"] = 1
	for node in sunburst[0]["children"]:
		try:
			node["color"] = colors[node["name"]]
		except KeyError:
			node["color"] = "#000000"

	if len(sunburst) == 0:
		raise Http404("No species {} for {}".format(search, browse_type))

	return JsonResponse(sunburst)

	root = Taxonomy.objects.get(name="root", type_name="scientific name")
	sunburst = create_sunburst(root)
	

def create_sunburst(root, sunburst=None, level=0):
	if sunburst is None:
		sunburst = {"name":root.name, "children":[]}

	for curr_taxa in root.direct_children.filter(type_name="scientific name").all():
		#print "{}{}".format(" "*level, curr_taxa.name),
		child = {"name":curr_taxa.name, "children":[]}
		#print child
		sunburst["children"].append(create_sunburst(curr_taxa, child, level=level+1))
	return sunburst





