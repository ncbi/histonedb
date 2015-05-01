from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
from django.shortcuts import get_list_or_404

from server.forms import SearchForm, FilterForm
from server.search import search

#Django libraires
from server.models import *
from phylocore_models import *

#BioPython
from Bio import SeqIO

def browse_types(request):
	"""Home"""
	return render(request, 'browse_types.html', {})

def browse_variants(request, histone_type):
	print "HI"
	data = {
		"histone_type": histone_type,
		"histone_description": "NOPE", #histone_description,
		"browse_section": "type",
		"name": histone_type,
		"filter_form": FilterForm(),
	}
	print data
	return render(request, 'browse_variants.html', data)

def browse_variant(request, histone_type):
	data = {
		"histone_type": histone_type,
		"histone_description": "NOPE", #histone_description,
	}
	return render(request, 'browse_varaint.html', {})

def search_view(request):
	if request.method == "POST":
		status, result = search(request.POST)
		if status:
			data = {'result':result}
		else:
			data = {'errors':result}
	else:
		data = {'search_form':SearchForm()}
	print render(request, 'search.html', data)
	return render(request, 'search.html', data)

def upload(request):
	if request.method == "POST":
		data = {'result':None}
	else:
		data = {}

	return render(request, 'upload.html', data)


def get_sequence_table_data(request, browse_type, seq_type):
	json = [{"GI": 1, "Variant":"H2A.Z", "Gene":1, "Splice":1, "Species":"Homo sapien", "header":"H2A.Z.1.s1 [Homo sapien]", "Score":600.2, "E-value":1e-100, "Program":"hmmer3.1b2"}]
	#json = [[1, "H2A.Z", 1, 1, "Homo sapien", "H2A.Z.1.s1 [Homo sapien]", 600.2, 1e-100, "hmmer3.1b2"]]
	return JsonResponse(json, safe=False)

	if browse_type == "type":
		sequences = Sequence.objects.filter(core_type=search)
	elif browse_type == "variant":
		sequences = Sequence.objects.filter(varaint=search)
	else:
		#404?
		return

	if seq_type == "all":
		#Downalod
		pass
	

def get_starburst_json(request, browse_type, search):
	"""
	"""
	if browse_type == "type":
		sequences = Sequence.objects.filter(core_type=search)
	elif browse_type == "variant":
		sequences = Sequence.objects.filter(variant=search)
	else:
		raise Http404("Must only search for core histone types ('type') or 'variants')")

	sunburst = []
	colors = {"eukaryota":"#6600CC", "prokaryota":"#00FF00", "archea":"#FF6600"}
	for sequence in sequences:
		taxa = sequence.taxonomy
		root = sunburst
		for i, curr_taxa in enumerate(reversed(taxa.parents.all())+taxa):
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





