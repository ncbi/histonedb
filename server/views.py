from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
from django.shortcuts import get_list_or_404

from server.forms import SearchForm, FilterForm

#Django libraires
from server.models import *
from phylocore_models import *

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

def search(request):
	print "search"
	if request.method == "POST":
		errors = Counter()
		if request.POST["id_id"]:
			try:
				result = Sequence.objects.get(id=request.POST["id_id"])
			except Sequence.DoesNotExist:
				errors["No sequence with GI: {}".format(request.POST["id_id"])] += 1
		else:
			result = Sequence.objects

			if request.POST["id_variant"]:
				result.filter(variant__search=request.POST["id_variant"])

			if request.POST["id_taxonomy"]:
				try:
					#Better way to search?
					taxonomy = Taxonomy.objects.get(name=request.POST["id_taxonomy"])
					result.filter(taxonomy=taxonomy)
				except Sequence.DoesNotExist:
					errors["No taxa with name: {}".format(request.POST["id_taxonomy"])] + 1

			if request.POST["id_gene"]:
				#Change to >, >=, <, <=, ==, range
				try:
					result.get(id=request.POST["id_gene"])
				except Sequence.DoesNotExist:
					errors["No gene with index: {}".format(request.POST["id_gene"])] += 1

			if request.POST["id_splice"]:
				#Change to >, >=, <, <=, ==, range
				try:
					result.get(id=request.POST["id_splice"])
				except Sequence.DoesNotExist:
					errors["No splice isoform with index: {}".format(request.POST["id_splice"])] += 1

			if request.POST["header"]:
				result.filter(header__search=request.POST["id_header"])


		if len(errors) == 0:
			data = {'result':None}
		else:
			data = {'errors':errors}
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
	#json = '[{"GI": 1, "Variant":"H2A.Z", "Gene":1, "Splice":1, "Species":"Homo sapien", "header":"H2A.Z.1.s1 [Homo sapien]", "Score":600.2, "E-value":1e-100, "Program":"hmmer3.1b2"}'
	json = '[[1, "H2A.Z", 1, 1, "Homo sapien", "H2A.Z.1.s1 [Homo sapien]", 600.2, 1e-100, "hmmer3.1b2"]]'
	return HttpResponse(json, content_type="application/json")

def get_starburst_json(requst, browse_type, search):
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




