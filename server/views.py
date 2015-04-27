from django.shortcuts import render
from django.http import HttpResponse
import json

from server.forms import SearchForm, FilterForm
# Create your views here.

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
	return render(request, 'browse_variants.html', data)

def browse_variant(request, histone_type):
	data = {
		"histone_type": histone_type,
		"histone_description": "NOPE", #histone_description,
	}
	return render(request, 'browse_varaint.html', {})

def search(request):
	if request.method == "POST":
		data = {'result':None}
	else:
		data = {'search_form':SearchForm()}

	return render(request, 'search.html', data)

def upload(request):
	if request.method == "POST":
		data = {'result':None}
	else:
		data = {}

	return render(request, 'upload.html', data)


def get_sequence_table_data(request, browse_type, seq_type):
	if seq_type == "all":
		#Downalod
		pass
	#json = '[{"GI": 1, "Variant":"H2A.Z", "Gene":1, "Splice":1, "Species":"Homo sapien", "header":"H2A.Z.1.s1 [Homo sapien]", "Score":600.2, "E-value":1e-100, "Program":"hmmer3.1b2"}'
	json = '[[1, "H2A.Z", 1, 1, "Homo sapien", "H2A.Z.1.s1 [Homo sapien]", 600.2, 1e-100, "hmmer3.1b2"]]'
	return HttpResponse(json, content_type="application/json")