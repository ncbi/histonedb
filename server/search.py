from server.models import *
from server.phylocore_models import *

def search(parameters, preset_query=None):
	errors = Counter()

	if preset_query is None:
		query = preset_query
	else:
		query = []

	if parameters["id"]:
		try:
			query["id"] = parameters["id"]
		except Sequence.DoesNotExist:
			errors["No sequence with GI: {}".format(parameters["id"])] += 1

	if parameters["variant"]:
		query["variant__search"] = parameters["variant"]

	if parameters["taxonomy"]:
		try:
			species = Taxonomy.objects.get(name=parameters["taxonomy"]).children.filter(rank="species")
			query["taxonomy__in"] = species
		except Sequence.DoesNotExist:
			errors["No taxa with name: {}".format(parameters["taxonomy"])] + 1

	if parameters["gene"]:
		search_type = parameters["gene_search_type"]
		try:
			if search_type == ">":
				query["gene__gt"] = int(parameters["gene"])
			elif search_type == ">=":
				query["gene__gte"] = int(parameters["gene"])
			elif search_type == "==":
				query["gene"] = int(parameters["gene"])
			elif search_type == "<":
				query["gene__lt"] = int(parameters["gene"])
			elif search_type == "<=":
				query["gene__lte"] = int(parameters["gene"])
			elif search_type == "range":
				if "-" in parameters["gene"]
					query["gene__range"] = map(int, parameters["gene"].split("-"))
				else:
					errors["Must include a dash if searching range"] += 1
			else:
				raise Http404()
		except ValueError:
			errors["Must include a number to seach gene"] += 1

	if parameters["splice"]:
		search_type = parameters["splice_search_type"]
		try:
			if search_type == ">":
				query["splice__gt"] = int(parameters["splice"])
			elif search_type == ">=":
				query["splice__gte"] = int(parameters["splice"])
			elif search_type == "==":
				query["splice"] = int(parameters["splice"])
			elif search_type == "<":
				query["splice__lt"] = int(parameters["splice"])
			elif search_type == "<=":
				query["splice__lte"] = int(parameters["splice"])
			elif search_type == "range":
				if "-" in parameters["splice"]
					query["splice__range"] = map(int, parameters["splice"].split("-"))
				else:
					errors["Must include a dash if searching range"] += 1
			else:
				raise Http404()
		except ValueError:
			errors["Must include a number to seach splice"] += 1

	if parameters["header"]:
		query["header__search"] = parameters["id_header"]

	if len(errors) > 0:
		return False, errors

	result = Sequences.filter(**query) 
	return True, result
