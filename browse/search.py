from django.db.models import Max
from browse.models import *
from djangophylocore.models import *

from django.http import HttpResponseRedirect
from django.http import JsonResponse

from collections import Counter

_digits = re.compile('\d')

search_types = {}
search_types[str] = {
    "is": "",
    "is (case-insesitive)": "__iexact",
    "contains": "__contains",
    "contains (case-insesitive)": "__icontains",
    "starts with": "__startswith",
    "starts with (case-insesitive)": "__istartswith",
    "ends with": "__endswith",
    "ends with (case-insesitive)": "__iendswith",
    "in (comma separated, case-sensitive)": "__in"
}
search_types[int] = {
    ">": "__gt",
    ">=": "__gte",
    "is": "",
    "<": "__lt",
    "<=": "__lte",
    "range (dash separator)": "__range",
    "in (comma separated)": "__in",
}
search_types[float] = search_types[int]

class format_query(dict):
    current_query = None
    errors = Counter()
    def format(field, search_type, value, conv_type, allow=None):
        if allow and search_type not in allow:
            errors["Invalid search type ({}) for field {}".format(search_type, field)] += 1
            return False, errors

        try:
            search_type = search_types[conv_type][search_type]
        except KeyError:
            self.errors["Invalid search type ({}) for field {}".format(search_type, field)] += 1
            return False, errors

        self.current_query = "{}{}".format(field, search_type)

        try:
            if search_type == "range":
                if "-" in value:
                    value = map(conv_type, value.split("-"))
                else:
                    self.errors["Must include a dash if searching range"] += 1
            elif search_type == "in":
                if "," in value:
                    value = map(conv_type, value.split(","))
                else:
                    self.errors["Must include a dash if searching range"] += 1
            else:
                value = conv_type(value)
        except:
            self.errors["Must include correct character to seach {}".format(field)] += 1
        
        if len(errors) > 0:
            return

        self[self.current_query] = value
        self.current_query = None

    def has_errors(self):
        return len(errors) > 0

def tax_sub_search(value):
    """
    """
    ids = set()
    for taxon in Taxonomy.objects.filter(**{format_query.current_query: value}):
        try:
            #Convert homonym into real taxon
            taxon = taxon.get_scientific_names()[0]
        except IndexError:
            #Already correct taxon
            pass
        ids |= set(taxon.children.values_list("id", flat=True).filter(rank=2))
    format_query.current_query = "taxonomy__in"
    return ids

def queryset_to_json(sequences):
    return JsonResponse({"count":len(sequences), "rows":sequences.all()})

def simple_search(parameters, source, query=None):
    """Search from simple text box in brand or sequence filter.
    """
    #Search all fields
    try:
        #If search is just a single digit, try to match id
        #Other int values don't make sense to search like this
        value = int(parameters["search"])
        sequence = Sequences.objects.filter(id=value)
        if len(sequence):
            return queryset_to_json(sequence)
    except ValueError:
        pass

    #search core_type, variant, old variant names, header if doens't match variant or core_type, taxonomy
    try:                
        core_type = Histone.objects.get(id=parameters["search"])
        if source == "navbar":
            # Go to histone browse page
            return HttpResponseRedirect("/type/{}/".format(core_type.id))
        else:
            return queryset_to_json(Sequences.objects.filter(variant__core_type=core_type))
    except Histone.DoesNotExist:
        pass
    
    try:
        variant = Variants.objects.get(id=parameters["search"])
        if source == "navbar":
            # Go to variant browse page
            return HttpResponseRedirect("/variant/{}/".format(variant.id))
        else:
            return queryset_to_json(Sequences.objects.filter(variant=variant))
    except Variant.DoesNotExist:
        pass
    
    try:
        variant = OldStyleVariants.objects.get(id=parameters["search"]).updated_variant
        if source == "navbar":
            # Go to vaiant browse page
            return HttpResponseRedirect("/variant/{}/".format(variant.id))
        else:
            return queryset_to_json(Sequences.objects.filter(variant=variant))
    except OldStyleVariant.DoesNotExist:
        pass

    try:
        taxon = Taxonomy.objects.filter(name=parameters["search"])
        try:
            #Convert homonym into real taxon
            taxon = taxon.get_scientific_names()[0]
        except IndexError:
            #Already correct taxon
            pass
    except Taxonomy.DoesNotExist:
        pass
        
    return queryset_to_json(Sequences.objects.filter(header__contains=parameters["search"]))

def search(parameters, query=None):
    errors = Counter()

    if query is None:
        query = []
              
    fields = [
        ("id", "id_search_type", int, ("is", "in")),
        ("core_type", "core_type_search_type", str, None),
        ("variant", "variant_search_type", str, None),
        ("gene", "gene_search_type", int, None),
        ("splice", "splice_search_type", int, None),
        ("header", "header_search_type", str, None),
        ("taxonomy", "taxonomy_search_type", tax_sub_search, None)
        #("specificity", "specificity_search_type", int, None),
        #("evalue", "evalue_search_type", float, None),
        #("score", "specificity_search_type", float, None),
        #("scoring_program", None, str, ("is")),
        #("show_lower_scoring_models", None, None, (""))
    ]

    query = format_query()
    for field, search_type, convert, allow in fields:
        if not parameters.get(field): continue
        search_type = parameters.get(search_type, "is")
        query.format(field, search_type, parameters[field], convert, allow)
            

    """specificity = parameters.get("specificity", 95):
    specificity_search_type = parameters.get("specificity_search_type", ">=")
    q, value = int_search_query("scores__specificty", specificity_search_type, specificity)
    query["scores__best"] = True
    query[q] = value
    """

    if query.has_errors():
        return False, query.errors

    if len(query) > 0:
        result = Sequences.objects.filter(**query)
    else:
        result = Sequences.objects.all()

    #Sort by best score. Using e-value so we can compare HMMER and SAM
    result = result.aggregate(evalue=Min("scores__evalue")).sort_by("-evalue")

    return True, result
