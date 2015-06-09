from django.db.models import Max, Min
from browse.models import *
import browse
from djangophylocore.models import *
from django.db.models import Q

from django.shortcuts import redirect
from django.http import JsonResponse

from collections import Counter

import django_filters

def tax_sub_search(value):
    """
    """
    ids = set()

    if not format_query.current_query.startswith("taxonomy"):
        return list(ids)
    
    search_type = format_query.current_query[len("taxonomy"):]
    queryName = {"name{}".format(search_type):value.lower()}
    try:
        queryID = {"id{}".format(search_type):int(value)}
    except ValueError:
        queryID = {}
    query = Q(**queryName)|Q(**queryID)
    for taxon in Taxonomy.objects.filter(query):
        if taxon.type_name != "scientific name":
            #Convert homonym into real taxon
            taxon = taxon.get_scientific_names()[0]
        ids.add(taxon.id)
        children = set(taxon.children.values_list("id", flat=True).filter(rank=2))
        ids |= children
    format_query.current_query = "taxonomy__in"
    return list(ids)

#Fields that are allowed: each row contains:
#    POST name, django model paramter, POST name for search type, input type (must be in seach_types)
allowable_fields = [
    ("id_id", "id", "id_id_search_type", str),
    ("id_core_histone", "variant__core_type__id", "id_core_type_search_type", str),
    ("id_variant", "variant_id", "id_variant_search_type", str),
    ("id_gene", "gene", "id_gene_search_type", int),
    ("id_splice","splice", "id_splice_search_type", int),
    ("id_header", "header", "id_header_search_type", str),
    ("id_taxonomy", "taxonomy", "id_taxonomy_search_type", tax_sub_search),
    ("id_evalue", "evalue", "id_evalue_search_type", float),
    ("id_score", "score", "id_specificity_search_type", float),
    #("show_lower_scoring_models", None, None, (""))
]

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
search_types["text"] = search_types[str]
search_types["int"] = search_types[int]

class HistoneSearch(object):
    """
    """

    def __init__(self, request, parameters, reset=False, navbar=False):
        """
        """
        assert isinstance(parameters, dict)

        if reset:
            HistoneSearch.reset(request)
        elif "query" in request.session and request.session["query"]:
            parameters.update(request.session["query"])
            self.ajax = True

        self.errors = Counter()
        self.navbar = navbar
        self.request = request
        self.sanitized = False
        self.query_set = None
        self.redirect = None
        self.count = 0

        self.sanitize_parameters(parameters, reset)
        self.create_queryset(reset)

    @classmethod
    def reset(cls, request):
        """Clears search from session"""
        try:
            del request.session["query"]
            del request.session["sort"]
        except KeyError:
            pass

    @classmethod
    def all(cls, request):
        return cls(request, {}, reset=True)

    def sanitize_parameters(self, parameters, reset=False):
        """
        """
        self.request.session["query"] = {field:parameters[field] for fields in allowable_fields \
            for field in (fields[0], fields[2]) if field in parameters}

        if "search" in parameters:
            self.request.session["search"] = parameters["search"]
        
        sort_parameters = {"limit": 10, "offset":0, "sort":"evalue", "order":"asc"}
        sort_query = {p:parameters.get(p, v) for p, v in sort_parameters.iteritems()}

        if reset or not "sort" in self.request.session:
            self.request.session["sort"] = sort_query
        else:
            self.request.session["sort"].update(sort_query)

        #if not reset:
            #raise RuntimeError(str(self.request.session["query"]))

        self.sanitized = True

    def create_queryset(self, reset=False):
        if not self.sanitized:
            raise RuntimeError("Parameters must be sanitized")

        parameters = self.request.session["query"]

        if "search" in parameters:
            self.simple_search(search_text=parameters["search"])
            return
        
        query = format_query()

        added = []
        for form, field, search_type, convert in allowable_fields:
            value = None
            if form in parameters:
                added.append((field, search_type, parameters.get(form)))
                value = parameters.get(form)

            if not value: 
                continue

            search_type = parameters.get(search_type, "is")
            
            query.format(field, search_type, value, convert)                

        """specificity = parameters.get("specificity", 95):
        specificity_search_type = parameters.get("specificity_search_type", ">=")
        q, value = int_search_query("scores__specificty", specificity_search_type, specificity)
        query["scores__best"] = True
        query[q] = value
        """

        if query.has_errors():
            raise RuntimeError(str(query.errors))
            return False
        
        self.query_set = Sequence.objects.annotate(evalue=Min("scores__evalue"), score=Max("scores__score")).filter(**query)
        #if hasattr(self, 'ajax'):
        #    assert 0, list(self.query_set)
        self.count = self.query_set.count()
        #if not reset:
        #    raise RuntimeError(str(query))
        
    def sorted(self):
        #Sort by best score. Using e-value so we can compare HMMER and SAM
        result = self.query_set
        
        sort_by = self.request.session["sort"]["sort"]
        sort_order = "-" if self.request.session["sort"]["order"] == "desc" else ""
        sort_key = "{}{}".format(sort_order, sort_by)
        result = result.order_by(sort_key)
        assert 0, list(result)
        
        try:
            page_size = int(self.request.session["sort"]["limit"])
        except ValueError:
            page_size = 10

        try:
            page_number = int(self.request.session["sort"]["offset"])
        except ValueError:
            page_number = 0

        start = page_size*page_number
        end = start+page_size

        result = result[start:end]


        return result

    def all(self):
        return self.query_set.all()

    def __len__(self):
        return self.query_set.count()

    def count(self):
        return self.query_set.count()

    def get_dict(self):
        sequences = self.sorted()
        result = [{"gi":r.id, "variant":r.variant_id, "gene":r.gene, "splice":r.splice, "species":r.taxonomy.name, "score":r.score, "evalue":r.evalue, "header":r.header} for r in sequences]
        return {"total":self.count, "rows":result}

    def simple_search(self, search_text=None):
        """Search from simple text box in brand or sequence filter.
        """
        assert search_text is not None or "search" in self.search_query
        if search_text is None:
            search_text = self.search_query["search"]

        #Search all fields
        try:
            #If search is just a single digit, try to match id
            #Other int values don't make sense to search like this
            value = int(search_text)
            sequence = self.query_set.filter(id=value)
            if len(sequence):
                self.query_set = sequence
                request.session["query"]["id"] = value
                request.session["query"]["id_search_type"] = "is"
                return
        except ValueError:
            pass

        #search core_type, variant, old variant names, header if doens't match variant or core_type, taxonomy
        try:
            core_type = Histone.objects.get(id=search_text)
            sequences = self.query_set.filter(variant__core_type_id=core_type.id)
            request.session["query"]["core_type"] = core_type.id
            request.session["query"]["core_type_search_type"] = "is"
            if self.navbar:
                self.redirect = redirect("browse.views.browse_variants", core_type.id)
            return
        except Exception as e:
            pass
        
        try:
            variant = Variant.objects.get(id=search_text)
            sequences = self.query_set.filter(variant_id=variant.id)
            request.session["query"]["variant"] = variant.id
            request.session["query"]["variant_search_type"] = "is"
            if self.navbar:
                self.redirct = redirect("browse.views.browse_variant", variant.id)
            return
        except Exception as e:
            pass
        
        try:
            variant = OldStyleVariant.objects.get(id=search_text).updated_variant
            sequences = self.query_set.filter(variant_id=variant.id)
            request.session["query"]["variant"] = variant.id
            request.session["query"]["variant_search_type"] = "is"
            if self.navbar:
                self.redirct = redirect("browse.views.browse_variant", variant.id)
            return
        except:
            pass

        try:
            #Search species
            sequences = self.query_set.filter(taxonomy__name__icontains=search_text)
            if sequences.count() > 0:
                request.session["query"]["taxonomy"] = search_text
                request.session["query"]["taxonomy_search_type"] = "contains (case-insesitive)"
                return
        except:
            pass

        try:
            taxon = Taxonomy.objects.filter(name=parameters["search"])
            try:
                #Convert homonym into real taxon
                taxon = taxon.get_scientific_names()[0]
            except IndexError:
                #Already correct taxon
                pass
            taxa = taxon.children.filter(rank__name="species").values_list("pk")
            sequences = self.query_set.filter(taxonomy__in=taxa)
            if sequences.count() > 0:
                request.session["query"]["taxonomy"] = search_text
                request.session["query"]["taxonomy_search_type"] = "contains (case-insesitive)"
                return
        except:
            pass
            
        headers = self.query_set.filter(header__icontains=search_text)

        if headers.count() > 0:
            request.session["query"]["header"] = search_text
            request.session["query"]["header_search_type"] = "contains (case-insesitive)"
        else:
            #Search sequence moetifs
            sequences = self.query_set.filter(sequence__contains=search_text)
            self.query_set = headers
            self.redirect = None

class format_query(dict):
    current_query = None
    errors = Counter()
    def format(self, field, search_type, value, conv_type):
        #if allow and search_type not in allow:
        #    format_query.errors["Invalid search type ({}) for field {}".format(search_type, field)] += 1
         #   return False, errors
        try:
            if conv_type.__class__.__name__ == "function":
                search_type = search_types[str][search_type]
            else:
                search_type = search_types[conv_type][search_type]
        except KeyError:
            format_query.errors["Invalid search type ({}; {}) for field {}".format(search_type,conv_type,field)] += 1
            return False, self.errors


        format_query.current_query = "{}{}".format(field, search_type)

        try:
            if search_type.endswith("range"):
                if "-" in value:
                    value = map(conv_type, value.split("-"))
                else:
                    self.errors["Must include a dash if searching range"] += 1
            elif search_type.endswith("in"):
                if "," in value:
                    value = map(conv_type, value.split(","))
                else:
                    self.errors["Must include a dash if searching range"] += 1
            else:
                value = conv_type(value)

        except ValueError:
            format_query.errors["Must include correct character to seach {}".format(field)] += 1
        
        if len(format_query.errors) > 0:
            return

        self[format_query.current_query] = value
        self.current_query = None

    def has_errors(self):
        return len(format_query.errors) > 0





