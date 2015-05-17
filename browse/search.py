from django.db.models import Max, Min
from browse.models import *
import browse
from djangophylocore.models import *

from django.http import HttpResponseRedirect
from django.http import JsonResponse

from collections import Counter

current_search = None

def search(parameters, initial_query=None):
    global current_search
    if current_search and current_search.parameters == parameters:
        return current_search
    else:
        current_search = HistoneSearch(parameters, initial_query)

class HistoneSearch(object):
    current_search = None
    def __init__(self, search_query, sort_query, navbar=False):
        self.search_query = search_query
        self.initial_query = (search_query, sort_query)
        self.sort_query = sort_query
        self.errors = Counter()
        self.navbar = navbar
        self.query_set = Sequence.objects.annotate(evalue=Min("scores__evalue")).filter()
        self.update_queryset(search_query)
        self.redirect = None
        HistoneSearch.current_search = self

    @classmethod
    def search(cls, parameters, navbar=False):
        search_parameters = [
            "id", "id_search_type", 
            "core_type", "core_type_search_type", 
            "variant", "variant_search_type", 
            "gene", "gene_search_type", 
            "splice", "splice_search_type", 
            "header", "header_search_type", 
            "taxonomy", "taxonomy_search_type",
            #"specificity", "specificity_search_type"
            #"evalue", "evalue_search_type"
            #"score", "specificity_search_type"
            #"scoring_program"
            #"show_lower_scoring_models"
            ]
        if "search" not in parameters:
            search_query = {p:v for p, v in parameters.iter_items() if p in search_parameters}
        else:
            search_query = {"search":parameters["search"]}

        sort_parameters = {"limit": 10, "offset":0, "sort":"evalue", "order":"asc"}
        sort_query = {p:parameters.get(p, v) for p, v in sort_parameters.iteritems()}

        if not HistoneSearch.current_search:
            print "new search"
            return cls(search_query, sort_query, navbar=navbar)
        elif parameters.get("reset", False):
            return cls(*HistoneSearch.initial_query)
        else:
            if search_query != HistoneSearch.current_search.search_query:
                #Add filters to current
                HistoneSearch.current_search.update_queryset(search_query)
            if not HistoneSearch.current_search.sort_query.get("limit", 10) == sort_query["limit"]:
                HistoneSearch.current_search.sort_query["limit"] = sort_query["limit"]
            if not HistoneSearch.current_search.sort_query.get("offset", 0) == sort_query["offset"]:
                HistoneSearch.current_search.sort_query["offset"] = sort_query["offset"]
            if not HistoneSearch.current_search.sort_query.get("sort", "evalue") == sort_query["evalue"]:
                HistoneSearch.current_search.sort_query["sort"] = sort_query["order"]
            if not HistoneSearch.current_search.sort_query.get("order", "asc") == sort_query["order"]:
                HistoneSearch.current_search.sort_query["order"] = sort_query["order"]

            HistoneSearch.current_search.navbar = navbar
            
            return HistoneSearch.current_search

    def update_queryset(self, parameters):
        if "search" in parameters:
            result = self.simple_search()
            if self.navbar and type(result) == HttpResponseRedirect:
                self.redirect = result
                return
        else:
            query = format_query()  
              
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
                return False

            self.query_set.filter(**query)

    def sorted(self):
        #Sort by best score. Using e-value so we can compare HMMER and SAM
        result = self.query_set.clone()
        sort_by = self.sort_query["sort"]
        sort_order = "-" if self.sort_query["order"] == "desc" else ""
        sort_key = "{}{}".format(sort_order, sort_by)
        result = result.sort_by(sort_key)

        try:
            page_size = int(parameters.get("limit", 10))
        except ValueError:
            page_size = 10

        try:
            page_number = int(parameters.get("offset", 0))
        except ValueError:
            page_number = 10

        start = page_size*(page_number+1)
        end = start+page_size+2

        result = result[start:end]

        return result

    def all(self):
        return self.query_set.all()

    def __len__(self):
        return self.query_set.count()

    def count(self):
        return self.query_set.count()

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
                return True
        except ValueError:
            pass

        #search core_type, variant, old variant names, header if doens't match variant or core_type, taxonomy
        try:                
            core_type = Histone.objects.get(id=search_text)
            if self.navbar:
                # Go to histone browse page
                return HttpResponseRedirect("/type/{}/".format(core_type.id))
            else:
                sequences = self.query_set.filter(variant__core_type_id=core_type.id)
                self.query_set = sequence
                return True
        except:
            pass
        
        try:
            variant = Variants.objects.get(id=search_text)
            if self.navbar:
                # Go to variant browse page
                return HttpResponseRedirect("/variant/{}/".format(variant.id))
            else:
                sequences = self.query_set.filter(variant_id=variant.id)
                self.query_set = sequence
                return True
        except:
            pass
        
        try:
            variant = OldStyleVariants.objects.get(id=search_text).updated_variant
            if self.navbar:
                # Go to vaiant browse page
                return HttpResponseRedirect("/variant/{}/".format(variant.id))
            else:
                sequences = self.query_set.filter(variant_id=variant.id)
                self.query_set = sequences
                return True
        except:
            pass

        try:
            """try:
                taxon = Taxonomy.objects.filter(name=parameters["search"])
            except:

            try:
                #Convert homonym into real taxon
                taxon = taxon.get_scientific_names()[0]
            except IndexError:
                #Already correct taxon
                pass"""
            sequences = self.query_set.filter(taxonomy__name__contains=search_text)
            if sequences.count() > 0:
                self.query_set = sequence
                return True
        except:
            pass
            
        headers = self.query_set.filter(header__contains=search_text)

        if headers.count() > 0:
            self.query_set = headers
            return True
        else:
            #Search seuence moetifs
            sequences = self.query_set.filter(sequence__contains=search_text)
            self.query_set = headers
            return True

    def get_dict(self):
        sequences = self.sorted()
        result = [{"gi":r.id, "variant":r.variant_id, "gene":r.gene, "splice":r.splice, "species":r.taxonomy.name, "evalue":r.evalue, "header":r.header} for r in sequences]
        return {"count":len(result), "rows":result}
        
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
    def format(self, field, search_type, value, conv_type, allow=None):
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




