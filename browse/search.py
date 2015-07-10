from django.db.models import Max, Min, Count
from browse.models import Histone, Variant, OldStyleVariant, Sequence, Score
import browse
from djangophylocore.models import Taxonomy
from django.db.models import Q

from django.shortcuts import redirect
from django.http import JsonResponse

from collections import Counter

import django_filters
from itertools import groupby

def tax_sub_search(value):
    """
    """
    ids = set()
    value = value.strip()
    
    if not format_query.current_query.startswith("taxonomy"):
        return list(ids)
    
    search_type = format_query.current_query[len("taxonomy"):]
    if search_type.endswith("iin"):
        search_type = "__iexact"
    elif search_type.endswith("in"):
        search_type = "__exact"
    queryName = {"name{}".format(search_type):value.lower()}
    try:
        queryID = {"id{}".format(search_type):int(value)}
    except ValueError:
        queryID = {}
    query = Q(**queryName)|Q(**queryID)
    
    taxons = []
    for taxon in Taxonomy.objects.filter(query):
        taxons.append(taxon)
        if taxon.type_name != "scientific name":
            #Convert homonym into real taxon
            taxon = taxon.get_scientific_names()[0]
        ids.add(taxon.id)
        children = set(taxon.children.values_list("id", flat=True).filter(rank=2))
        ids |= children

    format_query.current_query = "taxonomy__in"
    
    return list(ids)

def variant_sub_search(value):
    ids = set()

    if not format_query.current_query.startswith("variant__id"):
        return []

    try:
        #Search for current variant names or old variant names
        search_type = format_query.current_query[len("variant__id"):]
        queryCurrent = {"id{}".format(search_type):value}
        queryOld = {"old_names__name{}".format(search_type):value}
        query = Q(**queryCurrent)|Q(**queryOld)

        variant = Variant.objects.filter(query)
        if len(variant) > 0:
            format_query.current_query = "variant__id__in"
            return variant.values_list("id", flat=True)
        else:
            return variant.first().id
    except Variant.DoesNotExist:
        #Search new nomencalture with gene and splice isoform
        for histone in Histone.objects.all():
            for variant in histone.variants.all():
                if value.startswith(variant.id):
                    branch_points = value[len(variant.id):].split(".")
                    sequence_query = [("variant__id", variant.id)]
                    for branch_point in branch_points:
                        if branch_point.startswith("s"):
                            key = "splice"
                            branch_point = part[1:]
                        else:
                            key = "gene"
                        try:
                            branch_point = int(branch_point)
                            sequence_query.append((key, branch_point))
                        except ValueError:
                            pass
                    format_query.multiple = True
                    return sequence_query
        
  
    else:
        
    
        variants = Variant.objects.filter(query).values_list("id", flat=True)
    
        if len(variants) > 0:
            format_query.current_query = "variant__id__in"
            return list(variants)
        else:
            return []

#Fields that are allowed: each row contains:
#    POST name, django model paramter, POST name for search type, input type (must be in seach_types)
allowable_fields = [
    ("id_id", "id", "id_id_search_type", str),
    ("id_core_histone", "variant__core_type__id", "id_core_type_search_type", str),
    ("id_variant", "variant__id", "id_variant_search_type", variant_sub_search),
    ("id_gene", "gene", "id_gene_search_type", int),
    ("id_splice","splice", "id_splice_search_type", int),
    ("id_header", "header", "id_header_search_type", str),
    ("id_taxonomy", "taxonomy", "id_taxonomy_search_type", tax_sub_search),
    ("id_evalue", "evalue", "id_evalue_search_type", float),
    ("id_score", "score", "id_score_search_type", float),
    #("show_lower_scoring_models", None, None, (""))
]

search_types = {}
search_types[str] = {
    "is": "__exact",
    "is (case-insesitive)": "__iexact",
    "contains": "__contains",
    "contains (case-insesitive)": "__icontains",
    "starts with": "__startswith",
    "starts with (case-insesitive)": "__istartswith",
    "ends with": "__endswith",
    "ends with (case-insesitive)": "__iendswith",
    "in (comma separated)": "__in",
    "in (comma separated, case-insensitive)": "__iin"
}
search_types[int] = {
    ">": "__gt",
    ">=": "__gte",
    "is": "",
    "<": "__lt",
    "<=": "__lte",
    "range (dash separated)": "__range",
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

        self.errors = Counter()
        self.navbar = navbar
        self.request = request
        self.sanitized = False
        self.query_set = None
        self.redirect = None
        self.count = 0

        #assert 0, "search"+str("query" in request.session and request.session["query"])
        if reset:
            HistoneSearch.reset(request)
        elif "query" in request.session and request.session["query"]:
            parameters.update(request.session["query"])
            assert request.session["query"]=={u'id_variant': u'H2A.Z'}, parameters
            
            self.ajax = True

        if "search" in parameters:
            parameters.update(self.simple_search(parameters["search"]))

        if self.redirect:
            return

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

        if query.has_errors():
            raise RuntimeError(str(query.errors))
            return False
        
        #assert 0, query

        self.query_set = Sequence.objects.filter(
                (~Q(variant__id="Unknown") & Q(all_model_scores__used_for_classifiation=True)) | \
                (Q(variant__id="Unknown") & Q(all_model_scores__used_for_classifiation=False)) \
            ).annotate(
                num_scores=Count("all_model_scores"), 
                score=Max("all_model_scores__score"),
                evalue=Min("all_model_scores__evalue")
            ).filter(**query)
        #if hasattr(self, 'ajax'):
        
        self.count = self.query_set.count()
        #if not reset:
        #    raise RuntimeError(str(query))
        
    def sorted(self, unique=False):
        """Sort the query set and paginate results.
        """
        result = self.query_set
        
        sort_by = self.request.session["sort"]["sort"]
        sort_order = "-" if self.request.session["sort"]["order"] == "desc" else ""
        sort_key = "{}{}".format(sort_order, sort_by)
        result = result.order_by(sort_key)
        
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

        if unique:
            used = {}
            for seq in results:
                if not seq.sequence in used_taxa:
                    used[seq.sequence] = [seq.taxa.id]
                    yield seq
                elif seq.sequence in used_taxa and not seq.taxa.id in used[seq.sequence]:
                    used[seq.sequence].append(seq.taxa.id)
                    yield seq
                else:
                    pass
            """#Old method using group by, might be faster, but doesn;t preserver order
            result = sorted(self.query_set, key=lambda s:s.sequence)
            for sequence, same_sequences in groupby(result, key=lambda s:s.sequence):
                used_taxa = []
                for seq in same_sequences:
                    if seq.taxonomy.id not in used_taxa:
                        used_taxa.append(seq.taxonomy.id)
                        yield seq
            """
        else:
            for seq in result:
                yield seq

    def __len__(self):
        return self.query_set.count()

    def count(self):
        return self.query_set.count()

    def get_dict(self, unique=False):
        sequences = self.sorted(unique=unique)
        result = [{"id":r.id, "variant":r.variant_id, "gene":r.gene, "splice":r.splice, "taxonomy":r.taxonomy.name.capitalize(), "score":r.score, "evalue":r.evalue, "header":r.header[:80]} for r in sequences]
        return {"total":self.count, "rows":result}

    def simple_search(self, search_text):
        """Search from simple text box in brand or sequence filter.
        """
        parameters = {}
        #Search all fields
        try:
            #If search is just a single digit, try to match id
            #Other int values don't make sense to search like this
            value = int(search_text)
            sequence = self.query_set.filter(id=value)
            if len(sequence):
                parameters["id"] = value
                parameters["id_search_type"] = "is"
                return parameters
        except ValueError:
            pass

        #search core_type, variant, old variant names, header if doens't match variant or core_type, taxonomy
        try:
            core_type = Histone.objects.get(id=search_text)
            parameters["id_core_histone"] = core_type.id
            parameters["id_core_histone_search_type"] = "is"
            if self.navbar:
                self.redirect = redirect(core_type)
            return parameters
        except Histone.DoesNotExist:
            pass
        
        try:

            variant = Variant.objects.get(id=search_text)
            parameters["id_variant"] = variant.id
            parameters["id_variant_search_type"] = "is"
            if self.navbar:
                self.redirect = redirect(variant)
            return parameters
        except Variant.DoesNotExist:
            pass
        
        try:
            #Searches H2A.Z.1.s1
            sequences = Sequence.objects.filter(full_variant_name=search_text)
            if sequences.count() > 0:
                parameters["id_id"] = ",".join(sequences.values_list("id", flat=True))
                parameters["id_id_search_type"] = "in (comma separated, case-sensitive)"
                return parameters
        except:
            pass

        try:
            variant = OldStyleVariant.objects.get(name=search_text).updated_variant
            parameters["id_variant"] = variant.id
            parameters["id_variant_search_type"] = "is"
            if self.navbar:
                self.redirct = redirect(variant)
            return parameters
        except OldStyleVariant.DoesNotExist:
            pass

        try:
            #Search species
            sequences = Taxonomy.objects.filter(name__icontains=search_text)
            if sequences.count() > 0:
                parameters["id_taxonomy"] = search_text
                parameters["id_taxonomy_search_type"] = "contains (case-insesitive)"
                return parameters
        except Taxonomy.DoesNotExist:
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
                parameters["id_taxonomy"] = search_text
                parameters["id_taxonomy_search_type"] = "contains (case-insesitive)"
                return parameters
        except Taxonomy.DoesNotExist:
            pass
            
        headers = Sequence.objects.filter(header__icontains=search_text)

        if headers.count() > 0:
            parameters["id_header"] = search_text
            parameters["id_header_search_type"] = "contains (case-insesitive)"
        else:
            #Search sequence moetifs if everything else fails
            parameters["id_sequence"] = search_text
            parameters["id_sequence_search_type"] = "contains (case-insesitive)"

        return parameters

class format_query(dict):
    current_query = None
    multiple = False
    errors = Counter()
    def format(self, field, search_type, value, conv_type):
        #if allow and search_type not in allow:
        #    format_query.errors["Invalid search type ({}) for field {}".format(search_type, field)] += 1
         #   return False, errors
        search_type = search_type.replace("&gt;", ">").replace("&lt;", "<")
        
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
                values = []
                if "-" in value:
                    for v in value.split("-"):
                        values += conv_type(v)
                else:
                    self.errors["Must include a dash if searching 'range'"] += 1
                value = values
            elif search_type.endswith("in"):
                values = []
                if "," in value:
                    for v in value.split(","):
                        values += conv_type(v)
                else:
                    self.errors["Must include a comma if searching 'in'"] += 1
                value = values
            else:
                value = conv_type(value)

        except ValueError:
            format_query.errors["Must include correct character to seach {}".format(field)] += 1
        
        if len(format_query.errors) > 0:
            return

        if format_query.multiple and isinstance(value, list):
            for field, search in value:
                self[field] = search
        else:
            self[format_query.current_query] = value
        format_query.current_query = None
        format_query.multiple = False

    def has_errors(self):
        return len(format_query.errors) > 0





