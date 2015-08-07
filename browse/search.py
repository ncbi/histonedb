from django.db.models import Max, Min, Count, Sum
from browse.models import Histone, Variant, OldStyleVariant, Sequence, Score
import browse
from djangophylocore.models import Taxonomy
from django.db.models import Q

from django.shortcuts import redirect
from django.http import JsonResponse

from collections import Counter

import django_filters
from itertools import groupby
import re

####These subroutines are realization of searching in histone database


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

        variant = Variant.objects.filter(query).distinct()
        if len(variant) > 0:
            if len(variant) > 1:
                format_query.current_query = "variant__id__in"
                return variant.values_list("id", flat=True)
            else:
                return variant.first().id
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

bool_conv = lambda b: b in ["true", "on"]
#Fields that are allowed: each row contains:
#    POST name, django model paramter, POST name for search type, input type (must be in seach_types)
allowable_fields = {
    "id_id": ("id", "id_id_search_type", str),
    "id_hist_type": ("variant__hist_type__id", "id_core_type_search_type", str),
    "id_variant": ("variant__id", "id_variant_search_type", variant_sub_search),
    "id_gene": ("gene", "id_gene_search_type", int),
    "id_splice": ("splice", "id_splice_search_type", int),
    "id_header": ("header", "id_header_search_type", str),
    "id_taxonomy": ("taxonomy", "id_taxonomy_search_type", tax_sub_search),
    "id_evalue": ("evalue", "id_evalue_search_type", float),
    "id_score": ("score", "id_score_search_type", float),
    "id_sequence": ("sequence", "id_sequence_search_type", str),
    "id_reviewed": ("reviewed", "", bool_conv),
}

bool_fields = {
    "id_refseq":{"id_header":"|ref|", "id_header_search_type":"contains"},
}

search_types = {}
search_types[str] = {
    "is": "__exact",
    "is (case-insesitive)": "__iexact",
    "contains": "__contains",
    "contains (case-insensitive)": "__icontains",
    "starts with": "__startswith",
    "starts with (case-insensitive)": "__istartswith",
    "ends with": "__endswith",
    "ends with (case-insensitive)": "__iendswith",
    "in (comma separated)": "__in",
    "regex": "__regex",
    "regex (case-insensitive)": "__iregex"
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

import itertools
class Indexable(object):
    """Allow an iterable to be sliceable. We want this to be able slice a query
    if it is 'unique,' which is custom genrator mimicking query set"""
    def __init__(self,it):
        self.it = iter(it)
    def __iter__(self):
        for elt in self.it:
            yield elt
    def __getitem__(self,index):
        try:
            return next(itertools.islice(self.it,index,index+1))
        except TypeError:
            return list(itertools.islice(self.it,index.start,index.stop,index.step))

class HistoneSearch(object):
    """
    This class constructs query_sets - to search any histones in DB.
    The search parameters are in parameters - the dict should follow that defined in
    allowable_fields
    get_dict is used to evaluate the query_set and return the dict for rendering.
    """

    def __init__(self, request, parameters, navbar=False):
        """
        """
        assert isinstance(parameters, dict)

        self.errors = Counter()
        self.navbar = navbar
        self.request = request
        self.sanitized = False
        self.redirect = None
        self.query = format_query()
        self.count = 0
        self.sort = {}
        self.parameters = parameters
        self.unique = self.parameters.get("id_unique", False) in ["on", "true"]
        #this is the initial query set that is refined later by create_queryset

        # self.query_set = Sequence.objects.filter(all_model_scores__used_for_classifiation=True \
        #         ).annotate(
        #         num_scores=Count("all_model_scores"),
        #         score=Max("all_model_scores__score"),
        #         evalue=Min("all_model_scores__evalue")
        #     )

        self.query_set = Sequence.objects.all().annotate(\
                num_scores=Count("all_model_scores"),
                score=Max("all_model_scores__score"),
                evalue=Min("all_model_scores__evalue")
            )

        if "search" in self.parameters:
            self.simple_search(self.parameters["search"])

        if self.redirect:
            return

        self.update_bool_options()
        self.create_queryset()
        self.get_sort_options()

    @classmethod
    def all(cls, request):
        return cls(request, {})

    def get_sort_options(self):
        """
        """
        sort_parameters = {"limit": 10, "offset":0, "sort":"evalue", "order":"asc"}
        sort_query = {p:self.parameters.get(p, v) for p, v in sort_parameters.iteritems()}
        self.sort = sort_query

    def update_bool_options(self):
        for bool_field, updated_options in bool_fields.iteritems():
            if self.parameters.get(bool_field, False) in ["true", "on"]:
                self.parameters.update(updated_options)

    def create_queryset(self):
        #Sanitize parameters
        #parameters = {field:self.parameters[field] for fields in allowable_fields \
        #    for field in (fields[0], fields[2]) if field in self.parameters}
        #assert 0, allowable_fields

        added =[]
        for form, (field, search_type, convert) in allowable_fields.iteritems():
            value = self.parameters.get(form)

            if not value:
                continue

            search_type = self.parameters.get(search_type, "is")

            self.query.format(field, search_type, value, convert)

        if self.query.has_errors():
            return
        # print self.query
        self.query_set = self.query_set.filter(**self.query)

        if self.unique:
            unique_ids = self.query_set.values("sequence", "taxonomy").annotate(id=Max("id")).values_list("id", flat=True) #.annotate(num_unique=Sum("num_tax"))
            self.count = unique_ids.count()
            #assert 0, self.query_set
        else:
            self.count =  self.query_set.count()

    def sort_query_set(self, unique=False):
        """Sort the query set and paginate results.
        """
        sort_by = self.sort["sort"]
        sort_order = "-" if self.sort["order"] == "desc" else ""
        sort_key = "{}{}".format(sort_order, sort_by)
        self.query_set = self.query_set.order_by(sort_key)

    def paginate(self, sequences=None):
        if sequences is None:
            sequences
        page_size = int(self.sort["limit"])
        start = int(self.sort["offset"])
        end = start+page_size

        self.query_set = self.query_set[start:end]

    def unique(self):
        used_taxa = {}
        for seq in self.query_set:
            if not seq.sequence in used_taxa:
                used_taxa[seq.sequence] = [seq.taxonomy.id]
                yield seq
            elif seq.sequence in used_taxa and not seq.taxonomy.id in used_taxa[seq.sequence]:
                used_taxa[seq.sequence].append(seq.taxonomy.id)
                yield seq
            else:
                pass

    def __len__(self):
        return self.query_set.count()

    def count(self):
        return self.query_set.count()

    def get_dict(self):
        self.sort_query_set()

        if self.unique:
            sequences = Indexable(self.unique())
        else:
            sequences = self.query_set

        #Paginate
        page_size = int(self.sort["limit"])
        start = int(self.sort["offset"])
        end = start+page_size
        sequences = sequences[start:end]

        result = [{
            "id":r.id,
            "variant":r.variant.id,
            "core":r.variant.hist_type.id,
            "gene":r.gene, 
            "splice":r.splice, 
            "taxonomy":r.taxonomy.name.capitalize(), 
            "taxid":str(r.taxonomy.id),
            "score":r.score, 
            "evalue":r.evalue, 
            "header":r.header.replace(r.taxonomy.name.capitalize(), "").replace("[]", "")[:80]
            } for r in sequences]
        return {"total":self.count, "rows":result}

    def get_sunburst(self):
        from browse.management.commands.buildsunburst import build_sunburst
        self.sort_query_set()
        return build_sunburst(self.query_set)

    def get_score_range(self):
        scores = self.query_set.aggregate(scores_min=Min("score"), scores_max=Max("score"))
        return scores["scores_min"], scores["scores_max"]

    def simple_search(self, search_text):
        """Search from simple text box in brand or sequence filter.
        """
        #Search all fields
        try:
            #If search is just a single digit, try to match id
            #Other int values don't make sense to search like this
            value = int(search_text)
            sequence = self.query_set.filter(id=str(value))
            if len(sequence):
                self.query.format("id", "is", value, str)
                return
        except ValueError:
            pass

        #search core_type, variant, old variant names, header if doens't match variant or core_type, taxonomy
        try:
            core_type = Histone.objects.get(id=search_text)
            if self.navbar:
                self.redirect = redirect(core_type)
            else:
                self.query.format("variant__core_type__id", "is", core_type.id, str)
            return
        except Histone.DoesNotExist:
            pass
        
        try:
            variant = Variant.objects.get(id=search_text)
            if self.navbar:
                self.redirect = redirect(variant)
            else:
                self.query.format("variant_id", "is", variant.id, str)
            return
        except Variant.DoesNotExist:
            pass
        
        try:
            #Searches H2A.Z.1.s1
            sequences = Sequence.objects.filter(full_variant_name=search_text)
            if sequences.count() > 0:
                self.query_set &= sequences #Intersect query sets
                #self.query.format("id", "in (comma separated)", ",".join(sequences.values_list("id", flat=True), str)
                return
        except:
            pass

        try:
            variant = OldStyleVariant.objects.get(name=search_text).updated_variant
            
            if self.navbar:
                self.redirect = redirect(variant)
            else:
                self.query.format("variant_id", "is", variant.id, str)
            return
        except OldStyleVariant.DoesNotExist:
            pass

        try:
            #Search species
            taxas = Taxonomy.objects.filter(name__icontains=search_text)
            if taxas.count() > 0:
                self.query.format("taxonomy", "in (comma separated)", ",".join(map(str, taxas.values_list("id", flat=True))), str)
                return
        except Taxonomy.DoesNotExist:
            pass

        try:
            taxon = Taxonomy.objects.get(name=search_text)
            if taxon.type_name != "scientific name":
                #Convert homonym into real taxon
                taxon = taxon.get_scientific_names()[0]
            taxas = taxon.children.filter(rank__name="species").values_list("pk")
            if sequences.count() > 0:
                self.query.format("taxonomy", "in (comma separated)", ",".join(map(str, taxas.values_list("id", flat=True))), str)
                return
        except Taxonomy.DoesNotExist:
            pass
            
        headers = Sequence.objects.filter(header__icontains=search_text)

        if headers.count() > 0:
            self.query_set &= headers
        else:
            #Search sequence moetifs if everything else fails
            self.query.format("sequence", "contains (case-insensitive)", search_text, str)

        return

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
                        v = conv_type(v.strip())
                        if isinstance(v, list):
                            values += v
                        else:
                            values.append(v)
                else:
                    self.errors["Must include a dash if searching 'range'"] += 1
                value = values
            elif search_type.endswith("in"):
                values = []
                for v in value.split(","):
                    v = conv_type(v.strip())
                    if isinstance(v, list):
                        values += v
                    else:
                        values.append(v)
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





