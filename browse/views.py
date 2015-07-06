from django.shortcuts import render
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.http import JsonResponse
from django.shortcuts import redirect
from django.shortcuts import get_list_or_404

from browse.forms import AdvancedFilterForm, AnalyzeFileForm
from browse.search import HistoneSearch
from browse.process_upload import process_upload

from colour import Color

#Django libraires
from browse.models import *
from djangophylocore.models import *

#BioPython
from Bio import SeqIO

from django.db.models import Min, Max, Count

#Set2 Brewer, used in variant colors
colors = [
    "#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
    "#a6d854",
    "#ffd92f",
    "#e5c494"
    ]

def treemap(request):
    return render(request, 'circle.html', {"starburst_url": "data/type/json/{}/species/".format("H2A")})

def help(request):
    return render(request, 'help.html', {})

def browse_types(request):
    """Home"""
    return render(request, 'browse_types.html', {"filter_form":AdvancedFilterForm()})

def browse_variants(request, histone_type):
    try:
        core_histone = Histone.objects.get(id=histone_type)
    except:
        return "404"

    #Store sequences in session, accesed in get_sequence_table_data
    HistoneSearch(request, {"id_core_histone":histone_type}, reset=True)

    variants = core_histone.variants.annotate(num_sequences=Count('sequences')).order_by("id").all().values_list("id", "num_sequences")
    variants = [(id, num, color) for (id, num), color in zip(variants, colors)]


    data = {
        "histone_type": histone_type,
        "histone_description": core_histone.description,
        "browse_section": "type",
        "name": histone_type,
        "variants": variants,
        "tree_url": "browse/trees/{}.xml".format(core_histone.id),
        "seed_file":"browse/seeds/{}.fasta".format(core_histone.id),
        "filter_form": AdvancedFilterForm(),
    }

    return render(request, 'browse_variants.html', data)

def browse_variant(request, histone_type, variant):
    try:
        variant = Variant.objects.get(id=variant)
    except:
        return "404"

    HistoneSearch(request, {"id_variant":variant.id}, reset=True)

    green = Color("#66c2a5")
    red = Color("#fc8d62")
    color_range = map(str, red.range_to(green, 12))

    scores = Sequence.objects.filter(variant__id=variant).filter(all_model_scores__used_for_classifiation=True).annotate(score=Max("all_model_scores__score")).aggregate(max=Max("score"), min=Min("score"))

    data = {
        "core_type": variant.core_type.id,
        "variant": variant.id,
        "name": variant.id,
        "sunburst_url": "browse/sunbursts/{}/{}.json".format(variant.core_type.id, variant.id),
        "seed_file":"browse/seeds/{}/{}.fasta".format(variant.core_type.id, variant.id),
        "colors":color_range,
        "score_min":scores["min"],
        "score_max":scores["max"],
        "browse_section": "variant",
        "description": variant.description,
        "filter_form": AdvancedFilterForm(),
    }
    return render(request, 'browse_variant.html', data)

def search(request):
    data = {}
    if request.method == "POST": 
        result = HistoneSearch(
            request, 
            request.POST.copy(), 
            reset=request.POST.get("reset", True), 
            navbar="search" in request.POST)
    else:
        #Show all sequence
        result = HistoneSearch.all(request)
        
    if result.redirect: 
        return result.redirect
        
    return render(request, 'search.html', {"result":result, "filter_form": AdvancedFilterForm(),})

def analyze(request):
    data = {"filter_form":AdvancedFilterForm()}
    if request.method == "POST":
        type = request.POST.get("id_type_0")
        if request.POST.get("sequences"):
            format = "text"
            sequences = request.POST["sequences"]
        elif request.POST.get("file"):
            format="file"
            sequences = request.POST["file"]
        data["result"] = process_upload(type, sequences, format)
        data["search_type"] = type
    else:
        data["analyze_form"] = AnalyzeFileForm(initial={"type":"blastp"})
    return render(request, 'analyze.html', data)

def get_sequence_table_data(request):
    """Downloads the previos search and converts into json required by Bootstrap table
    """

    if request.method == "GET":
        parameters =  request.GET.dict()
    else:
        #Returning 'false' stops Bootstrap table
        parameters = []

    #Continues to filter previous search, unless paramters contains key 'reset'
    results = HistoneSearch(request, parameters)

    if len(results.errors) > 0:
        #Returning 'false' stops Bootstrap table
        return "false"
    
    return JsonResponse(results.get_dict())

def get_all_scores(request, ids=None):
    if ids is None and request.method == "GET" and "id" in request.GET:
        ids = request.GET.getlist("id")
    else:
        #Returning 'false' stops Bootstrap table
        return "false"
    
    variants = list(Variant.objects.all().order_by("id").values_list("id", "hmmthreshold"))
    indices = {variant: i for i, (variant, threshold) in enumerate(variants)}
    rows = [{} for _ in xrange(len(variants))]
    for i, (variant, threshold) in enumerate(variants):
        rows[i]["variant"] = "{} (T:{})".format(variant, threshold)
        for id in ids:
            rows[i][id] = "<0"
        rows[i]["data"] = {}
        rows[i]["data"]["above_threshold"] = {id:False for id in ids}
        rows[i]["data"]["this_classified"] = {id:False for id in ids}
    
    for i, id in enumerate(ids):
        try:
            sequence = Sequence.objects.get(id=id)
        except:
            return "404"
        classified_variant = sequence.variant.id
        
        scores = sequence.all_model_scores.all().order_by("variant__id")

        for j, score in enumerate(scores):
            if score.variant.id in indices:
                threshold = score.variant.hmmthreshold
                if rows[indices[score.variant.id]][id] == "<0" or score.score > rows[indices[score.variant.id]][id]:
                    rows[indices[score.variant.id]][id] = score.score
                    rows[indices[score.variant.id]]["data"]["above_threshold"][id] = score.score>=threshold
                    rows[indices[score.variant.id]]["data"]["this_classified"][id] = score.used_for_classifiation
            
    return JsonResponse(rows, safe=False)

def get_all_sequences(request, ids=None):
    if ids is None and request.method == "GET" and "id" in request.GET:
        ids = request.GET.getlist("id")
    else:
        #Returning 'false' stops Bootstrap table
        return "false"

    format = request.GET.get("format", "json")
    download = request.GET.get("download", "false") == "true"

    sequences = Sequence.objects.filter(id__in=ids[:50]) 
    if format == "fasta":
        response = HttpResponse(content_type='text')

        if download:
            response['Content-Disposition'] = 'attachment; filename="histone_variants.fasta"'

        for s in sequences:
            response.write(str(s))

        return response

    else:
        sequences = [s.to_dict() for s in sequences]
        return JsonResponse(sequences, safe=False)

def get_sequence_features(request, ids=None):
    if ids is None and request.method == "GET" and "id" in request.GET:
        ids = request.GET.getlist("id")
    else:
        #Returning 'false' stops Bootstrap table
        return "false"

    download = request.GET.get("download", "false") == "true"
    sequences = Sequence.objects.filter(id__in=ids[:50])

    response = HttpResponse(content_type='text')

    if download:
        response['Content-Disposition'] = 'attachment; filename="histone_annotations.fasta"'

    response.write(Features.gff_colors())
    for s in sequences:
        response.write(str(s.features))

    return response


def get_starburst_json(request, browse_type, search, debug=False):
    """
    """
    if ids is None and request.method == "GET" and "id" in request.GET:
        ids = request.GET.getlist("id")
    else:
        #Returning 'false' stops Bootstrap table
        return "false"
    

def create_sunburst(root, sunburst=None, level=0):
    if sunburst is None:
        sunburst = {"name":root.name, "children":[]}

    for curr_taxa in root.direct_children.filter(type_name="scientific name").all():
        #print "{}{}".format(" "*level, curr_taxa.name),
        child = {"name":curr_taxa.name, "children":[]}
        #print child
        sunburst["children"].append(create_sunburst(curr_taxa, child, level=level+1))
    return sunburst

def get_phyloxml_of_type(request):
    return None
