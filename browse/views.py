import sys
import json

from django.shortcuts import render
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.http import JsonResponse
from django.shortcuts import redirect
from django.core.urlresolvers import reverse
from django.shortcuts import get_list_or_404
from django.core.exceptions import ObjectDoesNotExist
from django.templatetags.static import static

from browse.forms import AdvancedFilterForm, AnalyzeFileForm
from browse.search import HistoneSearch
from browse.process_upload import process_upload, InvalidFASTA

from colour import Color

#Django libraires
from browse.models import *
from djangophylocore.models import *

#BioPython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from django.db.models import Min, Max, Count

#Set2 Brewer, used in variant colors
colors7 = [
    "#000000", #Fix for cononical
    "#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
    "#a6d854",
    "#ffd92f",
    "#e5c494"
    ]

colors = [
    "#8dd3c7",
    "#E6E600",
    "#bebada",
    "#fb8072",
    "#80b1d3",
    "#fdb462",
    "#b3de69",
    "#fccde5",
    "#d9d9d9",
    "#bc80bd",
    "#ccebc5",
    "#ffed6f",
]

def help(request):
    data = {
        "filter_form":AdvancedFilterForm(), 
        "original_query":{},
        "current_query":{}
    }
    return render(request, 'help.html', data)

def browse_types(request):
    """Home"""
    data = {
        "filter_form":AdvancedFilterForm(), 
        "original_query":{},
        "current_query":{}
    }
    return render(request, 'browse_types.html', data)

def browse_variants(request, histone_type):
    try:
        hist_type = Histone.objects.get(id=histone_type)
    except:
        return "404"

    variants = hist_type.variants.annotate(num_sequences=Count('sequences')).order_by("id").all().values_list("id", "num_sequences", "taxonmic_span")

    curated_variants = hist_type.variants.filter(sequences__reviewed=True).annotate(num_sequences=Count('sequences')).order_by("id").all().values_list("num_sequences", flat=True)

    variants = [(id, num_curated, num_all, ", ".join(Variant.objects.get(id=id).old_names.values_list("name", flat=True)), tax_span, color) \
        for (id, num_all, tax_span), num_curated, color in zip(variants, curated_variants, colors)]

    data = {
        "histone_type": histone_type,
        "histone_description": hist_type.description,
        "browse_section": "type",
        "name": histone_type,
        "variants": variants,
        "tree_url": "browse/trees/{}.xml".format(hist_type.id),
        "seed_url": reverse("browse.views.get_seed_aln_and_features", args=[hist_type.id]),
        "filter_form": AdvancedFilterForm(),
    }

    #Store sequences in session, accesed in get_sequence_table_data
    original_query = {"id_hist_type":histone_type}
    if request.method == "POST":
        query = request.POST.copy()
    else:
        query = original_query

    #Here we do the search
    result = HistoneSearch(query)

    data["original_query"] = original_query
    if result.errors:
        query = original_query
        data["filter_errors"] = result.errors
    data["current_query"] = query
    
    return render(request, 'browse_variants.html', data)

def browse_variant(request, histone_type, variant):
    try:
        variant = Variant.objects.get(id=variant)
    except:
        return "404"

    green = Color("#66c2a5")
    red = Color("#fc8d62")
    color_range = map(str, red.range_to(green, 12))

    scores = Sequence.objects.filter(variant__id=variant).filter(all_model_scores__used_for_classification=True).annotate(score=Max("all_model_scores__score")).aggregate(max=Max("score"), min=Min("score"))

    data = {
        "hist_type": variant.hist_type.id,
        "variant": variant.id,
        "name": variant.id,
        "sunburst_url": static("browse/sunbursts/{}/{}.json".format(variant.hist_type.id, variant.id)),
        "seed_url": reverse("browse.views.get_seed_aln_and_features", args=[variant.id]),
        "colors":color_range,
        "score_min":scores["min"],
        "score_max":scores["max"],
        "browse_section": "variant",
        "description": variant.description,
        "alternate_names": ", ".join(variant.old_names.values_list("name", flat=True)),
        "filter_form": AdvancedFilterForm(),
    }

    original_query = {"id_variant":variant.id}
    if request.method == "POST":
        query = request.POST.copy()
    else:
        query = original_query
    result = HistoneSearch(query)

    data["original_query"] = original_query
    if result.errors:
        query = original_query
        data["filter_errors"] = result.errors
    data["current_query"] = query

    return render(request, 'browse_variant.html', data)

def search(request):
    data = {"filter_form": AdvancedFilterForm()}
    
    if request.method == "POST": 
        query = request.POST.copy()
    else:
        query = request.GET.copy()
    result = HistoneSearch(
        query,
        navbar="search" in query.keys())

    if query.get("reset", True):
        request.session["original_query"] = query

    data["original_query"] = request.session.get("original_query", query)
    data["current_query"] = query

    if len(result.errors) == 0:
        data["result"] = True
    else:
        data["filter_errors"] = result.errors

    if result.redirect: 
        return result.redirect

    green = Color("#66c2a5")
    red = Color("#fc8d62")
    data["colors"] = map(str, red.range_to(green, 12))

    data["score_min"], data["score_max"] = result.get_score_range()

    return render(request, 'search.html', data)

def basket(request):
    data = {
        "filter_form":AdvancedFilterForm(), 
        "original_query":{},
        "current_query":{}
    }
    return render(request, 'basket.html', data)

def analyze(request):
    data = {
        "filter_form":AdvancedFilterForm(), 
        "original_query":{},
        "current_query":{}
    }
    if request.method == "POST":
        if request.POST.get("sequences"):
            format = "text"
            sequences = request.POST["sequences"]
        elif request.POST.get("file"):
            format="file"
            sequences = request.POST["file"]

        try:
            data["result"] = process_upload(sequences, format, request)
        except InvalidFASTA as e:
            data["error"] = "{}: {}".format(e.__class__.__name__, e.message)
            data["analyze_form"] = AnalyzeFileForm()

        data["search_type"] = type
    else:
        data["analyze_form"] = AnalyzeFileForm()

    return render(request, 'analyze.html', data)

def get_sequence_table_data(request):
    """Downloads the previous search and converts into json required by Bootstrap table
    """
    
    if request.method == "GET":
        parameters =  request.GET.dict()
    else:
        assert 0, request.method
        #Returning 'false' stops Bootstrap table
        parameters = []
    
    results = HistoneSearch(parameters)

    if len(results.errors) > 0:
        #Returning 'false' stops Bootstrap table
        return "false"

    result = results.get_dict()
    result["parameters"] = results.parameters
    result["method"] = request.method

    return JsonResponse(result)

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
        rows[i]["variant"] = variant #"{} (T:{})".format(variant, threshold)
        for id in ids:
            rows[i][id] = "n/a"
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
                if rows[indices[score.variant.id]][id] == "n/a" or score.score > rows[indices[score.variant.id]][id]:
                    rows[indices[score.variant.id]][id] = score.score
                    rows[indices[score.variant.id]]["data"]["above_threshold"][id] = score.score>=threshold
                    rows[indices[score.variant.id]]["data"]["this_classified"][id] = score.used_for_classification
            
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

def get_aln_and_features(request, ids=None):
    from tools.hist_ss import templ, get_hist_ss_in_aln, get_hist_ss
    from tools.L_shade_hist_aln import write_alignments
    import subprocess
    import StringIO
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.AlignInfo import SummaryInfo
    from Bio.SeqRecord import SeqRecord

    if ids is None and request.method == "GET" and "id" in request.GET:
        ids = request.GET.getlist("id")
        sequences = Sequence.objects.filter(id__in=ids[:50])
        download = False
        upload = False
    elif request.GET.get("download", False) == "true":
        download = True
        upload = False
    elif request.GET.get("upload", False) == "true":
        sequences = request.session.get("uploaded_sequences", [])
        sequences = [Sequence(
            id=s["id"], 
            variant=Variant.objects.get(id=s["variant"]),
            sequence=s["sequence"],
            taxonomy=Taxonomy.objects.get(name=s["taxonomy"])) for s in sequences]
        upload = True
        download = False
    else:
        #Returning 'false' stops Bootstrap table
        return "false"

    if not download:
        if len(sequences) == 0:
            return None, None
        elif len(sequences) == 1:
            #Already aligned to core histone
            seq = sequences[0]
            hist_type = seq.variant.hist_type.id
            #let's load the corresponding canonical
            try:
                canonical=Sequence.objects.filter(variant_id='canonical'+str(seq.variant.hist_type),reviewed=True,taxonomy=seq.taxonomy)[0]
            except:
                canonical = Sequence(id="0000|xenopus|canonical{}".format(hist_type), sequence=str(templ[hist_type]))
            sequences = [canonical, seq]
            
        else:
            seq = sequences[0]
            try:
                hist_type = max(
                   [(hist, sequences.filter(variant__hist_type_id=hist).count()) for hist in ["H2A", "H2B", "H3", "H4", "H1"]],
                   key=lambda x:x[1]
                   )[0]
            except ValueError:
                hist_type = "Unknown"
        
        muscle = os.path.join(os.path.dirname(sys.executable), "muscle")
        process = subprocess.Popen([muscle], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        sequences = "\n".join([s.format() for s in sequences])
        aln, error = process.communicate(sequences)
        seqFile = StringIO.StringIO()
        seqFile.write(aln)
        seqFile.seek(0)
        sequences = list(SeqIO.parse(seqFile, "fasta")) #Not in same order, but does it matter?
        msa = MultipleSeqAlignment(sequences)
        a = SummaryInfo(msa)
        cons = Sequence(id="Consensus", sequence=a.dumb_consensus(threshold=0.1, ambiguous='X').tostring())

        save_dir = os.path.join(os.path.sep, "tmp", "HistoneDB")
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        if not hist_type == "H1":
            #TODO: we need to test if gff annotation works correctly. MSA needs numbering with respect to MSA or individual seqs as TexShade???
            hv,ss = get_hist_ss_in_aln(msa, hist_type=hist_type, save_dir=save_dir, debug=False)
            # features = Features.from_dict(cons, ss)
            #Indeed the features are for MSA, we need just to use the name of last sequence:
            features = Features.from_dict(seq, ss)
        else:
            features = ""
        #A hack to avoid two canonical seqs
        unique_sequences = [sequences[0]] if sequences[0].id==sequences[1].id else sequences
        # doing the Sequence.short_description work
        #Note that the gffs are also generated with the short description not
        sequences = [{"name":Sequence.long_to_short_description(s.id), "seq":s.seq.tostring()} for s in unique_sequences]
        # sequences = [{"name":s.id, "seq":s.seq.tostring()} for s in sequences]
        #Uncomment to add consesus as first line
        # sequences.insert(0, cons.to_dict())

        request.session["calculated_msa_seqs"] = sequences
        request.session["calculated_msa_features"] = features.to_dict() if features else {}

        result = {"seqs":sequences, "features":features.full_gff() if features else ""}
        return JsonResponse(result, safe=False) 
    else:
        format = request.GET.get("format", "json")
        response = HttpResponse(content_type='text')
        response['Content-Disposition'] = 'attachment; filename="sequences.{}"'.format(format)

        
        sequences = request.session.get("calculated_msa_seqs", [])
        features_dict = request.session.get("calculated_msa_features", None)
        features = Features.from_dict(Sequence("Consensus"), features_dict) if features_dict else None

        if format == "fasta":
            for s in sequences:
                print >> response, ">{}\n{}".format(s["name"], s["seq"])
        elif format == "gff":
            response.write(features.full_gff() if features else "")
        elif format == "pdf":
            response.write("Sorry PDF Feature is still in development")
            #return redirect("{}.pdf".format(seed_file))
            save_dir = os.path.join(os.path.sep, "tmp", "HistoneDB")
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
                os.chmod(save_dir,0o777)

            aln = MultipleSeqAlignment([SeqRecord(Seq(s["seq"]), id=s["name"]) for s in sequences[1:]])
            result_pdf = write_alignments(
                [aln], 
                save_dir = save_dir
            )

            with open(result_pdf) as pdf:
                response.write(pdf.read())

            #Cleanup
            os.remove(result_pdf)
        else:
            #Default format is json
            result = {"seqs":sequences, "features":features.full_gff() if features else ""}
            response.write(json.dumps(result))

        return response

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

def get_seed_aln_and_features(request, seed):
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.AlignInfo import SummaryInfo

    seed_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "seeds")
    try:
        histone = Histone.objects.get(id=seed)
        seed_file = os.path.join(seed_file, "{}".format(histone.id))
    except Histone.DoesNotExist:
        try:
            variant = Variant.objects.get(id=seed)
            seed_file = os.path.join(seed_file, variant.hist_type.id, "{}".format(variant.id))
        except Variant.DoesNotExist:
            return HttpResponseNotFound('<h1>No histone variant with name {}</h1>'.format(seed))

    download = request.GET.get("download", False) == "true"

    format = request.GET.get("format", "json")

    try:
        limit = int(request.GET.get("limit", 0))
    except ValueError:
        limit = 0

    consensus = request.GET.get("consensus", False)

    if not consensus in ["limit", "all", False]:
        consensus = "all"

    response = HttpResponse(content_type='text')

    if download:
        response['Content-Disposition'] = 'attachment; filename="{}.{}"'.format(seed, format)

    sequences = SeqIO.parse("{}.fasta".format(seed_file), "fasta")

    if consensus:
        sequences = [s for i, s in enumerate(sequences) if consensus == "all" or (consensus == "limit" and i < limit)]
        msa = MultipleSeqAlignment(sequences)
        a = SummaryInfo(msa)
        sequences.insert(0, SeqRecord(id="Consensus", description="", seq=a.dumb_consensus(threshold=0.1, ambiguous='X')))
        limit = limit+1 if limit > 0 else 0

    def limited_seqs():
        for i, seq in enumerate(sequences):
            if not consensus or consensus == "limit" or (limit > 0 and i < limit):
                yield seq

    with open("{}.gff".format(seed_file)) as feature_file:
        features = feature_file.read()

    if format == "fasta":
        SeqIO.write(limited_seqs(), response, "fasta")
    elif format == "gff":
        response.write(features)
    elif format == "pdf":
        with open("{}.pdf".format(seed_file)) as pdf:
            response.write(pdf.read())
    else:
        #Default format is json
        sequences = [{"name":seq.id, "seq":seq.seq.tostring()} for seq in limited_seqs()]
        result = {"seqs":sequences, "features":features}
        response.write(json.dumps(result))

    return response

def get_sunburst_json(request, parameters=None):
    """
    """
    if parameters and isinstance(parameters, dict):
        query = parameters
    elif request.method == "POST": 
        query = request.POST.copy()
    else:
        query = request.GET.copy()
    
    if query:
        result = HistoneSearch(query)
        sunburst = result.get_sunburst()
        return JsonResponse(sunburst, safe=False)
    else:
        raise Http404("No taxonomy distribution sunburst for query") 
