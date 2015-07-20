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

from browse.forms import AdvancedFilterForm, AnalyzeFileForm
from browse.search import HistoneSearch
from browse.process_upload import process_upload

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
colors = [
    "#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
    "#a6d854",
    "#ffd92f",
    "#e5c494"
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
        core_histone = Histone.objects.get(id=histone_type)
    except:
        return "404"

    variants = core_histone.variants.annotate(num_sequences=Count('sequences')).order_by("id").all().values_list("id", "num_sequences")
    variants = [(id, num, color) for (id, num), color in zip(variants, colors)]

    data = {
        "histone_type": histone_type,
        "histone_description": core_histone.description,
        "browse_section": "type",
        "name": histone_type,
        "variants": variants,
        "tree_url": "browse/trees/{}.xml".format(core_histone.id),
        "seed_url": reverse("browse.views.get_seed_aln_and_features", args=[core_histone.id]),
        "filter_form": AdvancedFilterForm(),
    }

    #Store sequences in session, accesed in get_sequence_table_data
    original_query = {"id_core_histone":histone_type}
    if request.method == "POST":
        query = request.POST.copy()
    else:
        query = original_query
    
    result = HistoneSearch(request, query, reset=True)

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

    scores = Sequence.objects.filter(variant__id=variant).filter(all_model_scores__used_for_classifiation=True).annotate(score=Max("all_model_scores__score")).aggregate(max=Max("score"), min=Min("score"))

    data = {
        "core_type": variant.core_type.id,
        "variant": variant.id,
        "name": variant.id,
        "sunburst_url": "browse/sunbursts/{}/{}.json".format(variant.core_type.id, variant.id),
        "seed_url": reverse("browse.views.get_seed_aln_and_features", args=[variant.id]),
        "colors":color_range,
        "score_min":scores["min"],
        "score_max":scores["max"],
        "browse_section": "variant",
        "description": variant.description,
        "filter_form": AdvancedFilterForm(),
    }

    original_query = {"id_variant":variant.id}
    if request.method == "POST":
        query = request.POST.copy()
    else:
        query = original_query
    result = HistoneSearch(request, query, reset=True)

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
        result = HistoneSearch(
            request, 
            query,
            navbar="search" in request.POST)
        if request.POST.get("reset", True):
            request.session["original_query"] = query
        data["original_query"] = request.session.get("original_query", query)
        data["current_query"] = query
        if len(result.errors) == 0:
            data["result"] = True
        else:
            data["filter_errors"] = result.errors
    else:
        #Show all sequence
        result = HistoneSearch.all(request)
        data["current_query"] = data["original_query"] = {}
        data["result"] = True

    if result.redirect: 
        return result.redirect

    return render(request, 'search.html', data)

def analyze(request):
    data = {
        "filter_form":AdvancedFilterForm(), 
        "original_query":{},
        "current_query":{}
    }
    if request.method == "POST":
        type = request.POST.get("id_type_0")
        initial={"type":type}
        if request.POST.get("sequences"):
            format = "text"
            sequences = request.POST["sequences"]
            initial["sequences"] = sequences
        elif request.POST.get("file"):
            format="file"
            sequences = request.POST["file"]
        try:
            data["result"] = process_upload(type, sequences, format)
        except Exception as e:
            data["error"] = "{}: {}".format(e.__class__.__name__, e.message)
            data["analyze_form"] = AnalyzeFileForm(initial=initial)
        data["search_type"] = type
    else:
        data["analyze_form"] = AnalyzeFileForm(initial={"type":"blastp"})
    return render(request, 'analyze.html', data)

def get_sequence_table_data(request):
    """Downloads the previous search and converts into json required by Bootstrap table
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

    unique = "id_unique" in parameters
    
    return JsonResponse(results.get_dict(unique=unique))

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

def get_aln_and_features(request, ids=None):
    from tools.hist_ss import templ, get_hist_ss_in_aln
    from tools.L_shade_hist_aln import write_alignments
    import subprocess
    import StringIO
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.AlignInfo import SummaryInfo

    if ids is None and request.method == "GET" and "id" in request.GET:
        ids = request.GET.getlist("id")
        download = False
    elif request.GET.get("download", False) == "true":
        download = True
    else:
        #Returning 'false' stops Bootstrap table
        return "false"

    if not download:
        sequences = Sequence.objects.filter(id__in=ids[:50])
        if len(sequences) == 0:
            return None, None
        elif len(sequences) == 1:
            #Already aligned to core histone
            canonical = {"name":"Consensus".format(sequences.first().variant.core_type.id), "seq":str(templ[sequences.first().variant.core_type.id])}
            seq = sequences.first()
            sequences = [canonical, seq.to_dict()]
            try:
                features = seq.features 
            except ObjectDoesNotExist:
                features = ""
            finally:
                if "H1" in seq.variant.id:
                    features = ""
        else:
            try:
                hist_type = max(
                   [(hist, sequences.filter(variant__core_type_id=hist).count()) for hist in ["H2A", "H2B", "H3", "H4", "H1"]],
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
                hv,ss = get_hist_ss_in_aln(msa, hist_type=hist_type, save_dir=save_dir, debug=False)
                features = Features.from_dict(cons, ss)
            else:
                features = ""

            sequences = [{"name":s.id, "seq":s.seq.tostring()} for s in sequences]
            sequences.insert(0, cons.to_dict())

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

    seed_file = os.path.join(settings.STATIC_ROOT, "browse", "seeds")
    try:
        histone = Histone.objects.get(id=seed)
        seed_file = os.path.join(seed_file, "{}".format(histone.id))
    except Histone.DoesNotExist:
        try:
            variant = Variant.objects.get(id=seed)
            seed_file = os.path.join(seed_file, variant.core_type.id, "{}".format(variant.id))
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
