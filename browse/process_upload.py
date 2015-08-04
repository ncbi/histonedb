import sys
import os
import StringIO
import shlex

import uuid
from Bio import SeqIO, SearchIO
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import Alphabet

from browse.models import Histone, Variant, Sequence
from django.conf import settings

import subprocess

from django.db.models import Q
from django.db.models import Max, Min, Count

class TooManySequences(Exception):
    pass

class InvalidFASTA(Exception):
    pass

class InvalidSearchType(Exception):
    pass

def process_upload(type, sequences, format, request):
    assert format in ["file", "text"], "Invalid format: {}. Must be either 'file' or 'text'.".format(format)

    if format == "text":
        seq_file = StringIO.StringIO()
        seq_file.write(sequences)
        seq_file.seek(0)
        sequences = seq_file

    processed_sequences = []
        
    for i, seq in enumerate(SeqIO.parse(seq_file, "fasta", alphabet=IUPAC.protein)):
        if i >= 50:
            raise TooManySequences()
        elif not Alphabet._verify_alphabet(seq.seq):
            raise InvalidFASTA("Sequence {} is not a protein.".format(seq.id))

        processed_sequences.append(seq)

    if len(processed_sequences) == 0:
        raise InvalidFASTA("No sequences parsed.")

    sequences = "\n".join([seq.format("fasta") for seq in processed_sequences])

    if type == "blastp":
        result = upload_blastp(sequences, len(processed_sequences))
    elif type == "hmmer":
        result = upload_hmmer(processed_sequences, len(processed_sequences))
    else:
        raise InvalidSearchType("Muse search only 'blastp' or 'hmmer'")

    return result

def upload_blastp(seq_file, num_seqs):
    blastp = os.path.join(os.path.dirname(sys.executable), "blastp")
    output= os.path.join("/", "tmp", "{}.xml".format(seq_file[:-6]))
    blastp_cline = NcbiblastpCommandline(
        cmd=blastp,
        db=os.path.join(settings.STATIC_ROOT, "browse", "blast", "HistoneDB_sequences.fa"), 
        evalue=0.01, outfmt=5)
    out, err = blastp_cline(stdin=seq_file)
    blastFile = StringIO.StringIO()
    blastFile.write(out)
    blastFile.seek(0)
    
    results = []
    for i, blast_record in enumerate(NCBIXML.parse(blastFile)):
        result = []
        for alignment in blast_record.alignments:
            try:
                gi = alignment.hit_def.split("|")[0]
            except IndexError:
                continue
            for hsp in alignment.hsps:
                sequence = Sequence.objects.filter(
                        (~Q(variant__id="Unknown") & Q(all_model_scores__used_for_classifiation=True)) | \
                        (Q(variant__id="Unknown") & Q(all_model_scores__used_for_classifiation=False)) \
                    ).annotate(
                        num_scores=Count("all_model_scores"), 
                        score=Max("all_model_scores__score"),
                        evalue=Min("all_model_scores__evalue")
                    ).get(id=gi)
                search_evalue = hsp.expect
                result.append({
                    "id":str(sequence.id), 
                    "variant":str(sequence.variant_id), 
                    "gene":str(sequence.gene) if sequence.gene else "-", 
                    "splice":str(sequence.splice) if sequence.splice else "-", 
                    "taxonomy":str(sequence.taxonomy.name), 
                    "score":str(sequence.score), 
                    "evalue":str(sequence.evalue), 
                    "header":str(sequence.header), 
                    "search_e":str(search_evalue),
                })
        if not result:
            raise InvalidFASTA("No blast hits for {}.".format(blast_record.query))
        results.append(result)
    if not results:
        raise InvalidFASTA("No blast hits.")

    return results

def upload_hmmer(seq_file, num_seqs, evalue=10):
    """
    """
    save_dir = os.path.join(os.path.sep, "tmp", "HistoneDB")
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    temp_seq_path = os.path.join(save_dir, "{}.fasta".format(uuid.uuid4()))
    with open(temp_seq_path, "w") as seqs:
        for s in seq_file:
            SeqIO.write(s, seqs, "fasta");

    variantdb = os.path.join(settings.STATIC_ROOT, "browse", "hmms", "combined_variants.hmm")
    coredb = os.path.join(settings.STATIC_ROOT, "browse", "hmms", "combined_cores.hmm")
    hmmsearch = os.path.join(os.path.dirname(sys.executable), "hmmsearch")

    results = {}
    variants = []

    variants = list(Variant.objects.all().order_by("id").values_list("id", "hmmthreshold"))
    indices = {variant: i for i, (variant, threshold) in enumerate(variants)}
    seqs_index = {seq.id:i for i, seq in enumerate(seq_file)}
    ids = map(lambda s: s.id, seq_file)
    rows = [{} for _ in xrange(len(variants))]
    classifications = {s.id:"Unknown" for s in seq_file}
    for i, (variant, threshold) in enumerate(variants):
        rows[i]["variant"] = "{} (T:{})".format(variant, threshold)
        for id in ids:
            rows[i][id] = "n/a"
        rows[i]["data"] = {}
        rows[i]["data"]["above_threshold"] = {id:False for id in ids}
        rows[i]["data"]["this_classified"] = {id:False for id in ids}

    for i, db in enumerate((variantdb, coredb)):
        process = subprocess.Popen([hmmsearch, "-E", str(evalue), "--notextw", db, temp_seq_path], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        output, error = process.communicate()
        hmmerFile = StringIO.StringIO()
        hmmerFile.write(output)
        hmmerFile.seek(0)
        for variant_query in SearchIO.parse(hmmerFile, "hmmer3-text"):
            if i==1: 
                variant = "unclassified{}".format(variant_query.id)
            else:
                variant = variant_query.id

            print variant

            variants.append(variant)

            try:
                variant_model = Variant.objects.get(id=variant)
            except Variant.DoesNotExist:
                continue

            for hit in variant_query:
                print "Hit is", hit.id
                for hsp in hit:
                    if hsp.bitscore>=variant_model.hmmthreshold and \
                         (classifications[hit.id] == "Unknown" or \
                          float(hsp.bitscore) >= rows[indices[classifications[hit.id]]][hit.id]):
                        if i==1 and not (classifications[hit.id] == "Unknown" or "unclassified" in classifications[hit.id]):
                            #Skip canoninical score if already classfied as a variant
                            continue
                        classifications[hit.id] = variant

                        if not classifications[hit.id] == "Unknown":
                            for row in rows:
                                row["data"]["this_classified"][hit.id] = False
                        rows[indices[variant]]["data"]["this_classified"][hit.id] = True
                    
                    rows[indices[variant]]["data"]["above_threshold"][hit.id] = float(hsp.bitscore)>=variant_model.hmmthreshold
                    rows[indices[variant]][hit.id] = float(hsp.bitscore)
    
    classifications = [(id, classifications[id]) for id in ids]

    #Cleanup
    os.remove(temp_seq_path)
    
    return classifications, ids, rows
