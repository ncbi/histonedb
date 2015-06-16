import sys
import os
import StringIO
import shlex

import uuid
from Bio import SeqIO, SearchIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

from browse.models import *
from django.conf import settings

import subprocess

from django.db.models import Q
from django.db.models import Max, Min, Count

class TooManySequences(RuntimeError):
    pass

class InvalidFASTA(RuntimeError):
    pass

def process_upload(type, sequences, format):
    if format == "file":
        processed_sequences = list(SeqIO.parse(sequences, "fasta"))
        
    elif format == "text":
        seq_file = StringIO.StringIO()
        seq_file.write(sequences)
        seq_file.seek(0)
        processed_sequences = list(SeqIO.parse(seq_file, "fasta"))

    if len(processed_sequences) > 50:
        raise TooManySequences

    sequences = "\n".join([seq.format("fasta") for seq in processed_sequences])

    if type == "blastp":
        result = upload_blastp(sequences, len(processed_sequences))
    elif type == "hmmer":
        result = upload_hmmer(sequences)
    else:
        assert 0, type

    return result

def upload_blastp(seq_file, num_seqs):
    output= os.path.join("/", "tmp", "{}.xml".format(seq_file[:-6]))
    blastp_cline = NcbiblastpCommandline(
        db=os.path.join(settings.STATIC_ROOT, "browse", "blast", "HistoneDB_sequences.fa"), 
        evalue=0.01, outfmt=5)
    out, err = blastp_cline(stdin=seq_file)
    #cmd = str(blastp_cline)
    #process = subprocess.Popen(shlex.shlex(cmd), stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    #assert 0, "'{}'".format(seq_file)
    #output, error = process.communicate(seq_file)
    blastFile = StringIO.StringIO()
    blastFile.write(out)
    blastFile.seek(0)

    result = []
    for i, blast_record in enumerate(NCBIXML.parse(blastFile)):
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
                    "gi":str(sequence.id), 
                    "variant":str(sequence.variant_id), 
                    "gene":str(sequence.gene) if sequence.gene else "-", 
                    "splice":str(sequence.splice) if sequence.splice else "-", 
                    "species":str(sequence.taxonomy.name), 
                    "score":str(sequence.score), 
                    "evalue":str(sequence.evalue), 
                    "header":str(sequence.header), 
                    "search_e":str(search_evalue),
                })
    return result

def upload_hmmer(seq_file, evalue=10):
    """
    """
    temp_seq_path = "/tmp/{}.fasta".format(uuid.uuid4())
    with open(temp_seq_path, "w") as seqs:
        seqs.write(seq_file);

    variantdb = os.path.join(settings.STATIC_ROOT, "browse", "hmms", "combined_variants.hmm")
    coredb = os.path.join(settings.STATIC_ROOT, "browse", "hmms", "combined_cores.hmm")
    hmmsearch = os.path.join(os.path.dirname(sys.executable), "hmmsearch")

    results = {}
    variants = []

    for i, db in enumerate((variantdb, coredb)):
        process = subprocess.Popen([hmmsearch, "-E", str(evalue), "--notextw", db, temp_seq_path], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        output, error = process.communicate()
        hmmerFile = StringIO.StringIO()
        hmmerFile.write(output)
        hmmerFile.seek(0)
        for variant_query in SearchIO.parse(hmmerFile, "hmmer3-text"):
            if i==1: 
                variant = "canonical{}".format(variant_query.id)
            else:
                variant = variant_query.id

            variants.append(variant)

            try:
                variant_model = Variant.objects.get(id=variant)
            except:
                continue

            for hit in variant_query:
                for hsp in hit:
                    if not hit.id in results:
                        results[hit.id] = {"class":"Unknown", "best_score": 0, "scores":{}}
                    if hsp.bitscore >= variant_model.hmmthreshold and hsp.bitscore > results[hit.id]["best_score"]:
                        if i==1 and not (results[hit.id]["class"] == "Unknown" or "canonical" in results[hit.id]["class"]):
                            continue
                        model = variant_model.id

                        results[hit.id]["class"] = variant_model.id
                        results[hit.id]["best_score"] = hsp.bitscore

                    results[hit.id]["scores"][variant_model.id] = (hsp.bitscore, hsp.bitscore>=variant_model.hmmthreshold)

    return results
