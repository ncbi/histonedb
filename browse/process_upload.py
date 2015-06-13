import sys
import os
import StringIO

import uuid
from Bio import SeqIO, SearchIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

from browse.models import *
from django.conf import settings

import subprocess

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
        result = upload_blastp(sequences)
    elif type == "hmmer":
        result = upload_hmmer(sequences)
    else:
        assert 0, type

    return result

def upload_blastp(seq_file):
    output= os.path.join("/", "tmp", "{}.xml".format(seq_file[:-6]))
    blastx_cline = NcbiblastxCommandline(
        query=seq_file, db=os.path.join(settings.STATIC_ROOT, "static", "browse", "blast", "HistoneDB.db"), 
        evalue=0.01, outfmt=5, out=output)
    stdout, stderr = blastx_cline()
    result = {"total":0, "rows":[]}
    with open(output) as result_handle:
        for blast_record in NCBIXML.parse(result_handle):
            for alignment in blast_record.alignments:
                try:
                    gi = alignment.title.split("|")[0]
                except IndexError:
                    continue
                for hsp in alignment.hsps:
                    sequence = Sequence.objects.annotate(evalue=Min("scores__evalue"), score=Max("scores__score")).get(id=gi)
                    search_evalue = hsp.expect
                    result["rows"].append({
                        "gi":sequence.id, 
                        "variant":sequence.variant_id, 
                        "gene":sequence.gene, 
                        "splice":sequence.splice, 
                        "species":sequence.taxonomy.name, 
                        "score":sequence.score, 
                        "evalue":sequence.evalue, 
                        "header":sequence.header, 
                        "search_e":search_evalue
                    })
                    result["total"] += 1
    os.remove(output)
    return result

def upload_hmmer(seq_file, evalue=10):
    """
    """    
    variantdb = os.path.join(settings.STATIC_ROOT, "browse", "hmms", "combined_variants.hmm")
    coredb = os.path.join(settings.STATIC_ROOT, "browse", "hmms", "combined_cores.hmm")
    hmmsearch = os.path.join(os.path.dirname(sys.executable), "hmmsearch")

    assert 0, " ".join([hmmsearch, "-E", str(evalue), "--notextw", variantdb, "-"]) + " < '" + seq_file +"'"
    for db in (variantdb, coredb):
        hmmerFile = StringIO.StringIO()
        process = subprocess.Popen([hmmsearch, "-E", str(evalue), "--notextw", db, "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        output, error = process.communicate(seq_file)
        assert 0, process.poll()

        assert 0, hmmerFile.getvalue()

    results = {}
    variants = []
    for i, hmmer_string in enumerate((variant_output, core_output)):
        hmmerFile = StringIO.StringIO()
        hmmerFile.write(hmmer_string)
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
                        results[hit.id] = {"class":"Unknown", "best_score": 0, "Scores":[]}
                    if hsp.bitscore >= variant_model.hmmthreshold and hsp.bitscore > results[hit.id]["best_score"]:
                        if i==1 and not (results[hit.id]["class"] == "Unknown" or "canonical" in results[hit.id]["class"]):
                            continue
                        results[hit.id]["class"] = variant_model.id
                        results[hit.id]["best_score"] = hsp.score

                        results["scores"].append({variant_model.id: hsp.score})

    assert 0, variants

    return results


