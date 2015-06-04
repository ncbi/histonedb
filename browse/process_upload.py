import os

import uuid
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

from browse.models import *

class TooManySequences(RuntimeError):
    pass

class InvalidFASTA(RuntimeError):
    pass

def process_upload(type, sequences, format):
    seq_file = os.path.join("media", "{}.fasta".format(uuid.uuid4()))
    if format == "text":
        with open(seq_file, 'w+') as seqs: 
            seqs.write(sequences)
        sequences = seq_file
                
    processed_sequences = list(SeqIO.parse(sequences, "fasta"))

    if len(processed_sequences) > 50:
        raise TooManySequences

    if format == "file":
        os.remove(sequences)
        SeqIO.write(processed_sequences, seq_file)

    if type == "blastp":
        result = upload_blastp(seq_file)
    elif type == "hmmer":
        result = upload_hmmer(seq_file)

    os.remove(seq_file)

    return {format: result}
"""
def upload_blastp(seq_file):
    output="media/{}.xml".format(seq_file[:-6])
    blastx_cline = NcbiblastxCommandline(
        query=seq_file, db="static/browse/blast/HistoneDB.db", evalue=0.01,
        outfmt=5, out=output)
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
    variant_output= os.path.join("media", "{}_variant.out".format(seq_file[:-6]))
    core_output= os.path.join("media", "{}_core.out".format(seq_file[:-6]))
    variantdb = os.path("static", "browse", "hmms", "combined_variants.hmm")
    coredb = os.path("static", "browse", "hmms", "combined_cores.hmm")

    subprocess.call(["hmmsearch", "-o", variant_output, "-E", str(evalue), "--cpu", "4", "--notextw", variantdb, seq_file])
    subprocess.call(["hmmsearch", "-o", core_output, "-E", str(evalue), "--cpu", "4", "--notextw", variantdb, seq_file])

    results = {}
    for i, search_type in enumerate((variant_output, core_output)):
        if i==1: 
            variant = "canonical{}".format(core_query.id)
        else:
            variant = core_query.id
        for query in SearchIO.parse(variant_output, "hmmer3-text"):
            try:
                variant_model = Variant.objects.get(id=name)
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

    os.remove(variant_output)
    os.remove(core_output)

    return results
"""

