import re
from collections import defaultdict

#BioPython
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq

#Django libraires
from browse.models import *
from djangophylocore.models import Taxonomy
from django.db.models import Max
from django.db.utils import IntegrityError

#Custom librairies
from tools.hist_ss import get_hist_ss
from tools.taxonomy_from_gis import easytaxonomy_from_header, taxonomy_from_gis

def load_hmm_results(hmmerFile):
  """Save domain hits from a hmmer hmmsearch file into the Panchenko Histone
  Variant DB format.

  Parameters:
  ___________
  hmmerFile : string
    Path to HMMer hmmsearch output file.
  """
  #We need unknown Variant model - to assign to those that do not pass the threshold for analyzed models,
  #but are waiting if they will be pass threhold with other models.
  #at while searching
  try:
    unknown_model = Variant.objects.get(id="Unknown")
  except Variant.DoesNotExist:
    try:
      hist_unknown = Histone.objects.get(id="Unknown")
    except:
      hist_unknown = Histone("Unknown")
      hist_unknown.save()
    unknown_model = Variant(hist_type=hist_unknown, id="Unknown")
    unknown_model.save()

  for variant_query in SearchIO.parse(hmmerFile, "hmmer3-text"):
    print "Loading variant:", variant_query.id
    variant_model = Variant.objects.get(id=variant_query.id)
    for hit in variant_query:
      #Below we are fetching a list of headers if there are multiple headers for identical sequences
      #Technically HUMMER might put the second and on gis in description column.
      #The format should be strictly the genbank format: gi|343434|fafsf gdgfdg gi|65656|534535 fdafaf
      headers = "{}{}".format(hit.id, hit.description).split("gi|")[1:]

      ###Iterate through headers of identical sequences.
      for header in headers:
        gi = header.split("|")[0]
        ##Iterate through high scoring fragments.
        for i, hsp in enumerate(hit):
          seqs = Sequence.objects.filter(id=gi)
          if len(seqs)>0:
          #Sequence already exists. Compare bit scores, if loaded bit score is
          #greater than current, reassign variant and update scores. Else, append score
            seq=seqs.first()
            print "Already in database", seq
            best_scores = seq.all_model_scores.filter(used_for_classifiation=True)
            if len(best_scores)>0:
            ##Sequence have passed the threshold for one of previous models.
              best_score = best_scores.first()
              if (hsp.bitscore > best_score.score) and hsp.bitscore>=variant_model.hmmthreshold:
                #best scoring
                seq.variant = variant_model
                seq.sequence = str(hsp.hit.seq)
                best_score_2 = Score.objects.get(id=best_score.id)
                best_score_2.used_for_classification=False
                best_score_2.save()
                seq.save()
                print "UPDATED VARIANT"
                add_score(seq, variant_model, hsp, best=True)
              else:
                #A strange thing, bitscore is bigger, but threshold is not passed.
                add_score(seq, variant_model, hsp, best=False)
            else:
              #Did not pass threshold for previous models.
              if hsp.bitscore>=variant_model.hmmthreshold:
                seq.variant = variant_model
                seq.sequence = str(hsp.hit.seq)
                seq.save()
                add_score(seq, variant_model, hsp, best=True)
              else:
                add_score(seq, variant_model, hsp, best=False)
          else:
            ##A new sequence is found.
            taxonomy = taxonomy_from_header(header, gi)
            sequence = Seq(str(hsp.hit.seq))
            best = hsp.bitscore >= variant_model.hmmthreshold
            try:
              seq = add_sequence(
                gi,  
                variant_model if best else unknown_model, 
                taxonomy, 
                header, 
                sequence)
              add_score(seq, variant_model, hsp, best=best)
            except IntegrityError as e:
              print "Error adding sequence {}".format(seq)
              global already_exists
              already_exists.append(gi)
              continue
          print seq

  #Delete 'unknown' records that were found by HMMsearch but did not pass threshold.?
  #unknown_model.sequences.all().delete()

def add_sequence(gi, variant_model, taxonomy, header, sequence):
  """Add sequence into the database, autfilling empty Parameters"""
  seq = Sequence(
    id       = gi,
    variant  = variant_model,
    gene     = None,
    splice   = None,
    taxonomy = taxonomy,
    header   = header,
    sequence = str(sequence),
    reviewed = False,
    )
  seq.save()
  return seq

def add_score(seq, variant_model, hsp, best=False):
  """Add score for a given sequence"""
  score_num = Score.objects.count()+1
  score = Score(
    id                     = score_num,
    sequence               = seq,
    variant                = variant_model,
    score                  = hsp.bitscore,
    evalue                 = hsp.evalue,
    above_threshold        = hsp.bitscore >= variant_model.hmmthreshold,
    hmmStart               = hsp.query_start,
    hmmEnd                 = hsp.query_end,
    seqStart               = hsp.hit_start,
    seqEnd                 = hsp.hit_end,
    used_for_classifiation = best,
    )
  score.save()


