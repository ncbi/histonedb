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
from tools.taxonomy_from_gis import taxonomy_from_gis

already_exists = []


species_re = re.compile(r'\[(.*?)\]')
def taxonomy_from_header(header, gi):
  match = species_re.findall(header)
  if match:
    organism = match[-1]
  else:
    print "No match for {}: {}, searh NCBI".format(gi, header)
    while True:
      try:
        organism = taxonomy_from_gis([gi]).next()
        break
      except StopIteration:
        pass
  organism = organism.replace(":", " ")
  organism = re.sub(r'([a-zA-Z0-9]+)[\./#_:]([a-zA-Z0-9]+)', r'\1 \2', organism)
  organism = re.sub(r'([a-zA-Z0-9]+)\. ', r'\1 ', organism)
  organism = re.sub(r"['\(\)\.]", r'', organism)
  organism = organism.lower()
  try:
    return Taxonomy.objects.get(name=organism.lower())
  except:
    try:
      genus = Taxonomy.objects.get(name=organism.split(" ")[0].lower())
      if genus.type_name != "scientific name":
        genus = genus.get_scientific_names()[0]
    except:
      print header
      return Taxonomy.objects.get(name="unidentified")

    #Maybe the var is wrong
    if "var" in organism:
      try:
        org_end = organism.split("var ")[1]
      except IndexError:
        return Taxonomy.objects.get(name="unidentified")
      best_orgs = genus.children.filter(name__endswith=org_end).all()
      if len(best_orgs):
        return best_orgs[0]
      else:
        return Taxonomy.objects.get(name="unidentified")
    else:
      return Taxonomy.objects.get(name="unidentified")
      print organism
      print "Genus =", genus
      print "Are you sure there is no taxa that matches? Check again:"
      taxas = []
      for i, c in enumerate(genus.children.all()):
        print "{}: {} ==?== {}".format(i, c, organism)
        taxas.append(c)
      index = raw_input("Choose the index that best matches: ")

      try:
        index = int(save)
      except ValueError:
        return Taxonomy.objects.get(name="unidentified")

      try:
        return taxas[index]
      except IndexError:
        return Taxonomy.objects.get(name="unidentified")


def load_hmm_results(hmmerFile, reset=True):
  """Save domain hits from a hmmer hmmsearch file into the Panchenko Histone
  Variant DB format.

  Parameters:
  ___________
  hmmerFile : string
    Path to HMMer hmmsearch output file.
  """
#We need unknown Variant model - to assign to those that do not pass the threshold for analyzed models,
#but are waiting if they will be pass threhold with other models.
#  at while searching
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
                # update_features(seq)
                best_score_2 = Score.objects.get(id=best_score.id)
                best_score_2.used_for_classification=False
                best_score_2.save()
                seq.save()
                print "UPDATED VARIANT"
                add_score(seq, variant_model, hsp, best=True)
              else:
              ##A strange thing, bitscore is bigger, but threshold is not passed.
                add_score(seq, variant_model, hsp, best=False)
            else:
            ##Did not pass threshold for previous models.
              if hsp.bitscore>=variant_model.hmmthreshold:
                seq.variant = variant_model
                seq.sequence = str(hsp.hit.seq)
                seq.save()
                # update_features(seq)
                add_score(seq, variant_model, hsp, best=True)
              else:
                add_score(seq, variant_model, hsp, best=False)
          else:
          ##A new sequence is found.
            print "!!!!!!!!", header
            print "!!!!!!!!!", gi
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
  ##Theoretically we need to get rid of all the unknown records here.
  ##But we can leave them for now, these are that were found by HMMsearch but did not pass threshold.


def add_sequence(gi, variant_model, taxonomy, header, sequence):
  if not variant_model.hist_type.id == "H1" and not variant_model.hist_type.id == "Unknown":
    hist_identified, ss_position, sequence = get_hist_ss(sequence, variant_model.hist_type.id, save_alignment=True)
    sequence = str(sequence.seq)
  else:
    ss_position = defaultdict(lambda: (None, None))
    sequence = str(sequence)
  seq = Sequence(
    id       = gi,
    variant  = variant_model,
    gene     = None,
    splice   = None,
    taxonomy = taxonomy,
    header   = header,
    sequence = sequence,
    reviewed = False,
    )
  seq.save()
  return seq

def update_features(seq, ss_position=None, variant_model=None, return_not_save=False):
  if ss_position is None and variant_model is not None:
    hist_identified, ss_position, sequence = get_hist_ss(seq.to_biopython(ungap=True).seq, variant_model.core_type.id, save_alignment=True)
    sequence = str(sequence.seq)

  # if ss_position is None:
    # assert 0

  if hasattr(seq, "features") and seq.features:
    seq.features.delete()

  if variant_model is not None:
    if not variant_model.core_type.id == "H1":
      features = Features.from_dict(seq, ss_position)
    
    if return_not_save:
      return features

    try:
      features.save()
    except:
      pass


def add_score(seq, variant_model, hsp, best=False):
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


