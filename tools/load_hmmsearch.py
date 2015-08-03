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

def get_sequence(gi, sequence_file):
  """sequences is likely large, so we don't want to store it in memory.
  Read file in each time, saving only sequence with gi"""
  for s in SeqIO.parse(sequence_file, "fasta"):
    if gi in s.id:
      return s

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


def load_variants(hmmerFile, sequences, reset=True):
  """Save domain hits from a hmmer hmmsearch file into the Panchenko Histone
  Variant DB format.

  Parameters:
  ___________
  hmmerFile : string
    Path to HMMer hmmsearch output file
  sequences : string
    Path to of sequences used to search the HMM
  threshold : float
    Keep HSPS with scores >= threshold. Optional.
  """

  if reset:
    Sequence.objects.all().delete()
    Features.objects.all().delete()
    Score.objects.all().delete()

  try:
    unknown_model = Variant.objects.get(id="Unknown")
  except Variant.DoesNotExist:
    try:
      core_unknown = Histone.objects.get(id="Unknown")
    except:
      core_unknown = Histone("Unknown")
      core_unknown.save()
    unknown_model = Variant(core_type=core_unknown, id="Unknown")
    unknown_model.save()
  
  for variant_query in SearchIO.parse(hmmerFile, "hmmer3-text"):
    print "Loading variant:", variant_query.id
    try:
      variant_model = Variant.objects.get(id=variant_query.id)
    except:
      if "H2A" in variant_query.id:
        core_histone = Histone.objects.get(id="H2A")
      elif "H2B" in variant_query.id:
        core_histone = Histone.objects.get(id="H2B")
      elif "H3" in variant_query.id:
        core_histone = Histone.objects.get(id="H3")
      elif "H4" in variant_query.id:
        core_histone = Histone.objects.get(id="H4")
      elif "H1" in variant_query.id:
        core_histone = Histone.objects.get(id="H1")
      else:
        continue
      variant_model = Variant(id=variant_query.id, core_type=core_histone)
      variant_model.save()
    for hit in variant_query:
      #here is the problem, following line works only for NR formatted seqs, not for our seeds. we will try to replace it with a comprehensive version
      old_headers = "{}{}".format(hit.id, hit.description).split("gi|")[1:]
      #88888|ref|XP_dfdf|description - this is what we should get
      tmp = "{}{}".format(hit.id, hit.description).split("|")
      if(tmp[0]=="gi"):
        headers=["|".join(tmp[1:])]
      else:
        if(re.match("\d+",tmp[0])):
          headers=["|".join(tmp)]
        else:
          if(re.match("\d+",tmp[1])):
            headers=["|".join([tmp[1],tmp[0]]+tmp[2:])]
      # print headers
      # print old_headers
      # assert headers==old_headers

      for header in headers:
        gi = header.split("|")[0]
        for i, hsp in enumerate(hit):
          seqs = Sequence.objects.filter(id=gi)
          if len(seqs) > 0: 
            #Sequence already exists. Compare bit scores, if current bit score is 
            #greater than current, reassign variant and update scores. Else, append score
            seq = seqs.first()
            print "Repeated", seq
            best_scores = seq.all_model_scores.filter(used_for_classifiation=True)
            if len(best_scores)>0:
              best_score = best_scores.first()
              if (hsp.bitscore > best_score.score) and hsp.bitscore>=variant_model.hmmthreshold:
                #best scoring
                seq.variant = variant_model
                seq.sequence = str(hsp.hit.seq)
                update_features(seq)
                best_score_2 = Score.objects.get(id=best_score.id)
                best_score_2.used_for_classification=False
                best_score_2.save()
                seq.save()
                print "UPDATED VARIANT"
                add_score(seq, variant_model, hsp, best=True)
              else:
                add_score(seq, variant_model, hsp, best=False)
            else:
              #Was not classified, just test against this model
              if hsp.bitscore>=variant_model.hmmthreshold:
                seq.variant = variant_model
                seq.sequence = str(hsp.hit.seq)
                seq.save()
                update_features(seq)
                add_score(seq, variant_model, hsp, best=True)
              else:
                add_score(seq, variant_model, hsp, best=False)
          else:
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

def load_cores(hmmerFile, reset=True):
  unknown_model = Variant.objects.get(id="Unknown")
  
  if reset:
    variants = Variant.objects.filter(id__contains="cononical")
    for variant in variants:
      for sequence in variant.sequences:
        sequence.features.delete()
        sequence.scores.filter(variant__id=variant__id).delete()
      variant.sequences.delete()
    variants.delete()

  for core_query in SearchIO.parse(hmmerFile, "hmmer3-text"):
    print "Loading core:", core_query.id
    try:
      core_histone = Histone.objects.get(id=core_query.id)
    except:
      continue

    try:
      canonical_model = Variant.objects.get(id="canonical{}".format(core_query.id))
    except:
      canonical_model = Variant(id="canonical{}".format(core_query.id), core_type=core_histone)
      canonical_model.save()

    for hit in core_query:
      headers = "{}{}".format(hit.id, hit.description).split("gi|")[1:]
      for header in headers:
        gi = header.split("|")[0]
        for i, hsp in enumerate(hit):
          seqs = Sequence.objects.filter(id=gi)
          if len(seqs) > 0: 
            #Sequence already exists. Compare bit scores, if current bit score is 
            #greater than current, reassign variant and update scores. Else, append score
            seq = seqs.first()
            best_scores = seq.all_model_scores.filter(used_for_classifiation=True)
            if len(best_scores)>0:
              best_score = best_scores.first()
              if hsp.bitscore >= canonical_model.hmmthreshold and \
                  (seq.variant.id == "Unknown" or  \
                    ("canonical" in seq.variant.id and hsp.bitscore > best_score.score)) :
                seq.variant = canonical_model
                seq.sequence = str(hsp.hit.seq)
                update_features(seq)
                best_score_2 = Score.objects.get(id=best_score.id)
                best_score_2.used_for_classification=False
                best_score_2.save()
                seq.save()
                print "UPDATED VARIANT"
                add_score(seq, canonical_model, hsp, best=True)
              else:
                add_score(seq, canonical_model, hsp, best=False)
            else:
              #Was not classified, just test against this model
              if hsp.bitscore>=canonical_model.hmmthreshold:
                seq.variant = canonical_model
                seq.sequence = str(hsp.hit.seq)
                seq.save()
                update_features(seq)
                add_score(seq, canonical_model, hsp, best=True)
              else:
                add_score(seq, canonical_model, hsp, best=False)
          else:
            #only add sequence if it was not found by the variant models
            taxonomy = taxonomy_from_header(header, gi)
            sequence = Seq(str(hsp.hit.seq))
            best = hsp.bitscore >= canonical_model.hmmthreshold
            try:
              seq = add_sequence(
                gi,  
                canonical_model if best else unknown_model, 
                taxonomy, 
                header, 
                sequence)
              add_score(seq, canonical_model, hsp, best=best)
            except IntegrityError as e:
              print "Error adding sequence {}".format(seq)
              global already_exists
              already_exists.append(gi)
              continue
          print seq
  print already_exists

def add_sequence(gi, variant_model, taxonomy, header, sequence):
  if not variant_model.core_type.id == "H1" and not variant_model.core_type.id == "Unknown":
    hist_identified, ss_position, sequence = get_hist_ss(sequence, variant_model.core_type.id, save_alignment=True)
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
