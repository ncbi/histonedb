import re

#BioPython
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq

#Django libraires
from browse.models import *
from djangophylocore.models import Taxonomy
from django.db.models import Max

#Custom librairies
from tools.hist_ss import get_hist_ss
from tools.taxonomy_from_gis import taxonomy_from_gis

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
      org_end = organism.split("var ")[1]
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


def parseHmmer(hmmerFile, sequences, threshold=10.):
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
  Sequence.objects.all().delete()
  Features.objects.all().delete()
  Score.objects.all().delete()
  for variant_query in SearchIO.parse(hmmerFile, "hmmer3-text"):
    print "Loading variant:", variant_query.id
    variant_model = Variant.objects.get(id=variant_query.id)
    for hit in variant_query:
      headers = "{}{}".format(hit.id, hit.description).split("gi|")[1:]
      for header in headers:
        gi = header.split("|")[0]
        taxonomy = taxonomy_from_header(header, gi)
        #print "  Loading {} from {}".format(gi, taxonomy)
        splice_varaint = None
        for i, hsp in enumerate(hit):
          if hsp.bitscore < threshold: continue
          seqs = Sequence.objects.filter(id=gi)
          if len(seqs) == 1:
            #Sequence already exists. Compare bit scores, if current bit score is 
            #greater than current, reassign variant and update scores. Else, append score
            best_score = seqs.aggregate(score=Max("scores"))["score"]
            if hsp.bitscore > best_score:
              #best scoring
              seq.variant = variant_model
            """
            is_best_score = False
            elif hsp.bitscore == best_score.score and (hsp.hit_start-hsp.hitEnd+1)<(best_score.seqTo-best_score.seqFrom+1):
              #Edge case: if the scores are equal, choose the sequence that had the longest match in HMM
              is_best_score = False
            else:"""
          else:
            sequence = Seq(str(hsp.hit.seq))
            hist_identified, ss_position, sequence = get_hist_ss(sequence, variant_model.core_type.id, save_alignment=True)
            seq = Sequence(
              id       = gi,
              variant  = variant_model,
              gene     = None,
              splice   = None,
              taxonomy = taxonomy,
              header   = header,
              sequence = str(sequence.seq),
              reviewed = False,
              )
            seq.save()

            features = Features(
              sequence             = seq,
              alphaN_start         = ss_position["alphaN"][0],
              alphaN_end           = ss_position["alphaN"][1],
              alpha1_start         = ss_position["alpha1"][0],
              alpha1_end           = ss_position["alpha1"][1],
              alpha1ext_start      = ss_position["alpha1ext"][0],
              alpha1ext_end        = ss_position["alpha1ext"][1],
              alpha2_start         = ss_position["alpha2"][0],
              alpha2_end           = ss_position["alpha2"][1],
              alpha3_start         = ss_position["alpha3"][0],
              alpha3_end           = ss_position["alpha3"][1],
              alpha3ext_start      = ss_position["alpha3ext"][0],
              alpha3ext_end        = ss_position["alpha3ext"][1],
              alphaC_start         = ss_position["alphaC"][0],
              alphaC_end           = ss_position["alphaC"][1],
              beta1_start          = ss_position["beta1"][0],
              beta1_end            = ss_position["beta1"][1],
              beta2_start          = ss_position["beta2"][0],
              beta2_end            = ss_position["beta2"][1],
              loopL1_start         = ss_position["loopL1"][0],
              loopL1_end           = ss_position["loopL1"][1],
              loopL2_start         = ss_position["loopL2"][0],
              loopL2_end           = ss_position["loopL2"][1],
              mgarg1_start         = ss_position["mgarg1"][0],
              mgarg1_end           = ss_position["mgarg1"][1],
              mgarg2_start         = ss_position["mgarg2"][0],
              mgarg2_end           = ss_position["mgarg2"][1],
              mgarg3_start         = ss_position["mgarg3"][0],
              mgarg3_end           = ss_position["mgarg3"][1],
              docking_domain_start = ss_position["docking domain"][0],
              docking_domain_end   = ss_position["docking domain"][1],
              core                 = ss_position["core"],
              )
            features.save()
          
          score = Score(
            sequence       = seq,
            variant        = variant_model,
            score          = hsp.bitscore,
            evalue         = hsp.evalue,
            train_program  = "{}{}".format(variant_query.program, variant_query.version),
            search_program = "{}{}".format(variant_query.program, variant_query.version),
            specificity     = None,
            hmmStart       = hsp.query_start,
            hmmEnd         = hsp.query_end,
            seqStart       = hsp.hit_start,
            seqEnd         = hsp.hit_end
            )
          score.save()
