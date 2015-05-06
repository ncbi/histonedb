#BioPython
from Bio import SearchIO
from Bio import SeqIO

#Django libraires
from server.models import *
from phylocore_models import *

#Custom librairies
from hist_ss import get_hist_ss
from taxonomy_from_gis import taxonomy_from_header

def get_sequence(gi, sequence_file):
  """sequences is likely large, so we don't want to store it in memory.
  Read file in each time, saving only sequence with gi"""
  for s in SeqIO.parse(sequence_file, "fasta"):
    if gi in s.decription:
      return s


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
  for variant_query in SearchIO.parse(hmmerFile, "hmmer3-text"):
    variant_model = Sequences.objects.get(id=variant_query.id)
    for hit in variant_query:
      gi = hit.id.split("|")[1]
      splice_varaint = None
      header = hit.description
      for i, hsp in enumerate(hit):
        if hsp.bitscore < threshold: continue
        is_best_score = True
        hist_identified, ss_test, sequence = hist_ss(get_sequence(gi, sequences))
        try:
          #Sequence already exists. Compare bit scores, if current bit score is 
          #greater than current, reassign variant and update scores. Else, append score
          seq = Sequences.objects.get(id=gi)
          best_score = seq.scores.objects.filter(best=True)
          if hsp.bitscore < best_score.score:
            is_best_score = False
          elif hsp.bitscore == best_score.score and (hsp.hit_start-hsp.hitEnd+1)<(best_score.seqTo-best_score.seqFrom+1):
            #Edge case: if the scores are equal, choose the sequence that had the longest match in HMM
            is_best_score = False
          else:
            #best scoring
            seq.variant = variant_model
            best_score.best = False #Reset best after the exception
        except Entry.DoesNotExist:
          species_name = taxonomy_from_header(header)
          taxonomy = Taxonomy.objects.get(name=species_name)
          seq = Sequence(
            id       = gi,
            variant  = variant_model,
            gene     = None,
            splice   = None,
            reviewed = False,
            taxonomy = taxonomy,
            sequence = sequence
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
          id       = seq,
          score    = hsp.score,
          evalue   = hsp.evalue,
          program  = "{}{}".format(query.program, query.version),
          best     = is_best_score
          hmmStart = hsp.query_start,
          hmmEnd   = hsp.query_end,
          seqStart = hsp.hit_from,
          seqEnd   = hsp.hit_end
          )
        score.save()

