# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
This script determines secondary structure elements of histone,
by algining it to reference sequences.

"""
#Standard Library
import argparse
import csv
import collections
import os
import re
import subprocess
import sys
import uuid
import cPickle as pickle
from StringIO import StringIO
import logging

#Required libraires
#import pylab
#import pandas as pd
#import networkx as nx

#Django
from django.conf import settings
from browse.models import Feature, Sequence

#BioPython
from Bio import AlignIO
from Bio import Entrez
from Bio import ExPASy
from Bio import SeqIO
from Bio import SwissProt
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

Entrez.email = "l.singh@intbio.org"

# Logging info
log = logging.getLogger(__name__)

def get_variant_features(sequence, variants=None, save_dir="", save_not_found=False, save_gff=True):
    """Get the features of a sequence based on its variant.

    Parameters:
    -----------
    sequence : Sequence django model
        The seuqence to add get features for with identified variant
    variants : List of Variant models
        Anntate others variants. Optional.
    save_dir : str
        Path to save temp files.
    save_not_found : bool
        Add Features even if they weren't found. Indices will be (-1, -1)

    Return:
    -------
    A string containing the gff file of all features
    """
    #Save query fasta to a file to EMBOSS needle can read it
    n2=str(uuid.uuid4())
    test_record = sequence.to_biopython()
    query_file = os.path.join(save_dir, "query_{}.fa".format(n2))
    SeqIO.write(test_record, query_file, 'fasta')

    #A list of updated Features for the query
    variant_features = set()

    if not variants:
        variants = [sequence.variant]

    for variant in variants:
        for template_variant in [variant.id, "General{}".format(variant.hist_type.id)]:
            try:
                features = Feature.objects.filter(template__variant=template_variant)
            except:
                continue
            #Find features with the same taxonomy
            tax_features = features.filter(template__taxonomy=sequence.taxonomy)
            if len(tax_features) == 0:
                #Find features with closest taxonomy => rank class
                tax_features = features.filter(template__taxonomy__parent__parent__parent=sequence.taxonomy.parent.parent.parent)
            if len(tax_features) == 0:
                #Nothing, use unidentified which is the standard
                tax_features = features.filter(template__taxonomy__name="undefined")
            features = tax_features
            for updated_feature in transfer_features_from_template_to_query(features, query_file, save_dir=save_dir, save_not_found=save_not_found):
                variant_features.add(updated_feature)
   
    os.remove(query_file)
    if save_gff:
        return Feature.objects.gff(sequence.id, variant_features)
    return variant_features

def transfer_features_from_template_to_query(template_features, query_file, save_dir="", save_not_found=False):
    """Transfer features from template to query. Position are defined in the 
    template and we use needle to find the corresponding position in the template

    Parameters:
    -----------
    template_features : QuerySet of Feature django models
        The features that relate to the template. 
    query_file : str
        Path to FASTA file containing query sequence 
    save_dir : str
        Path to save temp files.
    save_not_found : bool
        Add Features even if they weren't found. Indices will be (-1, -1)

    Yeilds:
    -------
    A Feature django model with the name of the feature and position relative to the query
    """
    if len(template_features) == 0:
        return

    n2=str(uuid.uuid4())
    template = template_features.first().template
    template_file = template.path()
    needle_results = os.path.join(save_dir, "needle_{}.txt".format(n2))
    cmd = os.path.join(os.path.dirname(sys.executable), "needle")

    if not os.path.isfile(cmd):
        cmd = "needle"
    needle_cline = NeedleCommandline(
        cmd=cmd,
        asequence=template_file, 
        bsequence=query_file, 
        gapopen=10,
        gapextend=1, 
        outfile=needle_results)
    stdout, stderr = needle_cline()

    align = AlignIO.read(needle_results, "emboss")
    # print align.format("fasta")
    core_histone = align[0]
    query = align[1]
    
    corresponding_hist = range(len(template.get_sequence()))
    k=0
    for i, core_histone_postion in enumerate(core_histone):
        if core_histone_postion == "-":
            k += 1
        else:
            corresponding_hist[i-k]=i


    corresponding_test = range(len(SeqIO.parse(query_file, "fasta").next()))
    k=0
    for i, query_position in enumerate(query):
        if query_position == "-":
            k=k+1
        else:
            corresponding_test[i-k]=i


    for feature in template_features:
        start = feature.start
        stop = feature.end
        start_in_aln = corresponding_hist[start]
        end_in_aln = corresponding_hist[stop]
        start_in_test_seq = -1
        end_in_test_seq = -1

        for k in xrange(len(core_histone)):
            try:
                start_in_test_seq = corresponding_test.index(start_in_aln+k)
                break
            except ValueError:
                continue

        for k in xrange(len(core_histone)):
            try:
                end_in_test_seq = corresponding_test.index(end_in_aln-k)
                break
            except ValueError:
                continue

        if start_in_test_seq == -1 or end_in_test_seq == -1 or start_in_test_seq > end_in_test_seq:
            if save_not_found:
                yield Feature(
                    id          = "{}_{}".format(os.path.splitext(query_file)[0], feature.id),
                    name        = feature.name,
                    description = feature.description,
                    start       = -1,
                    end         = -1,
                    color       = feature.color,
                )
        else:
            yield Feature(
                id          = "{}_{}".format(os.path.splitext(query_file)[0], feature.id),
                name        = feature.name,
                description = feature.description,
                start       = start_in_test_seq,
                end         = end_in_test_seq,
                color       = feature.color,
            )

    #Cleanup
    os.remove(needle_results)

def get_features_in_aln(alignment, variant, save_dir="", save_gff=True):
    #Let's extract consensus
    a=SummaryInfo(alignment)
    cons=a.dumb_consensus(threshold=0.1, ambiguous='X')
    seq = Sequence(id="Consensus", variant_id=variant, taxonomy_id=1, sequence=str(cons))
    updated_features = get_variant_features(seq, save_dir=save_dir, save_gff=save_gff)
    return updated_features

def get_core_lendiff(query_features, template_features):
    """Get the ratio of core length for query sequence versus template sequence

    Parameters:
    -----------
    query_features : Iterable of Feature objects (QuerySet or list)
        Features with positions relative to query sequence
    template_features : Iterable of Feature objects (QuerySet or list)
        Features with positions relative to template sequence

    Return:
    -------
    ratio : float
        Ratio of core lengths
    """
    #check 640798122
    len_t_core=max([ f.end for f in template_features ])-min([ f.start for f in template_features ])
    len_core=max([ f.end for f in query_features ])-min([ f.end for f in query_features ])
    if(debug):
        log.info("Template core length {}".format(len_t_core))
        log.info("Testseq core length {}".format(len_core))
    ratio=float(len_core)/float(len_t_core)
    return ratio
