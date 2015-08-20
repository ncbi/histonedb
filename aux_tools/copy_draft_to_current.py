#!/usr/bin/env python
"""
Show differences between draft and current seeds
"""

from math import sqrt, log, e
from random import choice, random

import uuid
from Bio import ExPASy
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import csv
import collections
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess

import StringIO
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Emboss.Applications import NeedleCommandline


Entrez.email = "alexey.shaytan@nih.gov"


def get_draft_seeds():
        """
        Goes through aux_tools/gis
        """
        for i, (root, _, files) in enumerate(os.walk("draft_seeds/")):
            hist_type = os.path.basename(root)
            for f in files:
                if not f.endswith(".fasta"): continue
                yield root, f


if __name__ == '__main__':
    for root,f in get_draft_seeds():
        # print hist_var,hist_type,f
        seedpath=os.path.join("../static/browse/",root,f).replace("draft_","")
        print "Copying ",os.path.join(root,f)," to ",seedpath

        if not os.path.exists(seedpath):
            print "seed "+f+" does not exist"
        os.system("cp "+os.path.join(root,f)+" "+seedpath)