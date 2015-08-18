#!/usr/bin/env python
"""Gets a file with GIs, obtains FASTA from NCBI and aligns with MUSCLE.
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


def read_gis(file):
    """
    Reads a gis from file, one per line
    """
    with open(file,'r') as f:
        gis = [line.rstrip() for line in f]
    return gis



def get_prot_seqrec_by_gis(gi_list):
    """
    Download a dictionary of fasta SeqsRec from NCBI given a list of GIs.
    """

    print("Downloading FASTA SeqRecords by GIs from NCBI")
    num=len(gi_list)
    fasta_seqrec=dict()
    for i in range(int(num/1000)+1):
        while True:
            try:
                print("Fetching %d th thousands from %d"%(i,num))
                strn = ",".join(gi_list[i*1000:(i+1)*1000])
                request=Entrez.epost(db="protein",id=strn)
                result=Entrez.read(request)
                webEnv=result["WebEnv"]
                queryKey=result["QueryKey"]
                handle=Entrez.efetch(db="protein",rettype='fasta',retmode='text',webenv=webEnv, query_key=queryKey)
                for r in SeqIO.parse(handle,'fasta'):
                    fasta_seqrec[r.id.split('|')[1]]=r
            except:
                    continue
            break
    print("FASTA Records downloaded:")
    print(len(fasta_seqrec))
    return(fasta_seqrec)




def muscle_aln(seqreclist):
    """Align with muscle"""

    muscle = os.path.join(os.path.dirname(sys.executable), "muscle")
    process = subprocess.Popen([muscle], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    sequences = "\n".join([s.format("fasta") for key,s in seqreclist.iteritems()])
    print sequences
    aln, error = process.communicate(sequences)
    seqFile = StringIO.StringIO()
    seqFile.write(aln)
    seqFile.seek(0)
    sequences = list(SeqIO.parse(seqFile, "fasta")) #Not in same order, but does it matter?
    msa = MultipleSeqAlignment(sequences)
    return msa

def refactor_title(msa,variant):
    """
    refactors titles of sequence in format needed for histoneDB seeds
    """
    msa_r=MultipleSeqAlignment([])
    for i in msa:
        genus=re.search(r"\[(\S+)\s+\S+\]",i.description).group(1)
        gi=re.search(r"gi\|(\d+)\|",i.id).group(1)
        i.id=genus+"|"+gi+"|"+variant
        i.description=genus+"_"+variant+"_"+gi
        msa_r.append(i)
    return msa_r


if __name__ == '__main__':
    
    gis=read_gis('canonicalH2A.gis')
    seqrecs=get_prot_seqrec_by_gis(gis)
    msa=muscle_aln(seqrecs)
    print msa
    msa_r=refactor_title(msa,'canonicalH2A')
    print msa_r
    AlignIO.write(msa_r,"H2A_aln.fasta",'fasta')
