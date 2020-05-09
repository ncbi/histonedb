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

import io
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Emboss.Applications import NeedleCommandline


Entrez.email = "shaytanak@gmail.com"


def read_gis(file):
    """
    Reads a gis from file, one per line
    """
    with open(file,'r') as f:
        gis = [re.search('(\d+)',line).group(1) for line in f if (not line.startswith('#'))]
    return gis



def get_prot_seqrec_by_gis(gi_list):
    """
    Download a dictionary of fasta SeqsRec from NCBI given a list of GIs.
    """

    print("Downloading FASTA SeqRecords by GIs from NCBI")
    num=len(gi_list)
    while True:
        fasta_seqrec=dict()
        try:
            print("Fetching %d seqs"%(num))
            strn = ",".join(gi_list)
            request=Entrez.epost(db="protein",id=strn)
            result=Entrez.read(request)
            webEnv=result["WebEnv"]
            queryKey=result["QueryKey"]
            handle=Entrez.efetch(db="protein",rettype='fasta',retmode='text',webenv=webEnv, query_key=queryKey)
            for r in SeqIO.parse(handle,'fasta'):
                fasta_seqrec[r.id.split('|')[1]]=r
        except:
            continue
        if(len(fasta_seqrec)==num):
            break
        else:
            print("Mismatch:", num," ", len(fasta_seqrec))
    print("FASTA Records downloaded:")
    print(len(fasta_seqrec))
    return(fasta_seqrec)


def get_genus_by_gi(gi):
    org=0
    while True:
        try:
            print("Fetching gi %s genus"%str(gi))
            strn = str(gi)
            # request=Entrez.epost(db="protein",id=strn)
            # result=Entrez.read(request)
            # webEnv=result["WebEnv"]
            # queryKey=result["QueryKey"]
            # handle=Entrez.efetch(db="protein",rettype='gb',retmode='text',webenv=webEnv, query_key=queryKey)
            handle=Entrez.efetch(db="protein",rettype='gb',retmode='text',id=strn)
            for r in SeqIO.parse(handle,'gb'):
                org=r.annotations["organism"]
        except:
            continue
        if(org):
            break
    return org.split()[0]


def muscle_aln(seqreclist):
    """Align with muscle"""

    muscle = os.path.join(os.path.dirname(sys.executable), "muscle")
    process = subprocess.Popen([muscle], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    sequences = "\n".join([s.format("fasta") for key,s in seqreclist.items()])
    print(sequences)
    aln, error = process.communicate(sequences)
    seqFile = io.BytesIO()
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
        # print i.description
        gi=re.search(r"gi\|(\d+)\|",i.id).group(1)
        try:
            genus=re.search(r"\[(\S+)\s+.+\S+\]",i.description).group(1)
        except:
            genus=get_genus_by_gi(gi)

        i.id=genus+"|"+gi+"|"+variant
        i.description=genus+"_"+variant+"_"+gi
        msa_r.append(i)
    return msa_r

def refactor_title_allmsa(msa):
    """
    refactors titles of sequence in format needed for histoneDB seeds
    """
    msa_r=MultipleSeqAlignment([])
    for i in msa:
        print(i.description)
        # genus=re.search(r"\[(\S+)\s+.+\S+\]",i.description).group(1)
        text=re.search(r"(\S+)\|(\d+)\|(\S+)",i.id)
        i.id=text.group(3)+"|"+text.group(1)+"|"+text.group(2)
        # i.description=genus+"_"+variant+"_"+gi
        msa_r.append(i)
    return msa_r

def get_gis(pref=''):
        """
        Goes through aux_tools/gis
        """
        for i, (root, _, files) in enumerate(os.walk("gis/"+pref)):
            hist_type = os.path.basename(root)
            for f in files:
                if not f.endswith(".gis"): continue
                yield f[:-4], hist_type, f


if __name__ == '__main__':
    # print get_genus_by_gi(223590216)
    # exit()
    if not os.path.exists("draft_seeds"):
        os.makedirs("draft_seeds")
    for hist_var,hist_type,f in get_gis():
        print("##########Starting",hist_var,hist_type,f)
        if not os.path.exists(os.path.join("draft_seeds",hist_type)):
            os.makedirs(os.path.join("draft_seeds",hist_type))
        gis=read_gis(os.path.join("gis",hist_type,f))
        seqrecs=get_prot_seqrec_by_gis(gis)
        msa=muscle_aln(seqrecs)
        print(msa)
        msa_r=refactor_title(msa,hist_var)
        msa_r.sort()
        print(msa_r)
        AlignIO.write(msa_r,os.path.join("draft_seeds",hist_type,hist_var+".fasta"),'fasta')
    #combines MSA
    for hist_type in ['H2A','H2B','H3','H4','H1']:
        seqrecs=[]
        for hist_var,hist_type,f in get_gis(hist_type):
            seqrecs+=list(SeqIO.parse("draft_seeds/"+hist_type+"/"+hist_var+".fasta", "fasta"))
        ungseqrecs={}
        for s in seqrecs:
            ungseqrecs[s.id]=SeqRecord(id=s.id,description=s.description,seq=s.seq.ungap("-"))
        # print ungseqrecs
        msa=muscle_aln(ungseqrecs)
        msa_r=refactor_title_allmsa(msa)

        msa_r.sort()
        AlignIO.write(msa_r,os.path.join("draft_seeds",hist_type+".fasta"),'fasta')



