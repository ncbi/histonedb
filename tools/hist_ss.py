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

#Required libraires
import pylab
import pandas as pd
import networkx as nx

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

Entrez.email = "eli.draizen@nih.gov" 

def get_hist_ss(test_seq,hist_type='Unknown',debug=0,save_alignment=False):
    """Returns sequence elements in histone sequence, all numbers assume first element in seq has number 0!!! Not like in PDB"""

    #Let's define 1kx5 sequences
    templ_H3 = Seq("ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVALFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA", IUPAC.protein)
    templ_H4 = Seq("SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG", IUPAC.protein)
    templ_H2A = Seq("SGRGKQGGKTRAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLPKKTESSKSKSK", IUPAC.protein)
    templ_H2B = Seq("AKSAPAPKKGSKKAVTKTQKKDGKKRRKTRKESYAIYVYKVLKQVHPDTGISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK", IUPAC.protein)

    #'element_name':[start,stop], start stop - are inclusive as in PDB file
    #Numbering differes between symmetrical chains and 1kx5 vs 1aoi.
    #We simply take the minimum length of alpha helices over all chains in 1kx5
    #1 substructed from PDB values!!! because these values are in array index numberins starting from 0
    ss_templ_H3={'alphaN':[43,56],'alpha1':[62,76],'alpha2':[84,113],'alpha3':[119,130],'loopL1':[78,83],'loopL2':[114,118],'beta1':[82,83],'beta2':[117,118],'mgarg1':[62,62],'mgarg2':[82,82],'mgarg3':[48,48]}
    ss_templ_H4={'alpha1ext':[23,28],'alpha1':[29,40],'alpha2':[48,75],'alpha3':[81,92],'loopL1':[41,47],'loopL2':[76,81],'beta1':[44,45],'beta2':[79,80],'beta3':[95,97],'mgarg1':[44,44]}
    # ss_templ_H2A={'alpha1ext':[15,21],'alpha1':[25,36],'alpha2':[45,72],'alpha3':[78,88],'alpha3ext':[89,96],'loopL1':[37,44],'loopL2':[73,77],'beta1':[41,42],'beta2':[76,77],'beta3':[99,101],'docking domain':[91,107],'docking tail':[108,116],'mgarg1':[41,41],'mgarg2':[76,76]}
    #new def of docking domains as in Suto Luger 2000
    ss_templ_H2A={'alpha1ext':[15,21],'alpha1':[25,36],'alpha2':[45,72],'alpha3':[78,88],'alpha3ext':[89,96],'loopL1':[37,44],'loopL2':[73,77],'beta1':[41,42],'beta2':[76,77],'beta3':[99,101],'docking domain':[80,118],'mgarg1':[41,41],'mgarg2':[76,76]}
    
    ss_templ_H2B={'alpha1':[33,45],'alpha2':[51,80],'alpha3':[86,98],'alphaC':[99,119],'loopL1':[46,50],'loopL2':[81,85],'beta1':[49,50],'beta2':[84,85],'mgarg1':[29,29]}

    ss_templ={'H3':ss_templ_H3,'H4':ss_templ_H4,'H2A':ss_templ_H2A,'H2B':ss_templ_H2B}
    templ={'H3':templ_H3,'H4':templ_H4,'H2A':templ_H2A,'H2B':templ_H2B}

    #Lets use blast and see what histone is our query
    my_records=[SeqRecord(templ_H3,id='H3',name='H3'),SeqRecord(templ_H4,id='H4',name='H4'),SeqRecord(templ_H2A,id='H2A',name='H2A'),SeqRecord(templ_H2B,id='H2B',name='H2B')]
    
    n1=str(uuid.uuid4())
    n2=str(uuid.uuid4())
    test_record = SeqRecord(test_seq,id='Query',name='Query')
    SeqIO.write([test_record],n2+'.fasta','fasta')

    if(hist_type=='Unknown'):

            
        SeqIO.write(my_records, n1+".faa", "fasta")
        os.system('makeblastdb -dbtype prot -in %s.faa -out %s.db > /dev/null'%(n1,n1))



        blastp_cline = NcbiblastpCommandline(query=n2+".fasta", db=n1+".db", evalue=100,outfmt=5, out=n1+".xml")
        stdout, stderr = blastp_cline()

        blast_record = NCBIXML.read(open(n1+'.xml','r'))

        sname=list()
        evalue=list()
        hsp_list=list()
        # length_list=list()
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sname.append(alignment.title)
                evalue.append(hsp.expect)
                hsp_list.append(hsp)
                # length_list.append(alignment.length)
        hist_identified=sname[evalue.index(min(evalue))].split()[1]
        hsp=hsp_list[evalue.index(min(evalue))]
        # length=length_list[evalue.index(min(evalue))]
    else:
        hist_identified=hist_type

    if(debug): print('Most likely this is histone:')
    if(debug): print(hist_identified)
    if(debug): print('Blast alignment')
    #We need to determine secondary strucutre according to template using the alignment
    if(debug): print(hsp)

    SeqIO.write([SeqRecord(templ[hist_identified],id=hist_identified,name=hist_identified)],n1+'.fasta','fasta')

    #Now we will redo it with Needlman Wunsh - the global alignment
    needle_cline = NeedleCommandline(asequence=n1+".fasta", bsequence=n2+".fasta",gapopen=20, gapextend=1, outfile=n1+".txt")
    stdout, stderr = needle_cline()
    # print('Needle alignment')

    align = AlignIO.read(n1+".txt", "emboss")
    if(debug): print align
    # print(hsp.gaps)
    #Blast checking
    # ss_test=dict()
    # for key,value in ss_templ[hist_identified].iteritems():
    #     print('Checking %s'%key)
    #     if((hsp.sbjct_start<=value[1])&((hsp.sbjct_end)>=value[0])):
    #         print('Belongs')
    #     else:
    #         print('Not')
    
    #Now we will get correspondence

    ss_test=dict()
    hist=templ[hist_identified]

    corrsp_hist=range(len(hist))
    k=0
    for a,i in zip(align[0],range(len(align[0]))):
        if(a=='-'):
            k=k+1
        else:
            corrsp_hist[i-k]=i
    if(debug): print corrsp_hist


    corrsp_test=range(len(test_seq))
    k=0
    for a,i in zip(align[1],range(len(align[1]))):
        if(a=='-'):
            k=k+1
        else:
            corrsp_test[i-k]=i
    if(debug): print corrsp_test


    for key,value in ss_templ[hist_identified].iteritems():
        if(debug): print('Checking %s'%key)
        start_in_aln=corrsp_hist[value[0]]
        if(debug): print('Start in aln %d'%start_in_aln)

        end_in_aln=corrsp_hist[value[1]]
        if(debug): print('End in aln %d'%end_in_aln)

        for k in range(len(align[0])):
            try:
                start_in_test_seq=corrsp_test.index(start_in_aln+k)
            except:
                start_in_test_seq=-1
                if(debug): print("Trying to move start"),
                continue
            break
        # print('\n %d'%start_in_test_seq)

        for k in range(len(align[0])):
            try:
                end_in_test_seq=corrsp_test.index(end_in_aln-k)
            except:
                end_in_test_seq=-1
                if(debug): print('Trying to move end'),
                continue
            break

        # print('\n %d'%end_in_test_seq)
        if((start_in_test_seq==-1)|(end_in_test_seq==-1)|(start_in_test_seq>end_in_test_seq)):
            ss_test[key]=[-1,-1]
        else:
            ss_test[key]=[start_in_test_seq,end_in_test_seq]
        if(debug): print ss_test[key]

    ss_test["core"] = get_core_lendiff(ss_test, ss_templ[hist_identified])

    if(hist_type=='Unknown'):
        os.system("rm %s.faa %s.db.phr %s.db.pin %s.db.psq %s.fasta %s.xml %s.txt %s.fasta"%(n1,n1,n1,n1,n2,n1,n1,n1))
    else:
        os.system("rm   %s.fasta  %s.txt %s.fasta"%(n2,n1,n1))
        
    if save_alignment:
        return hist_identified,ss_test,align[1]

    return hist_identified,ss_test

def get_hist_ss_in_aln(alignment,hist_type='Unknown',debug=0):
    """Returns sequence elements in histone alignment, all numbers assume first element in seq has number 0!!! Not like in PDB"""

    #Let's extract consensus
    if(debug):
        print alignment
    a=SummaryInfo(alignment)
    cons=a.dumb_consensus(threshold=0.1, ambiguous='X')
    if(debug):
        print "Consensus"
        print cons
    hv,ss=get_hist_ss(cons,hist_type,debug)
    return hv,ss


def get_core_lendiff(test_ss,temp_ss,debug=0):
    """Returns ration of core length for test_seq versus template sequence"""
    #check 640798122
    len_t_core=max([ i[1] for i in temp_ss.values() ])-min([ i[0] for i in temp_ss.values() ])
    len_core=max([ i[1] for i in test_ss.values() ])-min([ i[0] for i in test_ss.values() ])
    if(debug):
        print "Template core length ", len_t_core
        print "Testseq core length ", len_core
    ratio=float(len_core)/float(len_t_core)
    return ratio

if __name__ == '__main__':

    H2A = Seq("SGRGKKKKKQGGKTRAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPPPPPVYLAAVLEYLTAEILELAGNARRRRARDNKTTTTTTKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLPKKTESSKSKSK", IUPAC.protein)
    H2At = Seq("SGRGKQGGKTRAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLPKKTESSKSKSK", IUPAC.protein)
    print get_hist_ss(H2At)
    # print get_core_lendiff(H2A,H2At)