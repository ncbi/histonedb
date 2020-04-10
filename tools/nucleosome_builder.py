"""This script homology models variant nucleosomes based on 1kx5 nucleosome 
structure

An example of this can be used to model a yeast nucleosome:
1) H2A: yeast HTA1 gene (gi 398366187)
2) H2B: HTB1 (gi 398366183)
3) H3: CSE4 (gi 27808712)
4) H4: HHF1 (gi 6319481)
5) DNA: A 601 nt DNA sequence, with truncated tails according to our NCP model simulations (see Shaytan et al. NAR(?) 2015)
   The DNA is also following the sequence
   aaGTCACATGATGATATTTGATTTTATTATATTTTTAAAAAAAGTAAAAAATAAAAAGTAG T TTATTTTTAAAAAATAAAATTTAAAATATTAGTGTATTTGATTTCCGAAAGTTAAAAaaga

Origianlly written by Alexey Shaytan
Modified by Eli Draizen to be reusable
"""

import os
from chimera import runCommand as rc  # use 'rc' as shorthand for runCommand
from chimera import replyobj  # for emitting status messages
from chimera import openModels
import sys

sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/ALIGNMENT_TOOLS')
sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/X3DNA_analysis_modeling')

sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/structure_analysis')
sys.path.append('/Library/modeller-9.14/modlib')
sys.path.append('/Library/modeller-9.14/lib/mac10v4')
# sys.path.append(';'.join(['', '/Users/alexeyshaytan/Library/Python/2.7/lib/python/site-packages/setuptools-0.9.6-py2.7.egg', '/Users/alexeyshaytan/Library/Python/2.7/lib/python/site-packages/PROPKA-3.1-py2.7.egg', '/Library/modeller-9.14/examples/atom_files', '/Applications/Bioinf/biana', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/readline', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/Users/alexeyshaytan/Library/Python/2.7/lib/python/site-packages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/PIL', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/PyObjC', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/gtk-2.0']))
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import L_fasta2pir
from L_aln_tools import muscle_aln

# from dna_tools_simple import change_dna_seq_in_pdb

os.environ['LD_LIBRARY_PATH'] = '/Library/modeller-9.14/lib/mac10v4'

from modeller import *  # Load standard Modeller classes
from modeller.automodel import *  # Load the automodel class


def create_nucleosome(outfile, save_chain=None, **chains):
    """Homology model a complete nucleosome particle

    Parameters:
    -----------
    outfile : File-like object
        Where to write the output.
    save_chain : str
        Save only the strucutre from chain with name save_chain.
    chains: key, value pairs
        key is chain name, value is sequence to model

    Return:
    -------
    A file is written to outfile
    """
    chain_names = "ABCDEFGH"
    template = {sequence.id.split("_")[1]: sequence for sequence in SeqIO.parse()}

    seq_aln = {}
    for chain_name, template_seq in template:
        template_seq.id = "template"
        if chain not in chains:
            chain = template[chain]

        # If chain was not specified in input, just use sequence from template
        input_chain = SeqRecord(Seq(chains.get(chain_name, template_seq.seq.tostr())), id="model")

        seqlist = [chain_seq, input_chain]
        aln_pir = L_fasta2pir.aln(muscle_aln(seqlist))
        aln_pir.add_pir_info('template', 'structureX', 'template_nucleosome', 'FIRST', i, 'LAST', i)
        aln_pir.add_pir_info('model', 'sequence', 'model_nucleosome')
        seq_aln[chain_name] = aln_pir

    mult_aln = L_fasta2pir.aln_mult_chains(list(seq_aln.values()))
    mult_aln.write('aln.pir')

    # Now let's do MODELLER
    env = environ()  # create a new MODELLER environment to build this model in
    env.io.atom_files_directory = [os.path.sep, "HistoneDB"]

    a = automodel(env,
                  alnfile='aln.pir',
                  knowns='template',
                  sequence='model')
    a.starting_model = 1
    a.ending_model = 1
    a.rand_method = None
    a.max_sc_mc_distance = 10
    a.max_sc_sc_distance = 10
    a.md_level = None  # We will do MD anyway, and the structurea are similar.

    a.make()


outfile = a.outputs[0]['name']

if save_chain:
    return chain,

    # To rewrite using cctbx
    """rc('open %s'%a.outputs[0]['name'])
    rc('open 1kx5.pdb')
    rc('open 1kx5_cen3.pdb')

    rc('matchmaker #1 #0')

    rc('match #2:@N1,N9 #1:@N1,N9')
    rc('changechains A,B I,J #2')

    rc('combine #0,2 newchainids false close false')
    rc('delete :MN')
    rc('delete :HOH')

    rc('resrenumber 127 #3:.A')
    rc('resrenumber 127 #3:.E')

    rc('resrenumber 16 #3:.B')
    rc('resrenumber 16 #3:.F')

    rc('resrenumber 13 #3:.C')
    rc('resrenumber 13 #3:.G')

    rc('resrenumber 36 #3:.D')
    rc('resrenumber 36 #3:.H')

    rc('resrenumber -73 #3:.I')
    rc('resrenumber -73 #3:.J')


    rc('write format pdb #3 temp.pdb')

    rc('close session')

    rc('open nucl_aln.pdb')
    rc('open temp.pdb')
    rc('matchmaker #0 #1')
    rc('write format pdb #1 yeast_nucl_cse4_trunc_aln_v1.pdb')

    os.system('cp yeast_nucl_cse4_trunc_aln_v1.pdb ../2_MD_relax/yeast_cse4/prep/')"""


