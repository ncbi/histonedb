from browse.models import Sequence
from django.conf import settings

import subprocess

from Bio import SeqIO

import os
import sys

def make_blast_db(force=False):
	"""Create new BLAST databse with seuqences in the HistoneDB. This will create a new subset of nr"""
	
	seqs_file = os.path.join(settings.STATIC_ROOT, "browse", "blast", "HistoneDB_sequences.fa")
	"""with open(seqs_file, "w") as seqs:
		for s in Sequence.objects.all():
			SeqIO.write(s.to_biopython(ungap=True), seqs, "fasta")"""

	print " ".join(["makeblastdb", "-in", seqs_file, "-dbtype", "'prot'","-title", "HistoneDB"])
	subprocess.call("makeblastdb", "-in", seqs_file, "-dbtype", "'prot'","-title", "HistoneDB"])

def from_nr(force=True):
	with open("gi_list.txt", "w") as gi_list:
		for gi in Sequence.objects.all().values_list("id", flat=True):
			print >> gi_list, gi
			print gi
	nr_db = os.path.join(settings.STATIC_ROOT, "blast", "nr.db")
	
	if force or not os.path.isfile(nr_db):
		nr = os.path.join(settings.BASE_DIR, "nr")
		if not os.path.isfile(nr):
			raise RuntimeError("Must have nr in projects base diretory, automaticallt placed when the server is built. Please rebuild.")
		subprocess.call([os.path.join(env, "makeblastdb"), "-in", nr, "-dbtype", "prot", "-title", "nr"])
	subprocess.call([os.path.join(env, "blastdbcmd"), "-db", nr, "-dbtype", "prot", "-entry_batch", "gi_list.txt", "-out", "HistoneDB_sequences.fa"])
	subprocess.call([os.path.join(env, "makeblastdb"), "-in", "HistoneDB_sequences.fa", "-dbtype", "prot","-title", "HistoneDB"])

make_blast_db()