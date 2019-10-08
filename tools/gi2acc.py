import os
import logging

from Bio import SearchIO, SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

seed_directory = "../static/browse/seeds"
seed_acc_dir = "../static/browse/seeds_accession"
Entrez.email = "l.singh@intbio.org"


def get_seeds():
    """
    Goes through static/browse/seeds directories and returns histone type names and fasta file name of variant (without path).
    """
    for i, (root, _, files) in enumerate(os.walk(seed_directory)):
        hist_type = os.path.basename(root)
        if hist_type == "seeds":  # means we are in top dir, we skip,
            # combinded alignmnents for hist types are their, but we do not use them in database constuction,
            # only in visualization on website
            continue
        for seed in files:
            if not seed.endswith(".fasta"): continue
            yield hist_type, seed

def update_seq_rec_from_file(seed_aln_file):
    for s in SeqIO.parse(seed_aln_file, "fasta"):
        gi = s.id.split("|")[1]
        if gi.startswith("NOGI"):
            print("NO GI detected {}".format(s.id))
        else:
            print('-----------------------------------------------------')
            print('GI {}'.format(gi))
            request = Entrez.epost(db="protein", id=str(gi))
            result = Entrez.read(request)
            webEnv = result["WebEnv"]
            queryKey = result["QueryKey"]
            handle = Entrez.efetch(db="protein", rettype='gb', retmode='text', webenv=webEnv, query_key=queryKey)
            for r in SeqIO.parse(handle, 'gb'):
                print r.id
                s.id = s.id.replace(str(gi), str(r.id), 1)
                s.description = s.description.split(' ')[1]
                print('new s.id: {}'.format(s.id))
            print '-------------------------------------------'
        yield s

def main():
    for hist_type, seed in get_seeds():
        variant_name = seed[:-6]
        print(variant_name, "===========")
        seed_aln_file = os.path.join(seed_directory, hist_type, seed)
        sequences = update_seq_rec_from_file(seed_aln_file)
        # for seq in sequences:
        #     print(seq)

        hist_type_dir = "{}/{}".format(seed_acc_dir, hist_type)
        if not os.path.exists(hist_type_dir): os.makedirs(hist_type_dir)
        with open("{}/{}.fasta".format(hist_type_dir, variant_name), "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")


main()