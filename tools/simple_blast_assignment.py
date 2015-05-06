import argparse
import csv
import os
import sys
from itertools import izip, groupby

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def assign(seqs, outfile, title):
    """
    """
    assigned = []
    sequences = list(SeqIO.parse(seqs, "fasta"))

    result_handle = NCBIWWW.qblast(
        "blastp", 
        "nr",
        "\n".join([">{}\n{}".format(sequence.id, str(sequence.seq.ungap("-"))) for sequence in sequences]))
    with open("{}_result.txt".format(title), 'w') as results:
        print >> results, result_handle.getvalue()
        
    simple_blast_assignment(seqs, "{}_result.txt".format(title), outfile, title=title)

def simple_blast_assignment(seqs, blast_results, outfile, title=""):
    """Use the best BLAST result to idenitfy GI and organism

    Parameters:
    -----------
    seqS : str
        FASTA filename
    blast_results : File-like object
        Results from BLAST in csv format
    outfile : File-like object
        Where to save the correct sequences
    title : str
        Overall title for multiple sequences. Optional.
    """
    ids = []
    sequences = SeqIO.parse(seqs, "fasta")
    blast_records = NCBIXML.parse(blast_results)
    for n, (sequence, blast_record) in enumerate(zip(sequences, blast_records)):
        print "Sequence", n
        print sequence.format("fasta")

        organism = sequence.id.split("_")[0]
        print organism
        try:
            orig_id = sequence.description.split(" ")[1].strip()
        except IndexError:
            orig_id = None

        best_hsp = None
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print alignment.title
                print organism in alignment.title
                print hsp.query == hsp.sbjct

                if organism not in alignment.title:
                    #Must be from the same organism
                    continue

                best_gi = None
                headers = alignment.title.split(">")
                for header in headers:
                    print header
                    if "{} ".format(organism) in header:
                        best_ids = header.split("|")[1:4:2]
                        if not orig_id:
                            best_gi = best_ids[0]
                        else:
                            #Sanity check if Talber did have the GI
                            if orig_id == best_ids[0]:
                                print "GIIII, yesss"
                                #GIs match
                                best_gi = best_ids[0]
                            elif orig_id in best_ids[1]:
                                print "ID, yes"
                                #Identifier (doesn't includ eperiod adn after)
                                best_gi = best_ids[0]
                            else:
                                #NO match
                                continue
                        break

                if best_gi is None:
                    continue

                if hsp.sbjct_start > 2 or hsp.sbjct_end < len(sequence)-2:
                    continue

                if hsp.query == hsp.sbjct or (orig_id and hsp.query.count("-") == 1 and hsp.query[:hsp.query.index("-")]+hsp.sbjct[hsp.query.index("-")]+hsp.query[hsp.query.index("-")+1:] == hsp.sbjct):
                    #Passed our tests
                    best_hsp = hsp
                    break

            if best_hsp is not None:
                break

        if best_hsp:
            print "****Best Alignment****"
            print "sequence:", alignment.title
            print "length:", alignment.length
            print "e value:", hsp.expect
            print hsp.query
            print hsp.match
            print hsp.sbjct
        
        else:
            print "NO BEST MATCH, pick the best:"
            for i, alignment in enumerate(blast_record.alignments):
                for j, hsp in enumerate(alignment.hsps):
                    print "****Alignment {}, {}****".format(i, j)
                    print "sequence:", alignment.title
                    print "length:", alignment.length
                    print "e value:", hsp.expect
                    print hsp.query
                    print hsp.match
                    print hsp.sbjct
                    print
            best_gi = raw_input("Best GI: ")
            for i, alignment in enumerate(blast_record.alignments):
                for j, hsp in enumerate(alignment.hsps):
                    headers = alignment.title.split(">")
                    for header in headers:
                        if best_gi in header:
                            break
            else:
                continue

        sequence.id = "|".join((organism, best_gi, title))
        sequence.description = sequence.description.replace(" ", "_")

        SeqIO.write(sequence, outfile, "fasta")
	ids.append(best_gi)

    with open("{}.ids.txt".format(outfile.name), "w") as id_file:
        print >> id_file, "\n".join(ids)


def parse_args():
    parser = argparse.ArgumentParser(description="")
     
    #Define Input
    parser.add_argument("sequences",
                        type=argparse.FileType('r'),
                        help="Sequences from FASTA file")
    parser.add_argument("-b", "--blast_results",
                        required=False,
                        type=argparse.FileType('r'),
                        help="Results from BLAST in csv format")
    parser.add_argument("-o", "--outfile",
                        type=argparse.FileType("w+"),
                        required=False,
                        default=None,
                        help="Where to save the correct sequences")
    parser.add_argument("-t", "--title",
                        default="",
                        help="Overall title for multiple sequences.")
    args = parser.parse_args()

    if not args.outfile:
        args.outfile = open("{}_updated.fasta".format(os.path.splitext(os.path.basename(args.sequences.name))[0]), "w")

    return args

if __name__ == "__main__":
    args = parse_args()
    if args.blast_results:
        simple_blast_assignment(args.sequences, args.blast_results, args.outfile, args.title)
    else:
        assign(args.sequences, args.outfile, args.title)



