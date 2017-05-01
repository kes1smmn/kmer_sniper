import sys
import os
from arupstl import MongoARUP
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(prog="add additional information to kmer_sniper output",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_file", help="the input file", required=True, type=argparse.FileType("r"))
    parser.add_argument("-f", "--fasta_file", help="the fasta file containing the sequences")
    args = parser.parse_args()
    sniper_output = args.input_file
    fasta_file = SeqIO.parse(args.input_file, "fasta")


    for line in sniper_output:
        print(line)



if __name__ == "__main__":
    main()