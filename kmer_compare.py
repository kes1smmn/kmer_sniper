from kmer_sniper import query_kmers
from Bio import SeqIO
import argparse
import sys

def reverse_complement(sequence):
    _d = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "Y": "R", "R": "Y", "M": "K", "K": "M",
          "S": "W", "W": "S", "V": "D", "D": "V", "H": "B", "B": "H"}
    return "".join([_d[char] for char in reversed(sequence)])


def kmerize(sequence, kmer_size=31):
    kmer_set = set([])
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence[i : i + kmer_size]
        if "N" in kmer:
            continue  # IGNORE KMERS WITH 'N'
        reverse_complement_kmer = reverse_complement(kmer)
        if kmer < reverse_complement_kmer:
            kmer_set.add(kmer)
        else:
            kmer_set.add(reverse_complement_kmer)
    return kmer_set

def get_sequence(file_path):
    sequence = None
    for i, s in enumerate(SeqIO.parse(file_path, "fasta")):
        if i > 0:
            raise ValueError("The fasta file should only contain one sequence.")
        sequence = str(s.seq).upper()
    if sequence is None:
        raise ValueError("Failed to load sequence [{0}]".format(file_path))
    return sequence




def main():
    parser = argparse.ArgumentParser(prog="compares kmers for two input files",
                                         formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-s", "--kmer_size", help="the input file", default=30, type=int)
    parser.add_argument("-f", "--file_format", help="Format of the input file [default=fa]",
                        choices={"bed", "fa", "fq"}, default="fa")
    parser.add_argument('files', nargs='+',help='fasta file of sequences')

    args = parser.parse_args()
    kmer_size = args.kmer_size
    sequence_files = args.files

    if len(sequence_files) != 2:
        raise ValueError("The number of files should be two.")

    seq_1 = get_sequence(sequence_files[0])
    seq_2 = get_sequence(sequence_files[1])

    kmer_1 = kmerize(seq_1, kmer_size=kmer_size)
    kmer_2 = kmerize(seq_2, kmer_size=kmer_size)

    total_kmers = len(kmer_1.union(kmer_2))
    shared_kmers = len(kmer_1.intersection(kmer_2))
    sys.stdout.write("{0} kmer size\n".format(kmer_size))
    sys.stdout.write("{0:.2f}% {1}[len={2}] : {3}[len={4}]\n".format(shared_kmers/ len(kmer_1) * 100, sequence_files[0],
                                                                     len(seq_1), sequence_files[1], len(seq_2)))



if __name__ == "__main__":
    main()


