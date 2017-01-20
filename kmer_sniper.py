from Bio import SeqIO
import jellyfish
import sys
import argparse
import numpy
import requests
import math
import pysam

VALID_CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8",
                     "9", "10", "11", "12", "13", "14", "15",
                     "16", "17", "18", "19", "20", "21", "22",
                     "X", "Y", "M", "MT"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def get_sequence_ucsc(chromosome, start, end):
    """
    sends request to ucsc to pull the sequence from hg19

    :param chromosome: can be formatted as 'chr1' or '1'
    :param start: the start position
    :param end: the end position
    :return: string of the sequence
    """
    r = requests.get('http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr{0}:{1},{2}'.format(
        chromosome, start, end), auth=('user', 'pass'))
    is_seq = False
    sequence = ""
    for line in r.text.split("\n"):
        if is_seq:
            if "</DNA>" in line:
                is_seq = False
            else:
                sequence += line.strip().upper()
        if "<DNA length" in line:
            is_seq = True
    return sequence


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def query_kmers(sequence, kmer_size, qf,):
    """
    takes string (dna sequence) and returns an array containing the kmer counts in a jellyfish database
    for a kmer size

    :param sequence:
    :param kmer_size:
    :param qf:
    :return:
    """
    values = []
    # print(sequence)
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence[i:i + kmer_size]
        mer = jellyfish.MerDNA(str(kmer))
        mer.canonicalize()
        values.append(qf[mer])

    distinct_kmers = len([v for v in values if v > 0])
    non_unique_kmers = len([v for v in values if v > 1])

    if distinct_kmers == 0:  # non of the kmers were found in the database don't report values:
        observational_rank_metric = 0
        non_unique_coverage = 0
    else:
        observational_rank_metric = math.log10(float(numpy.sum(values)) / distinct_kmers)
        non_unique_coverage = float(non_unique_kmers) / distinct_kmers

    return values, observational_rank_metric, non_unique_coverage

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def main():
    parser = argparse.ArgumentParser(prog="Returns kmer count of input sequences from jellyfish database",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_file", help="the input file", required=True, type=str)
    parser.add_argument("-f", "--file_format", help="Format of the input file [default=fa]",
                        choices={"bed", "fa", "fq"}, default="fa")
    parser.add_argument("-db", "--jellyfish_database", help="the jf database",
                        default="/Users/331-SimmonkLPTP/Documents/Projects/ARUP/Probe_analysis/kmer_analysis/"
                                "genome_hg19_30.jf",
                        type=str)
    parser.add_argument("-v", "--verbose", help="output sequence with only unique kmers too", default=False,
                        action='store_true')

    parser.add_argument("-os", help="output sequence from bed coords; suppress all other output", default=False,
                        action='store_true')

    parser.add_argument("-r", "--reference", help="Reference sequence for use with .BED file input")

    args = parser.parse_args()
    input_file = args.input_file
    kmer_size = None
    jf_db = args.jellyfish_database
    verbose = args.verbose
    file_format = args.file_format
    qf = jellyfish.QueryMerFile(jf_db)
    sequence_format = {"fa": "fasta", "fq": "fastq"}
    output_sequence = args.os

    try:
        mf = jellyfish.ReadMerFile(jf_db)
        for mer, count in mf:
            kmer_size = len(str(mer))
            break
        mf = None
    except:
        raise ValueError("Failed to infer kmer size")

    sys.stderr.write("{0}\tkmer size\n".format(kmer_size))
    if output_sequence is False:
        sys.stdout.write("# {0}\tkmer size\n".format(kmer_size))
        sys.stdout.write("# {0}\tjellyfish database\n".format(jf_db))

    # support for bed file [hg19 coordinates]
    if file_format == "bed":

        if args.reference:
            ref_fasta = pysam.FastaFile(args.reference)
        else:
            ref_fasta = None

        for line in open(input_file):
            if line[0] == "#":
                continue

            line = line.strip().split("\t")
            chromosome = line[0].replace("chr", "")

            if chromosome in ["MT", "M", "m", "mt"]:
                chromosome = "M"
            if chromosome in ["y"]:
                chromosome = "Y"
            if chromosome in "x":
                chromosome = "X"

            # if chromosome not in VALID_CHROMOSOMES:
            #     continue
            start, end = int(line[1]), int(line[2])
            if ref_fasta:
                sequence = ref_fasta.fetch(chromosome, start, end)
            else:
                sequence = get_sequence_ucsc(chromosome, start, end)

            if output_sequence:
                sys.stdout.write(">{0}\n{1}\n".format("chr{0}:{1}-{2}".format(chromosome, start, end), sequence))
                continue
            if len(sequence) >= kmer_size:
                values, score_1, score_2 = query_kmers(sequence, kmer_size, qf)
                if len(values) < numpy.sum(values):

                    sys.stdout.write("{0},{1},{2:.3f},{3:.3f},{4}\n".format(
                        "chr{0}:{1}-{2}".format(chromosome, start, end), len(sequence), score_1, score_2,
                        ",".join([str(i) for i in values])))
                else:
                    if verbose:
                        sys.stdout.write("{0},{1},{2:.3f},{3:.3f},{4}\n".format(
                            "chr{0}:{1}-{2}".format(chromosome, start, end), len(sequence), score_1, score_2,
                            ",".join([str(i) for i in values])))
            else:
                sys.stdout.write("{0},{1},{2},{3},{4}\n".format("chr{0}:{1}-{2}".format(chromosome, start, end),
                                                                len(sequence), "na", "na", ",".join(["na"])))

    # support for fasta/fastq files
    if file_format in ["fa", "fq"]:
        # check to see if valid fasta
        for line in open(input_file):
            if file_format == "fa":
                if line[0] != ">":
                    raise IOError("The file does not appear to be valid fasta")
            elif file_format == "fq":
                if line[0] != "@":
                    raise IOError("The file does not appear to be valid fastq")
            break

        for s in SeqIO.parse(input_file, sequence_format[file_format]):
            #  print(str(s))
            values, score_1, score_2 = query_kmers(str(s.seq), kmer_size, qf)

            if len(values) < numpy.sum(values):
                sys.stdout.write("{0},{1:.3f},{2:.3f},{3}\n".format(s.name, score_1, score_2,
                                                                    ",".join([str(i) for i in values])))
            else:
                if verbose:
                    sys.stdout.write("{0},{1:.3f},{2:.3f},{3}\n".format(s.name, score_1, score_2,
                                                                        ",".join([str(i) for i in values])))
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == "__main__":
    main()
