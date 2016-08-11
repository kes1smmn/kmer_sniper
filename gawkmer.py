from Bio import SeqIO
import jellyfish
import sys
import argparse
import numpy
import requests
import math


def get_sequence(chromosome, start, end):
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
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence[i:i + kmer_size]
        mer = jellyfish.MerDNA(str(kmer))
        mer.canonicalize()
        values.append(qf[mer])

    distinct_kmers = len([v for v in values if v > 0])
    non_unique_kmers = len([v for v in values if v > 1])
    observational_rank_metric = math.log10(float(numpy.sum(values)) / distinct_kmers)
    non_unique_coverage = float(non_unique_kmers) / distinct_kmers
    return values, observational_rank_metric, non_unique_coverage


def main():
    parser = argparse.ArgumentParser(prog="Returns kmer count of input sequences from jellyfish database",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_file", help="the input file", required=True, type=str)
    parser.add_argument("-f", "--file_format", help="Format of the input file [default=fa]",
                        choices={"bed", "fa", "fq"}, default="fa")
    parser.add_argument("-db", "--jellyfish_database", help="the jf database", required=True, type=str)
    parser.add_argument("-v", "--verbose", help="output sequence with only unique kmers too", default=False,
                        action='store_true')
    args = parser.parse_args()
    input_file = args.input_file
    kmer_size = None
    jf_db = args.jellyfish_database
    verbose = args.verbose
    file_format = args.file_format
    qf = jellyfish.QueryMerFile(jf_db)
    sequence_format = {"fa": "fasta", "fq": "fastq"}

    try:
        mf = jellyfish.ReadMerFile(jf_db)
        for mer, count in mf:
            kmer_size = len(str(mer))
            break
        mf = None
    except:
        raise ValueError("Failed to infer kmer size")

    sys.stderr.write("{0}\tkmer size".format(kmer_size))
    sys.stdout.write("# {0}\tkmer size".format(kmer_size))
    sys.stdout.write("# {0}\tjellyfish database".format(kmer_size))

    # support for bed file [hg19 coordinates]
    if file_format == "bed":
        for line in open(input_file):
            if line[0] != "#":
                line = line.strip().split("\t")
                chromosome, start, end = line[0].replace("chr", ""), int(line[1]), int(line[2])
                sequence = get_sequence(chromosome, start, end)
                if len(sequence) >= kmer_size:
                    values, score_1, score_2 = query_kmers(sequence, kmer_size, qf)

                    if len(values) < numpy.sum(values):
                        sys.stdout.write("{0},{1:.3f},{2:.3f},{3}\n".format(
                            "chr{0}:{1}-{2}".format(chromosome, start, end), score_1, score_2,
                            ",".join([str(i) for i in values])))
                    elif verbose:
                        sys.stdout.write("{0},{1:.3f},{2:.3f},{3}\n".format(
                            "chr{0}:{1}-{2}".format(chromosome, start, end), score_1, score_2,
                            ",".join([str(i) for i in values])))
                else:
                    sys.stdout.write("{0},{1},{2},{3}\n".format("chr{0}:{1}-{2}".format(chromosome, start, end), "na",
                                                                "na", ",".join(["na"])))

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
            values, score_1, score_2 = query_kmers(str(s.seq), kmer_size, qf)
            if len(values) < numpy.sum(values):
                sys.stdout.write("{0},{1:.3f},{2:.3f},{3}\n".format(s.name, score_1, score_2,
                                                                    ",".join([str(i) for i in values])))
            elif verbose:
                sys.stdout.write("{0},{1:.3f},{2:.3f},{3}\n".format(s.name, score_1, score_2,
                                                                    ",".join([str(i) for i in values])))

if __name__ == "__main__":
    main()
