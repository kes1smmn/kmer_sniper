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
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def get_sequence_ucsc(chromosome, start, end):
    """
    sends request to ucsc to pull the sequence from hg19

    :param chromosome: can be formatted as 'chr1' or '1'
    :param start: the start position
    :param end: the end position
    :return: string of the sequence
    """
    r = requests.get('http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr{0}:{1},{2}'.format(
        chromosome, start + 1, end), auth=('user', 'pass'))
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
def bin_reads(sequence, kmer_size, qf, bin_cutoff=1):
    found_count = 0
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence[i:i + kmer_size]
        mer = jellyfish.MerDNA(str(kmer))
        mer.canonicalize()
        if qf[mer] > 0:
            found_count += 1
        if found_count >= bin_cutoff:
            return sequence
    return None


def query_kmers(sequence, kmer_size, qf, db_contains_unique_kmer=False, return_non_unique_kmers=False,):
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
    kmer_list = []
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence[i:i + kmer_size]
        mer = jellyfish.MerDNA(str(kmer))
        mer.canonicalize()
        val = qf[mer]
        if db_contains_unique_kmer is False:
            if val == 0:
                val = 1
        values.append(val)
        if return_non_unique_kmers:
            if val > 1:
                kmer_list.append([kmer, val])
            else:
                kmer_list.append([])

    distinct_kmers = len([v for v in values if v > 0])
    non_unique_kmers = len([v for v in values if v > 1])

    if distinct_kmers == 0:  # none of the kmers were found in the database don't report values:
        observational_rank_metric = 0
        non_unique_coverage = 0
    else:

        observational_rank_metric = math.log10(float(numpy.sum(values)) / distinct_kmers)
        non_unique_coverage = float(non_unique_kmers) / distinct_kmers

    if return_non_unique_kmers:
        return values, observational_rank_metric, non_unique_coverage, kmer_list

    return values, observational_rank_metric, non_unique_coverage

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def get_kmer_size(jf_db):
    kmer_size = None
    try:
        mf = jellyfish.ReadMerFile(jf_db)
        for mer, count in mf:
            kmer_size = len(str(mer))
            break
        mf = None
    except:
        raise ValueError("Failed to infer kmer size")
    return kmer_size


def main():
    parser = argparse.ArgumentParser(prog="Returns kmer count of input sequences from jellyfish database",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_file", help="the input file", required=True, type=str)
    parser.add_argument("-f", "--file_format", help="Format of the input file [default=fa]",
                        choices={"bed", "fa", "fq"}, default="fa")
    parser.add_argument("-db", "--jellyfish_database", help="the jf database",
                        default="/Users/331-SimmonkLPTP/Documents/Projects/ARUP/Probe_analysis/kmer_analysis/"
                                "genome_hg19_30_non_unique.jf",
                        type=str)
    parser.add_argument("-v", "--verbose", help="output sequence with only unique kmers too", default=False,
                        action='store_true')

    parser.add_argument("-os", help="output sequence from bed coords; suppress all other output", default=False,
                        action='store_true')

    parser.add_argument("-r", "--reference", help="Reference sequence for use with .BED file input")

    parser.add_argument("--bin", help="will bin the sequence reads if read has a kmer in the jf database",
                        default=False, action='store_true')

    parser.add_argument("-bc", "--bin_cutoff", help="The number of kmer required to bin a read", default=1, type=int)
    parser.add_argument("-n", "--number_of_reads_to_bin", help="The number of reads to bin", default=None, type=int)

    args = parser.parse_args()
    input_file = args.input_file
    jf_db = args.jellyfish_database
    verbose = args.verbose
    file_format = args.file_format
    qf = jellyfish.QueryMerFile(jf_db)
    BIN = args.bin
    bin_cutoff = args.bin_cutoff
    sequence_format = {"fa": "fasta", "fq": "fastq"}
    output_sequence = args.os
    number_of_reads_to_bin = args.number_of_reads_to_bin


    kmer_size = get_kmer_size(jf_db)

    sys.stderr.write("{0}\tkmer size\n".format(kmer_size))
    if output_sequence is False and file_format == 'bed':
        sys.stdout.write("# {0}\tkmer size\n".format(kmer_size))
        sys.stdout.write("# {0}\tjellyfish database\n".format(jf_db))

    # support for bed file [hg19 coordinates]
    if file_format == "bed":
        if BIN:
            raise AttributeError("cannot use 'bin' argument with bed file")
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
    binned_count = 0
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

        if file_format == "fq" and BIN:
            for name, seq, qual in readfq(open(input_file)):

                binned_read_sequence = bin_reads(str(seq), kmer_size, qf, bin_cutoff=bin_cutoff)
                if binned_read_sequence is not None:
                    binned_count += 1
                    sys.stdout.write(">{0}\n{1}\n".format(name, binned_read_sequence))

                if binned_count is not None and binned_count >= number_of_reads_to_bin:
                    break

        else:
            for s in SeqIO.parse(input_file, sequence_format[file_format]):
                if BIN:
                    binned_read_sequence = bin_reads(str(s.seq), kmer_size, qf, bin_cutoff=bin_cutoff)
                    if binned_read_sequence is not None:
                        sys.stdout.write(">{0}\n{1}\n".format(s.name, binned_read_sequence))
                    continue

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
