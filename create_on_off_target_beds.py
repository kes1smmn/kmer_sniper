from os.path import expanduser as eu
from os.path import join
from intervaltree import IntervalTree
from collections import defaultdict, Counter
import os
import argparse
import sys



def parse_pos(pos):
    pos = pos.replace("*", "")
    chr, interval = pos.split(":")
    start, end = interval.split("-")
    return chr.strip(), int(start), int(end)

def output_bedfile_format(chromosome_interval_tree):
    order = [str(i) for i in range(1, 22)] + ["X", "Y", "MT", "M"]

    # for probe, it in chromosome_interval_tree.items():
    #     print(it)
    #     it.merge_overlaps(data_reducer=merge_data)
    #     for i in sorted(it):
    #
    #
    #         yield "{0}\t{1}\t{2}\t{3}\t{4}\n".format(probe.split(":")[0], i[0], i[1],
    #                                                  i[1] + 1 - i[0], ",".join(set(i[2])))
    for chromosome in order:
        if chromosome in chromosome_interval_tree:
            for name, it in chromosome_interval_tree[chromosome].items():
                it.merge_overlaps(data_reducer=merge_data)
                for i in sorted(it):
                    ordered_regions = dict(Counter(i[2]))
                    data = ["{0}={1}".format(k, ordered_regions[k]) for k in sorted(ordered_regions,
                                                                                    key=ordered_regions.get, reverse=True)]
                    yield "{0}\t{1}\t{2}\t{3}\t{4}\n".format(chromosome, i[0], i[1], i[1] + 1 - i[0], ",".join(data))


def arguments():
    parser = argparse.ArgumentParser(prog="add additional information to kmer_sniper output",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_file", help="the decorated output", type=argparse.FileType("r"))

    parser.add_argument("-o", "--output_prefix")
    return parser.parse_args()


base_path = eu("~/Documents/Projects/ARUP/psuedogene_comparsion/covered_bedfile_analysis/Capture-Bedfiles")
input_file = join(base_path, "0747991_Covered_v1_kmer.csv")

def main():
    args = arguments()
    create_off_on_target_beds(args.input_file, args.output_prefix)

def merge_data(data, value):
    data = data + value
    return data

def create_off_on_target_beds(input_file, output_prefix):
    on_target_locations = defaultdict(dict)
    off_target_locations = defaultdict(dict)


    for i, line in enumerate(input_file):
        line = line.strip().split(",")
        if line[1] != "MAP":
            continue


        # for it in on_target_locations.values():
        #     it.merge_overlaps(data_reducer=merge_data)
        # for it in off_target_locations.values():
        #     it.merge_overlaps(data_reducer=merge_data)
        # sys.stderr.write("{0}\n".format(i))


        name = line[0]
        for kmer_position in line[11:]:
            mapping_positions = kmer_position.split(";")
            if mapping_positions == ['']:
                continue
            for pos in mapping_positions[0:10]:
                chromosome, start, stop = parse_pos(pos)
                if "*" in pos:
                    if name not in on_target_locations[chromosome]:
                        on_target_locations[chromosome][name] = IntervalTree()
                    on_target_locations[chromosome][name].addi(start, stop, [name])

                else:
                    if name not in off_target_locations[chromosome]:
                        off_target_locations[chromosome][name] = IntervalTree()
                    off_target_locations[chromosome][name].addi(start, stop, [name])


    # for reg in output_bedfile_format(on_target_locations):
    #    sys.stderr.write(reg)

    with open(join(os.getcwd(), output_prefix + "_on_target_non-unique_kmer_regions.bed"), "w") as wf:
        for reg in output_bedfile_format(on_target_locations):
            wf.write(reg)

    with open(join(os.getcwd(), output_prefix + "_off_target_non-unique_kmer_regions.bed"), "w") as wf:
        for reg in output_bedfile_format(off_target_locations):
            wf.write(reg)
    return

if __name__ == "__main__":
    main()


