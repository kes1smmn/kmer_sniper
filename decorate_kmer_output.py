import sys
from arupstl import MongoARUP
import argparse
from Bio import SeqIO
from collections import OrderedDict, defaultdict
from transcript_bed_store.sub_commands.annotate_bed import annotate_bed
import pkg_resources
import jellyfish
from kmer_sniper import get_kmer_size, query_kmers
import sqlite3
from os.path import join, expanduser
from intervaltree import IntervalTree

conn = MongoARUP()
db = conn.db.transcripts

gene_group = pkg_resources.resource_filename('kmer_sniper', "/resources/gene_group_pseudogenes.txt")
gene_info = pkg_resources.resource_filename('kmer_sniper', "/resources/gene_info_9606.txt")

gene_symbol_to_id = {}
gene_id_to_symbol = {}

base_path = expanduser("~/Documents/Projects/ARUP/psuedogene_comparsion/")
db_path = join(base_path, "kmer_map_non_unique.sqlite3")
kmer_db = sqlite3.connect(db_path)
cur = kmer_db.cursor()

# load resource files
for line in open(gene_info):
    gene_name = line.strip().split("\t")[2]
    gene_id = int(line.strip().split("\t")[1])
    gene_id_to_symbol.update({gene_id: gene_name})
    gene_symbol_to_id.update({gene_name: gene_id})

related_pseudo_genes = defaultdict(set)
for line in open(gene_group):
    line = line.strip().split("\t")
    related_pseudo_genes[int(line[1])].add(int(line[4]))


def reverse_complement(sequence):
    _d = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "Y": "R", "R": "Y", "M": "K", "K": "M",
          "S": "W", "W": "S", "V": "D", "D": "V", "H": "B", "B": "H"}
    return "".join([_d[char] for char in reversed(sequence)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class DecoratedOutput:
    def __init__(self, name, description, rank_score, pct_non_unique, kmer_profile):
        self.name = name
        self.description = description
        self.rank_score = rank_score
        self.pct_non_unique = pct_non_unique
        self.kmer_profile = kmer_profile
        self.kmer_profile_decorated = None
        self.chr, self.start, self.end = self.__split_name(self.name)
        self.kmer_length = len(kmer_profile)
        self.length = self.end - self.start
        self.sequence = None
        self.coding_regions = []
        self.genes = set([])
        self.pseudo_genes = set([])
        self.multi_copy = set([])
        self.kmer_list = None
        self.interval = IntervalTree()
        self.interval.addi(int(self.start), int(self.end))


    def is_multi_copy(self):
        query = {"is_multi_copy": True, "gene_symbol" : {"$in" : list(self.genes)}}
        results = db.find(query)
        for r in results:
            self.multi_copy.add(r["gene_symbol"])


    def get_cds_sequence(self):
        bed_region = "{0},{1},{2}".format(self.chr.upper(), self.start, self.end)
        regions = annotate_bed(None, bed_region=bed_region, annotate_coding_exons=True)

        for reg in regions:
            reg = reg.strip().split("\t")

            if len(reg) <  10:
                continue

            if "cds" in reg[6]:
                self.coding_regions.append(reg[0:3])
                self.genes.update(set(reg[4].split("|")))

        for gene in self.genes:
            gene = gene.replace("MT-", "MT")

            try:
                for related_gene_id in related_pseudo_genes[gene_symbol_to_id[gene]]:
                    try:
                        self.pseudo_genes.add(gene_id_to_symbol[related_gene_id])
                    except KeyError:
                        self.pseudo_genes.add("GeneID:{0}".format(related_gene_id))
            except KeyError:
                pass


    def __split_name(self, name):
        """
        split a name if create from bed coordinate
        :param name:
        :return:
        """
        try:
            chr, coord = name.replace("chr","").split(":")

            chr = chr.upper()

            start, end = coord.split("-")
        except:
            raise ValueError("Could not parse name properly")
        return chr, int(start), int(end)

    def calculated_cds_length(self):
        cds_length = 0
        for cr in self.coding_regions:
            cds_start = int(cr[1])
            cds_end = int(cr[2])
            cds_length += cds_end - cds_start
        self.cds_length = cds_length

    def decorate_kmer_profile(self):
        profile = self.kmer_profile
        for cr in self.coding_regions:
            cds_start = int(cr[1])
            cds_end = int(cr[2])
            adjusted_coordinates_start, adjusted_coordinates_end = cds_start - self.start, cds_end - self.start
            for i in range(adjusted_coordinates_start, adjusted_coordinates_end):
                try:
                    profile[i] = str(profile[i]) + "*"
                except:
                    break
        return profile

    def mapping_output(self):
        s = ""
        # s += ",,,,,,,,,,"
        s += "{0},".format(self.name)
        s += "{0},".format("MAP")
        s += "{0},".format(self.rank_score)
        s += "{0},".format(self.pct_non_unique)
        s += "{0},".format(self.kmer_length)
        s += "{0},".format(self.length)
        s += "{0},".format(self.cds_length)
        s += "{0:.3f},".format(self.cds_length / self.length)
        s += "{0},".format(";".join(self.genes))
        s += "{0},".format(";".join(self.pseudo_genes))
        s += "{0},".format(";".join(self.multi_copy))
        _out = []
        for val in self.kmer_list:
            if val == []:
                _out.append("")
            else:
                _out.append("; ".join(val[2]))
        s += "{0},".format(",".join(_out))
        return s


    def sequence_output(self):
        s = ""
        s += "{0},".format(self.name)
        s += "{0},".format("SEQ")
        s += "{0},".format(self.rank_score)
        s += "{0},".format(self.pct_non_unique)
        s += "{0},".format(self.kmer_length)
        s += "{0},".format(self.length)
        s += "{0},".format(self.cds_length)
        s += "{0:.3f},".format(self.cds_length / self.length)
        s += "{0},".format(";".join(self.genes))
        s += "{0},".format(";".join(self.pseudo_genes))
        s += "{0},".format(";".join(self.multi_copy))
        s += "{0},".format(",".join(self.sequence))
        return s

    def cds_decocrated_output(self):
        s = ""
        #s += ",,,,,,,,,,"
        s += "{0},".format(self.name)
        s += "{0},".format("COUNT")
        s += "{0},".format(self.rank_score)
        s += "{0},".format(self.pct_non_unique)
        s += "{0},".format(self.kmer_length)
        s += "{0},".format(self.length)
        s += "{0},".format(self.cds_length)
        s += "{0:.3f},".format(self.cds_length / self.length)
        s += "{0},".format(";".join(self.genes))
        s += "{0},".format(";".join(self.pseudo_genes))
        s += "{0},".format(";".join(self.multi_copy))
        s += "{0},".format(",".join([str(i) for i in self.decorate_kmer_profile()]))
        return s

def arguments():
    parser = argparse.ArgumentParser(prog="add additional information to kmer_sniper output",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--fasta_file", help="the fasta file containing the sequences")
    parser.add_argument("-db", "--jellyfish_database", help="the jf database",
                        default="/Users/331-SimmonkLPTP/Documents/Projects/ARUP/Probe_analysis/kmer_analysis/"
                                "genome_hg19_30_non_unique.jf",
                        type=str)
    parser.add_argument("-v", "--verbose", help="output sequence with only unique kmers too", default=False,
                        action='store_true')
    parser.add_argument("-db_contains_unique_kmer", default=False, action='store_true')
    return parser.parse_args()



def main():
    args = arguments()
    fasta_file = SeqIO.parse(args.fasta_file, "fasta")
    jf_db = args.jellyfish_database
    verbose = args.verbose
    qf = jellyfish.QueryMerFile(jf_db)
    db_contains_unique_kmer = args.db_contains_unique_kmer

    kmer_size = get_kmer_size(jf_db)
    sniper_objects = OrderedDict()
    for i, s in enumerate(fasta_file):
        kmer_profile, rank_score, pct_non_unique, kmer_list= query_kmers(str(s.seq), kmer_size, qf,
                                                               db_contains_unique_kmer=db_contains_unique_kmer,
                                                               return_non_unique_kmers=True)



        if i % 50 == 0:
            sys.stderr.write("\r{0} sequences analyzed".format(i))
            sys.stderr.flush()
        if verbose is False and pct_non_unique == 0.0:
            continue

        do = DecoratedOutput(s.name, s.description, rank_score, pct_non_unique, kmer_profile)
        do.sequence = s.seq

        # find the positions
        for entry in kmer_list:
            if entry == []:
                continue
            kmer, count = entry
            positions = []
            cur.execute("""SELECT chr, pos FROM kmer_map WHERE kmer == (?)""", (min(kmer, reverse_complement(kmer)),))
            data = cur.fetchall()
            for d in data:
                decorater = ""
                if d[0] == do.chr and do.interval.overlaps(d[1]):
                    decorater = "*"
                if d[0] == "MT":
                    chromosome = "M"
                else:
                    chromosome = d[0]

                positions.append("{0}:{1}-{2}{3}".format(chromosome, d[1], d[1] + 29, decorater))
            entry.append(positions)

        do.kmer_list = kmer_list
        do.get_cds_sequence()
        do.is_multi_copy()
        do.calculated_cds_length()
        sys.stdout.write("{0}\n".format(do.sequence_output()))
        sys.stdout.write("{0}\n".format(do.cds_decocrated_output()))
        sys.stdout.write("{0}\n".format(do.mapping_output()))
    #sys.stderr.write("\r{0} sequences analyzed\n".format(i))
if __name__ == "__main__":
    main()