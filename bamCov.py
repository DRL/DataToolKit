#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: bamCov.py        -b <FILE> --bam_names <FILE> -f <FILE>
                        [--gff <FILE>] [--genelist <FILE>]
                        [-o <STR>] [--contig] [--exon] [--intron]
                        [-h|--help]

    Options:
        -h --help                               show this

        Input files
            -b, --bam_dir <DIR>                 Directory of BAM files
            --bam_names <FILE>                  File with "bamfile = name"
            --genelist <FILE>                   List of gene IDs
            -f, --fasta <FILE>                  Reference FASTA file
            --gff <FILE>                        GFF file
            -o, --outprefix <STR>               Output prefix
            --contig                            Get coverage for contigs [default: False]
            --intron                            Get coverage for introns [default: False]
"""


########################################################################
# Imports
########################################################################

from __future__ import division
from docopt import docopt
import sys
import os
from collections import Counter
import pysam
import pysamstats
import gffutils
import matplotlib as mat
mat.use("agg")
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt


def progress(iteration, steps, max_value, no_limit=False):
    if steps == 0:
        pass
    elif int(iteration) == max_value:
        if no_limit == True:
            sys.stdout.write('\r')
            print "[%%]\t%d%%" % (100),
        else:
            sys.stdout.write('\r')
            print "[%%]\t%d%%" % (100)
    elif int(iteration) % steps == 0:
        sys.stdout.write('\r')
        print "[%%]\t%d%%" % (float(int(iteration) / int(max_value)) * 100),
        sys.stdout.flush()
    else:
        pass

def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    print "[+] Parsing file %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            yield line.rstrip("\n")

class ContigObj():
    def __init__(self, contig_name, contig_seq):
        self.name = contig_name
        #self.seq = contig_seq
        self.length = len(contig_seq)
        #self.counter = Counter(contig_seq)
        #self.base_count = self.counter['A'] + self.counter['G'] + self.counter['C'] + self.counter['T']
        #self.n_count = self.counter['N']
        #if not self.length == (self.base_count + self.n_count):
        #    weird_count = {base: count for base, count in self.counter.items() if base not in {'A', 'G', 'C', 'T', 'N'}}
        #    sys.exit("[X] : Contig %s has non-AGCTN characters in its sequence\n\t%s" % (contig_name, weird_count))

class FeatureObj():
    def __init__(self, chrom, featuretype, name, start, end, strand):
        self.chrom = chrom
        self.type = featuretype
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.length = end - start

class CovObj():
    def __init__(self, featuretype):
        self.featuretype = featuretype
        self.feature_span = 0
        self.covered_span_by_read_cov = Counter()

    def add_coverage(self, covered_span_by_read_cov, span):
        self.covered_span_by_read_cov += covered_span_by_read_cov
        self.feature_span += span

class Main():
    def __init__(self, args):
        # input
        self.bam_dir = args['--bam_dir']
        self.bam_names_f = args['--bam_names']
        self.bam_file_by_name = {}
        self.bam_name_order = []
        self.parse_bam_names()
        self.fasta_f = args['--fasta']
        self.gff_f = args['--gff']
        self.outprefix = args['--outprefix']
        self.gene_list = args['--genelist']
        self.set_of_gene_ids = self.parse_gene_list()
        # setup
        self.contigObjs = []
        self.contigObjs_by_name = {}
        self.parse_fasta()
        # self.featuretypes_of_interest = set(['contig', 'exon'])
        self.featureObjs_by_type = {}
        self.featureObjs_length_by_type = {}
        featuretypes_of_interest = ['exon']
        if args['--contig']:
            featuretypes_of_interest.append('contig')
        if args['--intron']:
            featuretypes_of_interest.append('intron')
        self.featuretypes_of_interest = set(featuretypes_of_interest)
        if 'contig' in self.featuretypes_of_interest:
            self.add_contigs_as_features()
        self.parse_gff()
        self.covObjs_by_bamname_by_featuretype = {}
        self.parse_bams()
        self.plot_cummulative_cov()

    def parse_gene_list(self):
        gene_ids = []
        if self.gene_list:
            for line in read_file(self.gene_list):
                gene_ids.append(line)
        return set(gene_ids)

    def parse_gff(self):
        db_f = self.gff_f + '.db'
        db = ''
        if not os.path.exists(db_f):
            print("[+] - Building GFF database %s" % db_f)
            db = gffutils.create_db(self.gff_f, dbfn=db_f, force=True, keep_order=False, merge_strategy='merge', sort_attribute_values=False)
        else:
            print("[+] Loading GFF database %s" % db_f)
            db = gffutils.FeatureDB(db_f)
        if self.set_of_gene_ids:
            for gene_id in sorted(self.set_of_gene_ids):
                print("[+]\t%s ..." % gene_id)
                gene = db["gene:" + gene_id]
                for featuretype in self.featuretypes_of_interest:
                    for feature in db.children(gene, featuretype=featuretype):
                        featureObj = FeatureObj(feature.seqid, feature.featuretype, feature.id, feature.start, feature.end, feature.strand)
                        self.add_featureObj(featureObj)
        else:
            for featuretype in db.featuretypes():
                if featuretype in self.featuretypes_of_interest:
                    print("[+]\t%ss ..." % featuretype)
                    for feature in gffutils.FeatureDB(db_f).features_of_type(featuretype):
                        featureObj = FeatureObj(feature.seqid, feature.featuretype, feature.id, feature.start, feature.end, feature.strand)
                        self.add_featureObj(featureObj)

    def add_covObj(self, bam_name, covObj):
        if not covObj.featuretype in self.covObjs_by_bamname_by_featuretype:
            self.covObjs_by_bamname_by_featuretype[covObj.featuretype] = {}
        self.covObjs_by_bamname_by_featuretype[covObj.featuretype][bam_name] = covObj

    def plot_cummulative_cov(self):
        for feature_type in self.covObjs_by_bamname_by_featuretype:
            chart_f = "%s.cumulative_cov.pdf" % (feature_type)
            f, ax = plt.subplots(figsize=(10.0, 10.0))
            ax.set_xlabel('Read coverage', fontsize=12)
            ax.set_ylabel('Percentage of bases covered in reference', fontsize=12)
            ax.set_ylim([-0.05, 1.05])
            ax.set_xlim([-5.0, 200.0])
            # plt.margins(0.8)
            # plt.gca().set_ylim(bottom=0.8)
            # plt.gca().set_xlim(left=0.8)
            # ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            # ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            for bam_name in self.bam_name_order:
                covObj = self.covObjs_by_bamname_by_featuretype[feature_type][bam_name]
                x_array_all, y_array_all = self.get_arrays(covObj.covered_span_by_read_cov, covObj.feature_span)
                ax.plot(x_array_all, y_array_all, marker='.', alpha=0.5, label=bam_name)
            ax.legend()
            f.tight_layout()
            #f.suptitle("%s" % feature_type)
            ax.grid(True, linewidth=0.5, which="minor", color="lightgrey")
            ax.grid(True, linewidth=1, which="major", color="lightgrey")
            print "[STATUS] - Plotting %s" % (chart_f)
            f.savefig(chart_f, format='pdf')
            plt.close()

    def get_arrays(self, covDict, total_span):
        fraction = 1.0
        x_values = []
        y_values = []
        for read_cov, covered_span in sorted(covDict.items()):
            y_values.append(fraction)
            x_values.append(read_cov)
            fraction = fraction - (covered_span / total_span)
        return np.array(x_values), np.array(y_values)

    def parse_bams(self):
        for bam_name in self.bam_name_order:
            bam_f = self.bam_file_by_name[bam_name]
            pysamObj = pysam.AlignmentFile(bam_f)
            print("[+] Parsing BAM %s" % bam_f)
            for featuretype in self.featureObjs_by_type:
                count = 0
                print("[+]\t%ss ..." % featuretype)
                covObj = CovObj(featuretype)
                feature_count = len(self.featureObjs_by_type[featuretype])
                for featureObj in self.featureObjs_by_type[featuretype]:
                    count += 1  # debug
                    progress(count, 100, feature_count)
                    read_cov = []
                    for rec in pysamstats.stat_coverage(pysamObj, chrom=featureObj.chrom, start=featureObj.start, end=featureObj.end, truncate=True):
                        read_cov.append(rec['reads_all'])
                    covered_span_by_read_cov = Counter(read_cov)
                    covered_span_by_read_cov[0] = featureObj.length - len(read_cov)
                    covObj.add_coverage(covered_span_by_read_cov, featureObj.length)
                    #if count >= 10:
                    #    print()
                    #    break
                self.add_covObj(bam_name, covObj)

    def parse_fasta(self):
        header, seqs = '', []
        for line in read_file(self.fasta_f):
            if line[0] == '>':
                if header:
                    contigObj = ContigObj(header, ''.join(seqs))
                    self.add_contigObj(contigObj)
                header, seqs = line[1:-1].split()[0], []  # Header is split at first whitespace
            else:
                seqs.append(line[:-1])
        contigObj = ContigObj(header, ''.join(seqs))
        self.add_contigObj(contigObj)

    def parse_bam_names(self):
        for line in read_file(self.bam_names_f):
            col = line.split(" = ")
            bam_f = os.path.join(self.bam_dir, col[0])
            bam_name = col[1]
            self.bam_name_order.append(bam_name)
            self.bam_file_by_name[bam_name] = bam_f

    def add_contigs_as_features(self):
        for contigObj in self.contigObjs:
            featureObj = FeatureObj(contigObj.name, 'contig', contigObj.name, 1, contigObj.length, '+')
            self.add_featureObj(featureObj)

    def add_contigObj(self, contigObj):
        self.contigObjs.append(contigObj)
        self.contigObjs_by_name[contigObj.name] = contigObj

    def add_featureObj(self, featureObj):
        if featureObj.type not in self.featureObjs_by_type:
            self.featureObjs_by_type[featureObj.type] = []
            self.featureObjs_length_by_type[featureObj.type] = 0
        self.featureObjs_by_type[featureObj.type].append(featureObj)
        self.featureObjs_length_by_type[featureObj.type] += featureObj.length

if __name__ == "__main__":
    __version__ = "0.1"
    args = docopt(__doc__)
    Main(args)
