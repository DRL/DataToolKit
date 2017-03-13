#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: bamCov.py        -b <FILE> -f <FILE>
                        [--gff <FILE>] [-g <FILE>]
                        [-o <STR>]
                        [-h|--help]

    Options:
        -h --help                               show this

        Input files
            -b, --bam <FILE>                    BAM file
            -f, --fasta <FILE>                  Reference FASTA file
            --gff <FILE>                        GFF file
            -g, --gene_list <FILE>                           Gene-ID list
            -o, --outprefix <STR>               Output prefix
"""


########################################################################
# Imports
########################################################################

from __future__ import division
from docopt import docopt
import sys
from os.path import isfile, exists
from collections import Counter
import pysam
import pysamstats
import gffutils
import matplotlib as mat
mat.use("agg")
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt


def parse_gene_list(gene_list_f):
    gene_list = []
    if gene_list_f:
        with open(gene_list_f) as gene_list_fh:
            for l in gene_list_fh:
                gene_id = l.rstrip("\n")
                gene_list.append(gene_id)
    return gene_list


def parse_fasta(fasta_f):
    contigCollection = ContigCollection(fasta_f)
    with open(fasta_f) as fasta_fh:
        header, seqs = '', []
        for l in fasta_fh:
            if l[0] == '>':
                if header:
                    contigObj = ContigObj(header, ''.join(seqs))
                    contigCollection.add_contigObj(contigObj)
                header, seqs = l[1:-1].split()[0], []  # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        contigObj = ContigObj(header, ''.join(seqs))
        contigCollection.add_contigObj(contigObj)
    return contigCollection


def parse_gff(gff_f):
    db_fn = inputObj.gff_f + '.db'
    if not exists(db_fn):
        print("[STATUS] - Building GFF database")
        db = gffutils.create_db(inputObj.gff_f, dbfn=inputObj.gff_f + '.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    print("[STATUS] - Loading GFF database")
    return gffutils.FeatureDB(db_fn)


class ContigCollection():
    def __init__(self, fasta_f):
        self.fasta_f = fasta_f
        self.contigObjs = []
        self.contigObjs_by_name = {}

    def add_contigObj(self, contigObj):
        self.contigObjs.append(contigObj)
        self.contigObjs_by_name[contigObj.name] = contigObj


class ContigObj():
    def __init__(self, contig_name, contig_seq):
        self.name = contig_name
        self.seq = contig_seq
        self.length = len(contig_seq)
        self.counter = Counter(contig_seq)
        self.base_count = self.counter['A'] + self.counter['G'] + self.counter['C'] + self.counter['T']
        self.n_count = self.counter['N']
        if not self.length == (self.base_count + self.n_count):
            weird_count = {base: count for base, count in self.counter.items() if base not in {'A', 'G', 'C', 'T', 'N'}}
            sys.exit("[ERROR] : Contig %s has non-AGCTN characters in its sequence\n\t%s" % (contig_name, weird_count))


class FeatureObj():
    def __init__(self, chrom, featuretype, name, start, end, strand):
        self.chrom = chrom
        self.type = featuretype
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.length = end - start


class FeatureCollection():
    def __init__(self):
        self.featureObjs_by_type = {}
        self.featureObjs_length_by_type = {}

    def add_featureObj(self, featureObj):
        if featureObj.type not in self.featureObjs_by_type:
            self.featureObjs_by_type[featureObj.type] = []
            self.featureObjs_length_by_type[featureObj.type] = 0
        self.featureObjs_by_type[featureObj.type].append(featureObj)
        self.featureObjs_length_by_type[featureObj.type] += featureObj.length


class InputObj():
    def __init__(self, args):
        self.bam_f = args['--bam']
        self.check_file(self.bam_f)
        self.ref_f = args['--fasta']
        self.check_file(self.ref_f)
        self.gff_f = args['--gff']
        self.check_file(self.gff_f)
        self.gene_list_f = args['--gene_list']
        self.check_file(self.gene_list_f)
        self.outprefix = args['--outprefix']

    def check_file(self, infile):
        if infile:
            if not isfile(infile):
                sys.exit("[ERROR] : %s does not exist." % (infile))


class CovObj():
    def __init__(self, featuretype):
        self.featuretype = featuretype
        self.feature_span = 0
        self.covered_span_by_read_cov_all = Counter()
        self.covered_span_by_read_cov_pp = Counter()

    def add_coverage(self, covered_span_by_read_cov_all, covered_span_by_read_cov_pp, span):
        self.covered_span_by_read_cov_all += covered_span_by_read_cov_all
        self.covered_span_by_read_cov_pp += covered_span_by_read_cov_pp
        self.feature_span += span


class CovCollection():
    def __init__(self):
        self.covObjs_by_featuretype = {}

    def add_covObj(self, covObj):
        self.covObjs_by_featuretype[covObj.featuretype] = covObj

    def get_arrays(self, covDict, total_span):
        fraction = 1.0
        x_values = []
        y_values = []
        for read_cov, covered_span in sorted(covDict.items()):
            y_values.append(fraction)
            x_values.append(read_cov)
            fraction = fraction - (covered_span / total_span)
        return np.array(x_values), np.array(y_values)

    def plot_cummulative_cov(self, featuretype):
        covObj = self.covObjs_by_featuretype[featuretype]
        chart_f = "%s.cumulative_cov.pdf" % (featuretype)
        f, ax = plt.subplots(figsize=(10.0, 10.0))
        x_array_all, y_array_all = self.get_arrays(covObj.covered_span_by_read_cov_all, covObj.feature_span)
        ax.plot(x_array_all, y_array_all, marker='.', alpha=0.5, label="All reads")
        x_array_pp, y_array_pp = self.get_arrays(covObj.covered_span_by_read_cov_pp, covObj.feature_span)
        ax.plot(x_array_pp, y_array_pp, marker='.', alpha=0.5, label="PP reads")
        ax.set_xlabel('Read coverage', fontsize=12)
        ax.set_ylabel('Percentage of bases covered in reference', fontsize=12)
        ax.set_ylim([-0.05, 1.05])
        ax.set_xlim([-5.0, 200.0])
        #plt.margins(0.8)
        #plt.gca().set_ylim(bottom=0.8)
        #plt.gca().set_xlim(left=0.8)
        #ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax.legend()
        f.tight_layout()
        f.suptitle("%s" % featuretype)
        ax.grid(True, linewidth=0.5, which="minor", color="lightgrey")
        ax.grid(True, linewidth=1, which="major", color="lightgrey")
        print "[STATUS] - Plotting %s" % (chart_f)
        f.savefig(chart_f, format='pdf')
        plt.close()


if __name__ == "__main__":
    __version__ = "0.1"
    args = docopt(__doc__)
    inputObj = InputObj(args)

    # PARSING BAM
    print("[STATUS] - Loading BAM")
    pysamObj = pysam.AlignmentFile(inputObj.bam_f)
    # PARSING REF
    print("[STATUS] - Parsing Reference")
    ContigCollection = parse_fasta(inputObj.ref_f)
    featureCollection = FeatureCollection()
    print("[STATUS] - Adding contigs to data structure ...")
    for contigObj in ContigCollection.contigObjs:
        featureObj = FeatureObj(contigObj.name, 'contig', contigObj.name, 1, contigObj.length, '+')
        featureCollection.add_featureObj(featureObj)

    # PARSING GFF
    gffObj = ''
    if inputObj.gff_f:
        gffObj = parse_gff(inputObj.gff_f)
        print("[STATUS] - Adding features to data structure ...")
        for featuretype in gffObj.featuretypes():
            for feature in gffObj.features_of_type(featuretype):
                featureObj = FeatureObj(feature.seqid, feature.featuretype, feature.id, feature.start, feature.end, feature.strand)
                featureCollection.add_featureObj(featureObj)

    print("[STATUS] - Computing coverage for features ...")
    featuretypes_of_interest = set(['gene', 'contig', 'exon', 'CDS', 'intron'])
    covCollection = CovCollection()
    for featuretype in featureCollection.featureObjs_by_type:
        count = 0
        if featuretype in featuretypes_of_interest:
            print("[STATUS] - \t%s ..." % featuretype)
            covObj = CovObj(featuretype)
            for featureObj in featureCollection.featureObjs_by_type[featuretype]:
                count += 1  # debug
                read_cov_all = []
                read_cov_pp = []
                for rec in pysamstats.stat_coverage(pysamObj, chrom=featureObj.chrom, start=featureObj.start, end=featureObj.end, truncate=True):
                    read_cov_all.append(rec['reads_all'])
                    read_cov_pp.append(rec['reads_pp'])
                covered_span_by_read_cov_all = Counter(read_cov_all)
                covered_span_by_read_cov_all[0] = featureObj.length - len(read_cov_all)
                covered_span_by_read_cov_pp = Counter(read_cov_pp)
                covered_span_by_read_cov_pp[0] = featureObj.length - len(read_cov_pp)
                covObj.add_coverage(covered_span_by_read_cov_all, covered_span_by_read_cov_pp, featureObj.length)
                # if count >= 10:
                #   break
            covCollection.add_covObj(covObj)
            covCollection.plot_cummulative_cov(featuretype)

    # WRITE genes, CDS, EXONS, INTRONS as BEDs at the end
