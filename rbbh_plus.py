#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: rbbh.py                                  -1 <STR> --bed1 <FILE> [-blast1 <FILE>]
                                                -2 <STR> --bed2 <FILE> [-blast2 <FILE>]
                                                [-o <STRING>]
                                                [-h|--help]

    Options:
        -h --help                                   show this
        -1, --species1 <STR>                        Name of species 1
        --bed1 <FILE>                               BED of sp1. Protein_ids have to be in column 4, feature_type in column 8
        --blast1 <FILE>                             BLASTp result using sp1 as query (-max_target_seqs 1 -max_hsps 1 -outfmt '6 'std qlen slen qcovs qcovhsp')
        -2, --species2 <STR>                          Name of species 2
        --bed2 <FILE>                               BED of sp2. Protein_ids have to be in column 4, feature_type in column 8
        --blast2 <FILE>                             BLASTp result using sp2 as query (-max_target_seqs 1 -max_hsps 1 -outfmt '6 'std qlen slen qcovs qcovhsp')
        -o, --outprefix <STRING>                    Output prefix

"""
from __future__ import division
import re
import sys
import operator
from docopt import docopt
from os.path import isfile, join, exists, realpath, dirname, basename

def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    print "[+] Parsing file %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            if not line.startswith("#"):
                yield line.rstrip("\n")

class ProteinObj():
    def __init__(self, protein_id, contig_id, strand)
        self.protein_id = protein_id
        self.contig_id = contig_id
        self.strand = strand
        self.length = 0

class ContigObj():
    def __init__(self, contig_id)
        self.contig_id = contig_id
        self.proteins = []
        self.order_fw = {}
        self.order_rv = {}

    def add_protein(self, protein_id):
        self.proteins.append(protein_id)

    def compute_order(self):
        self.order_fw = {x:i for i,x in enumerate(self.proteins)}
        self.order_rv = {x:i for i,x in enumerate(self.proteins[::-1])}

class GenomeObj():
    def __init__(self, species_id, species_bed_f, species_blast_f):
        self.species_id = species_id
        self.contig_ids = []
        self.contigObjs_by_contig_id = {}
        self.protein_ids = []
        self.proteinObjs_by_protein_id = {}
        self.bed_f = species_bed_f
        self.blast_f = species_blast_f
        self.parse_bed_f(species_bed_f)
        self.parse_blast_f(species_blast_f)

    def parse_bed_f(self, species_bed_f):
        for line in read_file(species_bed_f):
            col = line.split("\t")
            if col[7] == 'mRNA':
                contig_id = col[0]
                protein_id = col[3]
                strand = col[5]
                proteinObj = ProteinObj(protein_id, contig_id, strand)
                self.protein_ids.append(protein_id)
                self.proteinObjs_by_protein_id[protein_id] = proteinObj
                if not contig_id in self.contigObjs_by_contig_id:
                    contigObj = ContigObj(contig_id)
                    self.contig_ids.append(contig_id)
                    self.contigObjs_by_contig_id[contig_id] = contigObj
                self.contigObjs_by_contig_id[contig_id].add_protein(protein_id)
        for contig_id, contigObj in contigObjs_by_contig_id.items():
            contigObj.compute_order()

    def parse_blast_f(self, species_blast_f):
        if species_blast_f:
            for line in read_file(species_bed_f):
                col = [x.strip() for x in line.split("\t")] # cleaning leading and trailing whitespaces, and turning it into a list ...
                if not col[0] in self.proteinObjs_by_protein_id:
                    sys.exit("[X] - qseqid '%s' in BLAST file %s is not part of BED file %s." % (qseqid, self.bed_f, self.blast_f))
                qseqid = col[0]
                sseqid = col[1]
                pident = float(col[2])
                evalue = float(col[10])
                bitscore = int(col[11])
                qlen = int(col[12])
                slen = int(col[13])
                qcov = int(col[14])




if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    species1_id = args['--species1']
    species1_bed_f = args['--bed1']
    species1_blast_f = args['--blast1']
    species2_id = args['--species2']
    species2_bed = args['--bed2']
    species2_blast_f = args['--blast2']
    outprefix = args['--outprefix']

    genomeObj1 = GenomeObj(species1_id, species1_bed_f, species1_blast_f)
    genomeObj2 = GenomeObj(species2_id, species2_bed_f, species2_blast_f)
    rbbh_flag = 0
    if species1_blast_f and species2_blast_f:
        rbbh_flag = 1


