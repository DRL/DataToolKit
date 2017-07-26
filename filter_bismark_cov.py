#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: filter_bismark_cov.py                    -b <FILE> -v <FILE>
                                                [-o <STRING>]
                                                [-h|--help]

    Options:
        -h --help                                   show this
        -b, --bismark <FILE>                        Bismark coverage
        -v, --vcf <FILE>                            VCF of SNPs to be excluded
        -o, --outprefix <STRING>                    Outprefix
"""
from __future__ import division
import sys
from docopt import docopt
import os


def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    #print "[+] Parsing file %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            if not line.startswith("#"):
                yield line.rstrip("\n")


class BismarkCov():
    def __init__(self, args):
        self.bismark_f = args['--bismark']
        self.vcf_f = args['--vcf']
        self.outprefix = args['--outprefix']
        self.sites_by_contig_id = {}
        self.parse_vcf_f()
        self.parse_bismark_f()

    def parse_vcf_f(self):
        for line in read_file(self.vcf_f):
            col = line.split("\t")
            contig_id = col[0]
            site = int(col[1])
            if contig_id not in self.sites_by_contig_id:
                self.sites_by_contig_id[contig_id] = set()
            self.sites_by_contig_id[contig_id].add(site)

    def parse_bismark_f(self):
        for line in read_file(self.bismark_f):
            col = line.split("\t")
            contig_id = col[0]
            site = int(col[1])
            if contig_id in self.sites_by_contig_id:
                if site not in self.sites_by_contig_id[contig_id]:
                    print line
            else:
                print line


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    bismarkCov = BismarkCov(args)


