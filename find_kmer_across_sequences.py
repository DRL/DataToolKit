#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage:
    find_kmer_across_sequences.py                    -f <FILE> [-m <INT>] [-h|--help]

    Options:
        -h, --help                              show this
        -f, --fasta <FILE>                      FASTA file
        -m, --max_kmer <INT>                    max kmer (default: 1000)
"""

from __future__ import division
from docopt import docopt
import os
import sys


def read_fasta(fasta_f):
    header, seq = '', ''
    with open(fasta_f) as fasta_fh:
        for line in fasta_fh:
            if line.startswith(">"):
                if (seq):
                    yield header, seq
                    seq = ''
                header = line[1:-1]
            else:
                seq += line[:-1]
        yield header, seq


def write_file(out_f, outprefix, header, strings):
    if outprefix:
        if outprefix.endswith("/"):
            if not os.path.exists(outprefix):
                os.mkdir(outprefix)
            out_f = "%s" % os.path.join(outprefix, out_f)
        else:
            out_f = "%s.%s" % (outprefix, out_f)
    print "[+] \t Writing file %s ..." % (out_f)
    with open(out_f, 'w') as out_fh:
        out_fh.write("%s\n" % (header))
        out_fh.write("%s\n" % "\n".join(strings))

def revcom(sequence):
    complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    revcom_seq = "".join([complement.get(nt.upper(), '') for nt in sequence[::-1]])
    return revcom_seq

class SeqObj():
    def __init__(self, header, seq):
        self.header = header
        self.seq = seq
        self.length = len(seq)


class DataCollection():
    def __init__(self, args):
        self.fasta_f = args['--fasta']
        self.max_kmer = int(args['--max_kmer'])
        self.seqObjs = self.parse_fasta_f()
        self.kmers = self.generate_kmers()

    def parse_fasta_f(self):
        seqObjs = []
        for header, seq in read_fasta(self.fasta_f):
            seqObj = SeqObj(header, seq)
            seqObjs.append(seqObj)
        return seqObjs

    def generate_kmers(self):
        for k in xrange(self.max_kmer, 0, -1):
            print k
            headers_by_kmer = {}
            for seqObj in self.seqObjs:
                kmers = []
                for i in xrange(seqObj.length - k + 1):
                    kmer = seqObj.seq[i:i + k]
                    kmers.append(kmer)
                    kmer_rev = revcom(kmer)
                    kmers.append(kmer_rev)
                for kmer in kmers:
                    if kmer not in headers_by_kmer:
                        headers_by_kmer[kmer] = []
                    headers_by_kmer[kmer].append(seqObj.header)
                    if len(headers_by_kmer[kmer]) == len(self.seqObjs):
                        print kmer
                        print headers_by_kmer[kmer]
                        sys.exit("Done")


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    print args
    print "[+] Start ..."
    dataCollection = DataCollection(args)

