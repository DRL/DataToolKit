#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: change_ids_in_gff3.py        -i <FILE> -g <FILE>
                                    [-o <STRING>]
                                                [-h|--help]

    Options:
        -h, --help                          show this
        -i, --id_mapping <FILE>             rename_log.txt
        -g, --gff3 <FILE>                   gff3 file
        -o, --outprefix <STRING>            Output prefix

"""

import sys
import operator
from docopt import docopt
from os.path import isfile, join, exists, realpath, dirname, basename

def read_file(infile):
    if not infile or not exists(infile):
        sys.exit("[ERROR] - File '%s' does not exist." % (infile))
    with open(infile) as fh:
        for line in fh:
            line = line.replace(r'\r','\n')
            yield line.rstrip("\n")

def parse_id_mapping_f(id_mapping_f):
    new_id_by_old_id = {}
    for line in read_file(id_mapping_f):
        temp = line.replace(">", "").split()
        new_id_by_old_id[temp[0]] = temp[1]
    return new_id_by_old_id

def parse_gff3_f(gff3_f, new_id_by_old_id):
    new_lines = []
    for line in read_file(gff3_f):
        temp = []
        new_line = ''
        if line.startswith("#"):
            temp = line.split(" ")
            for idx, field in enumerate(temp):
                if field in new_id_by_old_id:
                    temp[idx] = new_id_by_old_id[field]
            new_line = " ".join(temp)
        else:
            temp = line.split("\t")
            for idx, field in enumerate(temp):
                if field in new_id_by_old_id:
                    temp[idx] = new_id_by_old_id[field]
            new_line = "\t".join(temp)
        new_lines.append(new_line)
    return new_lines

def write_output(out_prefix, new_lines):
    out_f = "renamed.gff3"
    if out_prefix:
        out_f = "%s.renamed.gff3" % (out_prefix)
    with open(out_f, 'w') as out_fh:
        out_fh.write("%s\n" % ("\n".join(new_lines)))

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    id_mapping_f = args['--id_mapping']
    gff3_f = args['--gff3']
    out_prefix = args['--outprefix']

    print "[+] Start ..."
    new_id_by_old_id = parse_id_mapping_f(id_mapping_f)
    new_lines = parse_gff3_f(gff3_f, new_id_by_old_id)
    write_output(out_prefix, new_lines)
