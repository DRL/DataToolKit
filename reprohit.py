#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Usage:
        genome_analyses.py analyse              -d <DIR> [rbbh] [ss] [-o <STRING>]
        genome_analyses.py plot                 -d <DIR> [rbbh] [ss] [protein] [contig] [genome] [-o <STRING>]

Options:
        -h --help                                show this
        -d, --input_dir <DIR>                    Input directory
        -o, --outprefix <STRING>                 Output prefix

"""
from __future__ import division
from docopt import docopt
from collections import Counter
import sys
import os
import matplotlib as mat
mat.use("agg")
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    print "[+] Parsing file %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            yield line.rstrip("\n")


def write_file(out_f, outprefix, header, string):
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
        out_fh.write("%s\n" % "\n".join(string))


def getRegion(seq, start, stop, strand):
    region = seq[start:stop]
    if strand == '-':
        complement = {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A',
            'N': 'N'}
        region = "".join([complement.get(nt.upper(), '') for nt in region[::-1]])
    elif strand == '+':
        pass
    else:
        sys.exit("[ERROR] - strand should be +/-, not : %s" % (strand))
    return region


class FeatureObj():
    def __init__(self, protein_id, contig_id, feature_type, seq, start, stop, strand):
        self.protein_id = protein_id
        self.contig_id = contig_id
        self.feature_type = feature_type
        self.start = start
        self.stop = stop
        self.strand = strand
        self.seq = seq
        self.length = len(seq)
        self.N_count = seq.count("N")


class ProteinObj():
    def __init__(self, protein_id, contig_id, strand):
        self.protein_id = protein_id
        self.contig_id = contig_id
        self.strand = strand
        self.length = 0
        self.cdsObjs = []
        self.CDS_count = 0
        self.intronObjs = []
        self.intron_count = 0

    def add_featureObj(self, featureObj):
        if featureObj.feature_type == "CDS":
            self.cdsObjs.append(featureObj)
            self.CDS_count += 1
        elif featureObj.feature_type == "intron":
            self.intronObjs.append(featureObj)
            self.intron_count += 1
        else:
            sys.exit("[ERROR] - unknown feature_type %s" % (featureObj.feature_type))

    def get_gc(self, feature_type):
        featureObjs = self.get_featureObjs(feature_type)
        if featureObjs:
            gcs = [featureObj.seq.count("G") + featureObj.seq.count("C") for featureObj in featureObjs]
            lengths = [(featureObj.length - featureObj.N_count) for featureObj in featureObjs]
            return sum(gcs) / sum(lengths)
        else:
            return "-"

    def get_featureObjs(self, feature_type):
        featureObjs = []
        if feature_type == "CDS":
            featureObjs = self.cdsObjs
        elif feature_type == "intron":
            featureObjs = self.intronObjs
        else:
            sys.exit("[ERROR] - unknown feature_type %s" % (feature_type))
        if self.strand == '-':
            featureObjs = featureObjs[::-1]
        return featureObjs

    def get_splicesites(self):
        intronObjs = self.get_featureObjs('intron')
        acceptors = []
        donors = []
        for intronObj in intronObjs:
            donors.append(intronObj.seq[0:NUCCOUNT])
            acceptors.append(intronObj.seq[-NUCCOUNT:])
        return ["%s/%s" % (d, a) for d, a in zip(donors, acceptors)]


class ContigObj():
    def __init__(self, contig_id, seq):
        self.contig_id = contig_id
        self.seq = seq
        self.length = len(seq)
        self.N_count = seq.count("N")
        self.protein_ids = []
        self.protein_count = 0
        self.CDS_count = 0
        self.intron_count = 0
        self.order_fw = {}
        self.order_rv = {}

    def add_protein_id(self, protein_id):
        self.protein_ids.append(protein_id)
        self.protein_count += 1

    def get_gc(self):
        return (self.seq.count("G") + self.seq.count("C")) / (self.length - self.N_count)

    def yield_protein_ids(self):
        for protein_id in self.protein_ids:
            yield protein_id

    def compute_order(self):
        self.order_fw = {x: i for i, x in enumerate(self.protein_ids)}
        self.order_rv = {x: i for i, x in enumerate(self.protein_ids[::-1])}


class GenomeObj():
    def __init__(self, genome_id, species_fasta_f, species_bed_f, species_blast_f):
        self.genome_id = genome_id
        self.length = 0
        self.N_count = 0
        self.contig_count = 0
        self.CDS_count = 0
        self.intron_count = 0
        self.protein_count = 0
        self.contig_ids = []
        self.contigObjs_by_contig_id = {}
        self.protein_ids = []
        self.proteinObjs_by_protein_id = {}
        self.bed_f = species_bed_f
        self.hitObjs = []
        self.hit_count = 0
        self.blast_f = species_blast_f
        self.parse_fasta_f(species_fasta_f)
        self.parse_bed_f(species_bed_f)
        self.parse_blast_f(species_blast_f)

    def add_contigObj(self, contigObj):
        self.contig_count += 1
        self.N_count += contigObj.N_count
        self.length += contigObj.length
        self.contig_ids.append(contigObj.contig_id)
        self.contigObjs_by_contig_id[contigObj.contig_id] = contigObj

    def parse_fasta_f(self, species_fasta_f):
        header, seqs = '', []
        for line in read_file(species_fasta_f):
            if line[0] == '>':
                if header:
                    contigObj = ContigObj(header, ''.join(seqs))
                    self.add_contigObj(contigObj)
                header, seqs = line[1:-1].split()[0], []  # Header is split at first whitespace
            else:
                seqs.append(line)
        contigObj = ContigObj(header, ''.join(seqs))
        self.add_contigObj(contigObj)
        print "[+] \t %s contigs parsed ..." % (self.contig_count)

    def add_proteinObj(self, proteinObj):
        self.protein_count += 1
        self.protein_ids.append(proteinObj.protein_id)
        self.proteinObjs_by_protein_id[proteinObj.protein_id] = proteinObj
        contigObj = self.contigObjs_by_contig_id[proteinObj.contig_id]
        contigObj.add_protein_id(proteinObj.protein_id)

    def add_featureObj(self, featureObj):
        proteinObj = self.proteinObjs_by_protein_id[featureObj.protein_id]
        contigObj = self.contigObjs_by_contig_id[featureObj.contig_id]
        if featureObj.feature_type == "CDS":
            self.CDS_count += 1
            contigObj.CDS_count += 1
            proteinObj.add_featureObj(featureObj)
        elif featureObj.feature_type == "intron":
            self.intron_count += 1
            contigObj.intron_count += 1
            proteinObj.add_featureObj(featureObj)
        else:
            sys.exit("[ERROR] - unknown feature_type %s" % (featureObj.feature_type))

    def parse_bed_f(self, species_bed_f):
        features_to_parse = set(['mRNA', 'CDS', 'intron'])
        for line in read_file(species_bed_f):
            col = line.split("\t")
            contig_id = col[0]
            protein_id = col[3]
            start = int(col[1])
            stop = int(col[2])
            strand = col[5]
            if col[7] in features_to_parse:
                contigObj = self.contigObjs_by_contig_id[contig_id]
                feature_type = col[7]
                if feature_type == 'mRNA':
                    if protein_id not in self.proteinObjs_by_protein_id:
                        self.add_proteinObj(ProteinObj(protein_id, contig_id, strand))
                elif feature_type == 'CDS':
                    if protein_id not in self.proteinObjs_by_protein_id:
                        self.add_proteinObj(ProteinObj(protein_id, contig_id, strand))
                    cds_seq = getRegion(self.contigObjs_by_contig_id[contig_id].seq, start, stop, strand)
                    cdsObj = FeatureObj(protein_id, contig_id, 'CDS', cds_seq, start, stop, strand)
                    self.add_featureObj(cdsObj)
                elif feature_type == 'intron':
                    if protein_id not in self.proteinObjs_by_protein_id:
                        self.add_proteinObj(ProteinObj(protein_id, contig_id, strand))
                    intron_seq = getRegion(self.contigObjs_by_contig_id[contig_id].seq, start, stop, strand)
                    intronObj = FeatureObj(protein_id, contig_id, 'intron', intron_seq, start, stop, strand)
                    self.add_featureObj(intronObj)
                else:
                    pass
        for contig_id, contigObj in self.contigObjs_by_contig_id.items():
            contigObj.compute_order()
        print "[+] \t %s mRNAs parsed ..." % (self.protein_count)
        print "[+] \t %s CDSs parsed ..." % (self.CDS_count)
        print "[+] \t %s introns parsed ..." % (self.intron_count)

    def parse_blast_f(self, species_blast_f):
        if species_blast_f:
            for line in read_file(species_blast_f):
                col = [x.strip() for x in line.split("\t")]  # cleaning leading and trailing whitespaces, and turning it into a list ...
                if not col[0] in self.proteinObjs_by_protein_id:
                    sys.exit("[X] qseqid '%s' in BLAST file %s is not part of BED file %s." % (col[0], self.blast_f, self.bed_f))
                try:
                    hitObj = HitObj(col)
                except TypeError:
                    sys.exit("[X] BLAST outfmt should bed '6 std qlen slen qcovs qcovhsp'")
                self.hitObjs.append(hitObj)
                self.hit_count += 1
            print "[+] \t %s hits parsed ..." % (self.hit_count)


class HitObj():
    def __init__(self, col):
        self.qseqid = col[0]
        self.sseqid = col[1]
        self.pair = frozenset([col[0], col[1]])
        self.pident = float(col[2])
        self.evalue = float(col[10])
        self.bitscore = float(col[11])
        self.qlen = int(col[12])
        self.slen = int(col[13])
        self.qcov = int(col[14])


class ResultsObj():
    def __init__(self, analysis, genome_ids, outprefix):
        self.analysis = analysis
        self.genome_ids = genome_ids
        self.outprefix = outprefix
        self.headers_by_level = {}
        self.genome_results_by_genome_id = {genome_id: '' for genome_id in self.genome_ids}
        self.protein_results_by_genome_id = {genome_id: [] for genome_id in self.genome_ids}
        self.contig_results_by_genome_id = {genome_id: [] for genome_id in self.genome_ids}

    def write_result(self):
        genome_out_f = "%s.genomes.txt" % (self.analysis)
        genome_string = []
        for genome_id in self.genome_ids:
            genome_string.append(self.genome_results_by_genome_id[genome_id])
            contig_out_f = "%s.%s.contig.txt" % (genome_id, self.analysis)
            write_file(contig_out_f, self.outprefix, self.headers_by_level['contig'], self.contig_results_by_genome_id[genome_id])
            protein_out_f = "%s.%s.protein.txt" % (genome_id, self.analysis)
            write_file(protein_out_f, self.outprefix, self.headers_by_level['protein'], self.protein_results_by_genome_id[genome_id])
        write_file(genome_out_f, self.outprefix, self.headers_by_level['genome'], genome_string)


class AnalysisCollection():
    def __init__(self, genomeObjs):
        self.genomeObjs_by_genome_id = {genomeObj.genome_id: genomeObj for genomeObj in genomeObjs}
        self.genome_ids = sorted([genomeObj.genome_id for genomeObj in genomeObjs])
        self.genome_count = len(self.genome_ids)

    def analyse_genomes(self, analysis):
        if analysis == 'splicesites':
            self.analyse_splice_sites()
        elif analysis == 'rbbh':
            if not self.genome_count == 2:
                sys.exit("[X] - RBBH analyses only implemented for two genomes. Analysis contains %s genomes." % (self.genome_count))
            self.analyse_rbbhs()
        else:
            sys.exit("[X] - Analysis %s is not implemented." % (analysis))

    def analyse_rbbhs(self):
        print "[+] Calculate RBBHs (reciprocal best BLAST hits) ..."
        hitObjs = []
        for genome_id in self.genome_ids:
            for hitObj in self.genomeObjs_by_genome_id[genome_id].hitObjs:
                hitObjs.append(hitObj)
        seen_pairs = set()
        seen_proteins = set()
        rbbh_protein_id_pairs = []
        length_difference_by_pair = {}
        bitscores_by_pair = {}
        for hitObj in sorted(hitObjs, key=lambda x: x.bitscore, reverse=True):
            if hitObj.evalue <= EVALUE and hitObj.qcov >= QCOV:
                if hitObj.pair in seen_pairs:
                    rbbh_protein_id_pairs.append(hitObj.pair)
                    bitscores_by_pair[hitObj.pair].append(hitObj.bitscore)
                    length_difference_by_pair[hitObj.pair] = abs(hitObj.qlen - hitObj.slen)
                else:
                    if hitObj.pair.intersection(seen_proteins):
                        seen_proteins.add(hitObj.qseqid)
                        seen_proteins.add(hitObj.sseqid)
                    else:
                        bitscores_by_pair[hitObj.pair] = [hitObj.bitscore]
                        seen_pairs.add(hitObj.pair)
        header = "#%s\tlength_diff\tmean_bitscore" % "\t".join(self.genome_ids)
        body = []
        for rbbh_protein_id_pair in rbbh_protein_id_pairs:
            protein_a, protein_b = rbbh_protein_id_pair
            if protein_a in self.genomeObjs_by_genome_id[self.genome_ids[0]].proteinObjs_by_protein_id:
                body.append("%s\t%s\t%s\t%s" % (protein_a, protein_b, length_difference_by_pair[rbbh_protein_id_pair], sum(bitscores_by_pair[rbbh_protein_id_pair]) / 2))
            else:
                body.append("%s\t%s\t%s\t%s" % (protein_b, protein_a, length_difference_by_pair[rbbh_protein_id_pair], sum(bitscores_by_pair[rbbh_protein_id_pair]) / 2))
        rbbh_out_f = "%s.%s.%s.%s.txt" % (".vs.".join(self.genome_ids), EVALUE, QCOV, 'rbbh')
        write_file(rbbh_out_f, outprefix, header, body)

    def analyse_splice_sites(self):
        print "[+] Analysing splice sites ..."
        resultsObj = ResultsObj('splicesites', self.genome_ids, outprefix)
        resultsObj.headers_by_level['genome'] = "#%s" % ("\t".join([
            "genome_id",
            "genome_length",
            "N_count",
            "contig_count",
            "protein_count",
            "CDS_count",
            "intron_count",
            "non_canonical_ss_frac",
            "ss_counts"]))
        resultsObj.headers_by_level['contig'] = "#%s" % ("\t".join([
            "genome_id",
            "contig_id",
            "contig_length",
            "N_count",
            "contig_gc"
            "protein_count",
            "CDS_count",
            "intron_count",
            "non_canonical_ss_frac",
            "ss_counts"]))
        resultsObj.headers_by_level['protein'] = "#%s" % ("\t".join([
            "genome_id",
            "protein_id",
            "contig_id",
            "CDS_count",
            "CDS_gc",
            "intron_count",
            "intron_gc",
            "non_canonical_ss_frac",
            "ss_counts",
            "ss_string"]))
        for genome_id in self.genome_ids:
            genomeObj = self.genomeObjs_by_genome_id[genome_id]
            genome_splicesites = []
            for contig_id in genomeObj.contig_ids:
                contigObj = genomeObj.contigObjs_by_contig_id[contig_id]
                contig_splicesites = []
                for protein_id in contigObj.protein_ids:
                    proteinObj = genomeObj.proteinObjs_by_protein_id[protein_id]
                    protein_line = []
                    protein_splicesites = proteinObj.get_splicesites()
                    protein_line.append(genomeObj.genome_id)
                    protein_line.append(proteinObj.protein_id)
                    protein_line.append(proteinObj.contig_id)
                    protein_line.append(proteinObj.CDS_count)
                    protein_line.append(proteinObj.get_gc('CDS'))
                    protein_line.append(proteinObj.intron_count)
                    protein_line.append(proteinObj.get_gc('intron'))
                    if proteinObj.intron_count > 0:
                        protein_splicesites_counter = Counter(protein_splicesites)
                        protein_line.append("%s" % ("{:.2}".format(protein_splicesites_counter["GT/AG"] / proteinObj.intron_count)))
                        protein_line.append("%s" % (";".join(["%s=%s" % (k, v) for k, v in protein_splicesites_counter.most_common()])))
                        protein_line.append("%s" % (",".join(protein_splicesites)))
                        contig_splicesites[0:0] = protein_splicesites  # fastest way of concatenating lists
                        genome_splicesites[0:0] = protein_splicesites  # fastest way of concatenating lists
                    else:
                        protein_line.append("-")
                        protein_line.append("-")
                        protein_line.append("-")
                    resultsObj.protein_results_by_genome_id[genomeObj.genome_id].append("%s" % "\t".join([str(x) for x in protein_line]))
                contig_line = []
                contig_line.append(genomeObj.genome_id)
                contig_line.append(contigObj.contig_id)
                contig_line.append(contigObj.length)
                contig_line.append(contigObj.N_count)
                contig_line.append(contigObj.get_gc())
                contig_line.append(contigObj.protein_count)
                contig_line.append(contigObj.CDS_count)
                contig_line.append(contigObj.intron_count)
                if contigObj.intron_count > 0:
                    contig_splicesites_counter = Counter(contig_splicesites)
                    contig_line.append("%s" % ("{:.2}".format(contig_splicesites_counter["GT/AG"] / contigObj.intron_count)))
                    contig_line.append("%s" % (";".join(["%s=%s" % (k, v) for k, v in contig_splicesites_counter.most_common()])))
                else:
                    contig_line.append("-")
                    contig_line.append("-")
                resultsObj.contig_results_by_genome_id[genomeObj.genome_id].append("%s" % "\t".join([str(x) for x in contig_line]))
            genome_line = []
            genome_line.append(genomeObj.genome_id)
            genome_line.append(genomeObj.length)
            genome_line.append(genomeObj.N_count)
            genome_line.append(genomeObj.contig_count)
            genome_line.append(genomeObj.protein_count)
            genome_line.append(genomeObj.CDS_count)
            genome_line.append(genomeObj.intron_count)
            genome_splicesites_counter = Counter(genome_splicesites)
            genome_line.append("%s" % ("{:.2}".format(genome_splicesites_counter["GT/AG"] / genomeObj.intron_count)))
            genome_line.append("%s" % (";".join(["%s=%s" % (k, v) for k, v in genome_splicesites_counter.most_common()])))
            resultsObj.genome_results_by_genome_id[genomeObj.genome_id] = "%s" % "\t".join([str(x) for x in genome_line])
        resultsObj.write_result()


def generate_genomeObjs(input_dir):
    dir_idx = 0
    genomeObjs = []
    genome_ids = []
    for root, dirs, files in os.walk(input_dir):
        if dir_idx == 0:
            genome_ids = dirs
        else:
            genome_id = genome_ids[dir_idx - 1]
            paths = {'bed_f': '', 'fasta_f': '', 'blast_f': ''}
            for infile in files:
                if infile.endswith(".bed"):
                    paths['bed_f'] = os.path.join(root, infile)
                elif infile.endswith(".sl.fa"):
                    paths['fasta_f'] = os.path.join(root, infile)
                elif infile.endswith(".out"):
                    paths['blast_f'] = os.path.join(root, infile)
                else:
                    pass
            genomeObj = GenomeObj(genome_id, paths['fasta_f'], paths['bed_f'], paths['blast_f'])
            genomeObjs.append(genomeObj)
        dir_idx += 1
    return genomeObjs


class PlottingCollection():
    def __init__(self, args):
        self.input_dir = args['--input_dir']
        self.rbbh_flag = args['rbbh']
        self.ss_flag = args['ss']
        self.protein_flag = args['protein']
        self.contig_flag = args['contig']
        self.genome_flag = args['genome']
        if self.rbbh_flag:
            self.plot_rbbh()
        if self.ss_flag:
            self.plot_ss()

    def plot_rbbh(self):
        xy_by_data_id = {}
        out_f = 'rbbh.pdf'
        x_col, x_label, x_log = 3, '', True
        y_col, y_label, y_log = 2, '', False
        for infile in os.listdir(self.input_dir):
            if infile.endswith("rbbh.txt"):
                data_id = ".".join(infile.split(".")[0:5])
                rbbh_f = os.path.join(self.input_dir, infile)
                x, y = [], []
                for line in read_file(rbbh_f):
                    col = line.split()
                    if line.startswith("#"):
                        x_label = col[x_col]
                        y_label = col[y_col]
                    else:
                        x.append(float(col[x_col]))
                        y.append(float(col[y_col]))
                xy_by_data_id[data_id] = (x, y)
        PlotObj(xy_by_data_id, x_label, y_label, y_log, x_log, out_f, 'scatter')


class PlotObj():
    def __init__(self, xy_by_data_id, x_label, y_label, y_log, x_log, out_f, plot_type):
        self.xy_by_data_id = xy_by_data_id
        self.x_label = x_label
        self.y_label = y_label
        self.y_log = y_log
        self.x_log = x_log
        self.out_f = out_f
        self.plot_type = plot_type
        self.plot()

    def plot(self):
        f, ax = plt.subplots(figsize=(10, 5))
        if self.plot_type == 'scatter':
            for data_id, xy in self.xy_by_data_id.items():
                x, y = xy
                ax.scatter(x, y, label=data_id, marker='o', alpha=1, s=5)
        elif self.plot_type == 'hist':
            for data_id, xy in self.xy_by_data_id.items():
                x = xy
                ax.hist(x_values, histtype='stepfilled', align='mid', bins=np.arange(0.0, 1.0 + 0.1, 0.1))
        else:
            sys.exit("[X] Plot '%s' not implemented." % (self.plot_type))
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        if self.y_log:
            ax.set_yscale('log')
        if self.x_log:
            ax.set_xscale('log')
        print "[+] - Plotting %s" % (self.out_f)
        ax.legend()
        f.savefig(self.out_f, format="pdf")
        plt.close()

if __name__ == "__main__":
    __version__ = 0.1

    NUCCOUNT = 2
    EVALUE = 1e-25
    QCOV = 70
    args = docopt(__doc__)
    analysis_flag = args['analyse']
    rbbh_flag = args['rbbh']
    plot_flag = args['plot']
    input_dir = args['--input_dir']
    outprefix = args['--outprefix']
    if analysis_flag:
        genomeObjs = generate_genomeObjs(input_dir)
        analysisCollection = AnalysisCollection(genomeObjs)
        analysisCollection.analyse_genomes('splicesites')
        if rbbh_flag:
            analysisCollection.analyse_genomes('rbbh')
    if plot_flag:
        plottingCollection = PlottingCollection(args)


