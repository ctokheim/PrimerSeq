import csv
import re
import networkx as nx
import itertools as it
import argparse
import json
import algorithms as algs
from pygr.seqdb import SequenceFileDB
import matplotlib.pyplot as plt
import gtf
from utils import get_chr, get_start_pos, get_end_pos, get_pos
import sam
import sys
from bed import Bed
from wig import Wig
from exon_seek import ExonSeek

# logging imports
import logging


def gene_annotation_reader(file_path, FILTER_FACTOR=2):
    """
    *creates two data structures from gtf:*

        gene_dict is a dictionary with gene_id's as keys:
            * gene_dict['My_favorite_gene']['graph'] = list of list of exons
            * gene_dict['My_favorite_gene']['chr'] = gene chromosome (chrXX)
            * gene_dict['My_favorite_gene']['strand'] = gene strand (+, -)
            * gene_dict['My_favorite_gene']['start'] = start of gene
            * gene_dict['My_favorite_gene']['end'] = end of gene
            * gene_dict['My_favorite_gene']['exons'] = the set of exons (nodes)

        gene_lookup is a dictionary to lookup the gene that an exon coordinate belongs to:
            * gene_lookup['(strand)(chr):(start)-(end)'] = ['My_favorite_gene', ...]
    """
    # iterate through each gtf feature
    file_input = open(file_path)
    gene_dict, gene_lookup = {}, {}
    for tx_id, tx in it.groupby(gtf.gtf_reader(file_input, delim='\t'), lambda x: x.attribute['transcript_id']):
        tx = list(tx)
        if len(tx) == 0: continue  # no 'exon' feature case

        gene_id = tx[0].attribute['gene_id']
        strand = tx[0].strand

        # sort exons
        tx_path = sorted([(exon.start, exon.end) for exon in tx],
                         key=lambda x: (x[0], x[1]))  # needs to be sorted because gtf files might not have them in proper order

        # add info to gene_dict
        gene_dict.setdefault(gene_id, {})  # add the gene key if it doesn't exist
        gene_dict[gene_id].setdefault('chr', tx[0].seqname)  # add chr if doesn't exist
        gene_dict[gene_id].setdefault('strand', strand)  # add strand if doesn't exist
        gene_dict[gene_id].setdefault('graph', []).append(tx_path)  # append the tx path
        gene_dict[gene_id].setdefault('start', float('inf'))
        gene_dict[gene_id]['start'] = min(gene_dict[gene_id]['start'], tx_path[0][0])  # change start if this tx has lower start position
        gene_dict[gene_id].setdefault('end', 0)
        gene_dict[gene_id]['end'] = max(gene_dict[gene_id]['end'], tx_path[-1][1])
        gene_dict[gene_id].setdefault('exons', set())
        for ex in tx_path:
            gene_dict[gene_id]['exons'].add(ex)  # hold a set of non-redundant exons

        # lookup dict for which exon is in which gene
        for exon in tx:
            pos_string = '%s%s:%d-%d' % (strand, exon.seqname, exon.start, exon.end)
            gene_lookup.setdefault(pos_string, []).append(gene_id)

    file_input.close()
    return gene_dict, gene_lookup


class SpliceGraph(object):
    '''
    The SpliceGraph class is meant to provide a common repository for creation
    of splice graphs for a single gene.
    '''

    def __init__(self, annotation, chr, strand, read_threshold=5):
        self.chr = chr
        self.strand = strand
        self.READ_THRESHOLD = read_threshold
        self.annotation = []  # set value using set_graph_as_annotation
        if annotation is not None:
            self.set_graph_as_annotation(annotation)
        else:
            self.graph = None

    def get_graph(self):
        return self.graph

    def set_graph_as_annotation(self, annotation, FILTER_FACTOR=2):
        """
        Create a nx DiGraph from list of tx in gene. FILTER_FACTOR defines a cutoff for
        using a tx of a gene. A tx must have x num of exon where x > MAX tx exon num / FILTER_FACTOR.
        """
        # filter low exon tx
        max_exons = max(map(len, annotation))  # figure out max num exons
        self.annotation = map(lambda y: sorted(y, key=lambda z:(z[0], z[1])),  # make sure exons are sorted by position
                              filter(lambda x: len(x) > max_exons / FILTER_FACTOR, annotation))  # filter based on max num exons criteria

        # create graph
        graph = nx.DiGraph()
        for tx in self.annotation:
            graph.add_path(tx)
        self.graph = graph  # set graph attribute

    def set_annotation_edge_weights(self, weights):
        """
        Only try to find weights for already existing edges in the graph.
        This function is intended to add weight values to edges defined
        in the gtf annotation.
        """
        # add edge weights to edges from annotation
        for u, v in self.graph.edges():
            try:
                #tmpChr = get_chr(exon_forms[u])
                start = u[1]  # get_start_pos(exon_forms[u])
                end = v[0]  # get_end_pos(exon_forms[v])
                tmpWeight = weights[self.chr][start][end]
                self.graph[u][v]['weight'] = tmpWeight
            except KeyError:
                self.graph[u][v]['weight'] = 1

    def set_graph_as_nodes_only(self, exons):
        """
        Simple function that makes a DAG (nx.DiGraph) with only nodes and no edges.
        Meant to be used to before add_all_possible_edge_weights.
        """
        G = nx.DiGraph()
        G.add_nodes_from(exons)
        self.graph = G

    def add_all_possible_edge_weights(self, weights):  # use to have exon_forms rather than chr
        """
        Add edge/weights to graph if supported by atleast READ_THRESHOLD number of reads
        """
        # add novel edges if well supported
        sorted_nodes = sorted(self.graph.nodes())
        for i in range(len(sorted_nodes) - 1):
            for j in range(i + 1, len(sorted_nodes)):
                try:
                    # tmpChr = get_chr(exon_forms[sorted_nodes[i]])
                    start = sorted_nodes[i][1]   # get_start_pos(exon_forms[sorted_nodes[i]])
                    end = sorted_nodes[j][0]    # get_end_pos(exon_forms[sorted_nodes[j]])
                    if weights[self.chr][start][end] >= self.READ_THRESHOLD:
                        self.graph.add_edge(sorted_nodes[i], sorted_nodes[j])
                        self.graph[sorted_nodes[i]][sorted_nodes[
                            j]]['weight'] = weights[self.chr][start][end]
                except KeyError:
                    pass


def get_flanking_biconnected_exons(name, target, sGraph, genome):
    graph = sGraph.get_graph()  # nx.DiGraph
    # search through each biconnected component
    for component in algs.get_biconnected(graph):
        component = sorted(component, key=lambda x: (x[0], x[1]))  # ensure first component is first exon, etc
        if target in component[1:-1]:
            # define upstream/downstream flanking exon
            if sGraph.strand == '+':
                upstream = component[0]
                downstream = component[-1]
            else:
                upstream = component[-1]
                downstream = component[0]

            # get possible lengths
            all_paths = algs.AllPaths(graph, component, target,
                                      chr=sGraph.chr, strand=sGraph.strand)
            inc_length, skip_length = all_paths.get_all_path_lengths()

            # get sequence of upstream/target/downstream combo
            genome_chr = genome[sGraph.chr]  # chr object from pygr
            upstream_seq, target_seq, downstream_seq = genome_chr[upstream[0]:upstream[1]], genome_chr[target[0]:target[1]], genome_chr[downstream[0]:downstream[1]]
            if sGraph.strand == '-':
                upstream_seq, target_seq, downstream_seq =  \
                    -upstream_seq, -target_seq, -downstream_seq

            return [sGraph.strand, name[1:],
                    sGraph.chr + ':' + '-'.join(map(str, upstream)),
                    sGraph.chr + ':' + '-'.join(map(str, downstream)),
                    inc_length, skip_length, str(upstream_seq).upper(),
                    str(target_seq).upper(), str(downstream_seq).upper()]
    return [name + ' was not found in a biconnected component']


def get_sufficient_psi_exons(name, target, sGraph, genome):
    # find appropriate flanking "constitutive" exon for primers
    # upstream, downstream, component, (psi_target, psi_upstream, psi_downstream) = find_fuzzy_constitutive(target, sGraph)
    exon_seek_obj = ExonSeek(target, sGraph)
    upstream, downstream, component, psi_target, psi_upstream, psi_downstream = exon_seek_obj.get_info()

    # lack of successor/predecessor nodes
    if upstream is None or downstream is None:
        logging.debug("%s does not have an upstream exon, downstream exon, or possibly both" % str(component))
        return ["%s does not have an upstream exon, downstream exon, or possibly both" % str(component)]

    # get possible lengths
    all_paths = algs.AllPaths(sGraph.get_graph(), component, target,
                              chr=sGraph.chr, strand=sGraph.strand)
    inc_length, skip_length = all_paths.get_all_path_lengths()

    # get sequence of upstream/target/downstream combo
    genome_chr = genome[sGraph.chr]  # chr object from pygr
    upstream_seq, target_seq, downstream_seq = genome_chr[upstream[0]:upstream[1]], genome_chr[target[0]:target[1]], genome_chr[downstream[0]:downstream[1]]
    if sGraph.strand == '-':
        upstream_seq, target_seq, downstream_seq =  \
            -upstream_seq, -target_seq, -downstream_seq

    return [sGraph.strand, name[1:], psi_target,
            sGraph.chr + ':' + '-'.join(map(str, upstream)), psi_upstream,
            sGraph.chr + ':' + '-'.join(map(str, downstream)), psi_downstream,
            inc_length, skip_length, str(upstream_seq).upper(),
            str(target_seq).upper(), str(downstream_seq).upper()]


def main(options, args_output='tmp/debug.json'):
    """
    The gtf main function is the function designed to be called from other scripts. It iterates through each target
    exons and returns the necessary information for primer design.
    """
    fasta, args_big_bed, args_target = options['fasta'], options['big_bed'], options['target']
    genome = SequenceFileDB(fasta)  # load fasta file from file
    bed = Bed(args_big_bed, ext='bed')

    # the sam object interfaces with the user specified BAM/SAM file!!!
    sam_obj = sam.Sam(options['rnaseq'])

    # iterate through each target exon
    output = []  # output from program
    for line in args_target:  # was line in handle
        line = line.strip()
        strand = line[0]
        tmp_start, tmp_end = get_pos(line)
        chr = get_chr(line[1:])

        # This try block is to catch assertions made about the graph. If an
        # assertion is raised it only impacts a single target for primer design
        # so complete exiting of the program is not warranted.
        try:
            # get gene annotation from bigBed file
            bed.extractBigRegion(strand, chr, tmp_start, tmp_end)
            bed.load_bed_file()
            gene_dict = bed.get_annotation()

            if options['both_flag']:
                splice_graph = SpliceGraph(annotation=gene_dict['graph'],  # use junctions from annotation
                                           chr=chr,
                                           strand=strand,
                                           read_threshold=options['read_threshold'])
                edge_weights = sam_obj.extractSamRegion(chr, gene_dict['start'], gene_dict['end'])
                splice_graph.set_annotation_edge_weights(edge_weights)  # set edge weights supported from annotation
                splice_graph.add_all_possible_edge_weights(edge_weights)  # also use junctions from RNA-Seq
            elif options['annotation_flag']:
                splice_graph = SpliceGraph(annotation=gene_dict['graph'],  # use annotation
                                           chr=chr,
                                           strand=strand,
                                           read_threshold=options['read_threshold'])
            elif options['rnaseq_flag']:
                splice_graph = SpliceGraph(annotation=None,  # set to None to not use annotation
                                           chr=chr,
                                           strand=strand,
                                           read_threshold=options['read_threshold'])
                splice_graph.set_graph_as_nodes_only(list(gene_dict['exons']))
                edge_weights = sam_obj.extractSamRegion(chr, gene_dict['start'], gene_dict['end'])
                splice_graph.add_all_possible_edge_weights(edge_weights)  # add all edges supported by rna seq

            # default case
            if options['psi'] > .9999:
                tmp = get_flanking_biconnected_exons(line, gene_dict['target'],
                                                     splice_graph,
                                                     genome)  # note this function ignores edge weights
            # user specified a sufficient psi value to call constitutive exons
            else:
                tmp = get_sufficient_psi_exons(line, gene_dict['target'],
                                               splice_graph,
                                               genome)
            output.append(tmp)
        except AssertionError:
            t, v, trace = sys.exc_info()
            output.append([str(v)])  # just append assertion msg

    # if they specify a output file then write to it
    if args_output:
        with open(args_output, 'wb') as write_handle:
            json.dump(output, write_handle, indent=4)

    return output


if __name__ == '__main__':
    """Running this script directly is only for debug purposes"""
    # process command line arguments
    parser = argparse.ArgumentParser(description='Get flanking constitutive exons')
    parser.add_argument('-b', '--big-bed', action='store', dest='big_bed', required=True,
                        help='annotation file with legitimate gene_id\'s')
    parser.add_argument('-t', '--target', action='store', dest='target', required=True,
                        help='file of list of coordinate targets')
    parser.add_argument('-f', '--fasta', action='store', dest='fasta', required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--annotaton', dest='annotation_flag', action='store_true')
    group.add_argument('--rnaseq', dest='rnaseq_flag', action='store_true')
    group.add_argument('--both', dest='both_flag', action='store_true')
    parser.add_argument('--psi', dest='psi', action='store', type=float)
    parser.add_argument('--read-threshold', dest='read_threshold', type=int, action='store')
    parser.add_argument('-o', '--output', action='store', dest='output', required=True)
    options = vars(parser.parse_args())

    options['target'] = options['target'].replace('dash', '-').split(',')  # fix bug with - as input for strand

    # call main function
    main(options, options['output'])
