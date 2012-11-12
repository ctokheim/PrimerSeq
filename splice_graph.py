import csv
import re
import networkx as nx
import itertools as it
import argparse
import json
import algorithms as algs
import matplotlib.pyplot as plt
import gtf
from utils import get_chr, get_start_pos, get_end_pos, get_pos, merge_list_of_dicts
import sys
from bed import Bed
from wig import Wig
from exon_seek import ExonSeek

# logging imports
import logging


class SpliceGraph(object):
    '''
    The SpliceGraph class is meant to provide a common repository for creation
    of splice graphs for a single gene.
    '''

    def __init__(self, annotation, chr, strand, read_threshold=5, filter_factor=2, min_count=1):
        self.chr = chr
        self.strand = strand
        self.READ_THRESHOLD = read_threshold
        self.MIN_COUNT = min_count
        self.FILTER_FACTOR = filter_factor
        self.annotation = []  # set value using set_graph_as_annotation
        if annotation is not None:
            self.set_graph_as_annotation(annotation)
        else:
            self.graph = None

    def get_graph(self):
        return self.graph

    def set_graph_as_annotation(self, annotation):
        """
        Create a nx DiGraph from list of tx in gene. FILTER_FACTOR defines a cutoff for
        using a tx of a gene. A tx must have x num of exon where x > MAX tx exon num / FILTER_FACTOR.
        """
        # filter low exon tx
        max_exons = max(map(len, annotation))  # figure out max num exons
        self.annotation = map(lambda y: sorted(y, key=lambda z: (z[0], z[1])),  # make sure exons are sorted by position
                              filter(lambda x: len(x) > max_exons / self.FILTER_FACTOR, annotation))  # filter based on max num exons criteria

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
                #tmpWeight = weights[self.chr][start][end]
                tmpWeight = weights[(self.chr, start, end)]
                self.graph[u][v]['weight'] = tmpWeight
            except KeyError:
                self.graph[u][v]['weight'] = 1  # set dummy value
            self.graph[u][v]['weight'] = max(self.graph[u][v]['weight'], self.MIN_COUNT)  # set read count to at least a user-defined value

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
                        #self.graph[sorted_nodes[i]][sorted_nodes[
                        #    j]]['weight'] = weights[self.chr][start][end]
                        self.graph[sorted_nodes[i]][sorted_nodes[
                            j]]['weight'] = weights[(self.chr, start, end)]
                except KeyError:
                    pass


def get_from_gtf_using_gene_name(gtf, strand, chr, start, end):
    '''
    This function finds the first gene in the gtf that completely contains the target interval.
    I should really think about checking for multiple genes instead of just returning.
    '''
    for gene_key in gtf[chr]:
        if gtf[chr][gene_key]['strand'] == strand and gtf[chr][gene_key]['start'] <= start and gtf[chr][gene_key]['end'] >= end:
            for ex in gtf[chr][gene_key]['exons']:
                if start >= ex[0] and end <= ex[1]:
                    gtf[chr][gene_key]['target'] = ex  # this line needed for compatability reasons
                    return gtf[chr][gene_key]
    assert 1 == 0, "Did not find an appropriate gtf annotation"  # not the most elegant way to say I don't expect to get to this line


def get_weakly_connected_tx(gtf, strand, chr, start, end, plus_or_minus=1000000):
    '''
    This function is meant to handle tx annotations without gene ids.
    Currently this is a function outside of the SpliceGraph class but it may be beneficial
    to later include this as a method.
    '''
    # compile all tx paths that are reasonably close
    tmp_tx = []
    for gene_key in gtf[chr]:
        if gtf[chr][gene_key]['strand'] == strand and gtf[chr][gene_key]['start'] <= (start + plus_or_minus) and gtf[chr][gene_key]['end'] >= (end - plus_or_minus):
            tmp_tx += gtf[chr][gene_key]['graph']

    # get the weakly connected subgraph that contains the target exon
    sg = SpliceGraph(tmp_tx, chr, strand, filter_factor=1000)
    G = sg.get_graph()
    weakly_con_subgraphs = nx.weakly_connected_component_subgraphs(G)
    assert len(weakly_con_subgraphs) > 0, 'No annotations were even near your target'
    target_graph = None
    for weak_subgraph in weakly_con_subgraphs:
        for node_start, node_end in weak_subgraph.nodes():
            if node_start <= start and node_end >= end:
                target_graph = weak_subgraph
                start, end = node_start, node_end
    assert target_graph is not None, 'Target was not contained in a tx'

    # filter tmp_tx to tx that contain atleast one node in subgraph
    filtered_tmp_tx = []
    for tx in tmp_tx:
        for exon in tx:
            if exon in target_graph.nodes():
                filtered_tmp_tx.append(tx)
                break
    assert len(filtered_tmp_tx) > 0, 'Your target was not contained in a tx.'

    ### convert info to dict ###
    g_dict = {}

    # get a unique set of all exons
    exons = set()
    for t in filtered_tmp_tx:
        exons |= set(t)
    g_dict['exons'] = sorted(exons, key=lambda x: (x[0], x[1]))
    g_dict['start'] = g_dict['exons'][0][0]
    g_dict['end'] = g_dict['exons'][-1][1]
    g_dict['chr'] = chr
    g_dict['graph'] = filtered_tmp_tx
    g_dict['target'] = (start, end)

    return g_dict


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
            all_paths.set_all_path_lengths()
            all_paths.set_all_path_coordinates()

            # get sequence of upstream/target/downstream combo
            genome_chr = genome[sGraph.chr]  # chr object from pygr
            upstream_seq, target_seq, downstream_seq = genome_chr[upstream[0]:upstream[1]], genome_chr[target[0]:target[1]], genome_chr[downstream[0]:downstream[1]]
            if sGraph.strand == '-':
                upstream_seq, target_seq, downstream_seq =  \
                    -upstream_seq, -target_seq, -downstream_seq

            return [sGraph.strand, name[1:],
                    sGraph.chr + ':' + '-'.join(map(str, upstream)),
                    sGraph.chr + ':' + '-'.join(map(str, downstream)),
                    all_paths, str(upstream_seq).upper(),
                    str(target_seq).upper(), str(downstream_seq).upper()]
    return [name + ' was not found in a biconnected component']


def get_sufficient_psi_exons(name, target, sGraph, genome, ID):
    # find appropriate flanking "constitutive" exon for primers
    # upstream, downstream, component, (psi_target, psi_upstream, psi_downstream) = find_fuzzy_constitutive(target, sGraph)
    exon_seek_obj = ExonSeek(target, sGraph, ID)
    all_paths, upstream, downstream, component, psi_target, psi_upstream, psi_downstream = exon_seek_obj.get_info()

    # lack of successor/predecessor nodes
    if upstream is None or downstream is None:
        logging.debug("%s does not have an upstream exon, downstream exon, or possibly both" % str(component))
        return ["%s does not have an upstream exon, downstream exon, or possibly both" % str(component)]

    # get possible lengths
    #all_paths = algs.AllPaths(sGraph.get_graph(), component, target,
    #                          chr=sGraph.chr, strand=sGraph.strand)
    # all_paths.set_all_path_lengths()

    # get sequence of upstream/target/downstream combo
    genome_chr = genome[sGraph.chr]  # chr object from pygr
    upstream_seq, target_seq, downstream_seq = genome_chr[upstream[0]:upstream[1]], genome_chr[target[0]:target[1]], genome_chr[downstream[0]:downstream[1]]
    if sGraph.strand == '-':
        upstream_seq, target_seq, downstream_seq =  \
            -upstream_seq, -target_seq, -downstream_seq

    return [sGraph.strand, name[1:], psi_target,
            sGraph.chr + ':' + '-'.join(map(str, upstream)), psi_upstream,
            sGraph.chr + ':' + '-'.join(map(str, downstream)), psi_downstream,
            all_paths, str(upstream_seq).upper(),
            str(target_seq).upper(), str(downstream_seq).upper()]


def main(options, args_output='tmp/debug.json'):
    """
    The gtf main function is the function designed to be called from other scripts. It iterates through each target
    exons and returns the necessary information for primer design.
    """
    genome, args_big_bed, args_gtf, args_target = options['fasta'], options['big_bed'], options['gtf'], options['target']
    if args_big_bed: bed = Bed(args_big_bed, ext='bed')

    # the sam object interfaces with the user specified BAM/SAM file!!!
    sam_obj_list = options['rnaseq']

    # iterate through each target exon
    output = []  # output from program
    for line in args_target:  # was line in handle
        name, line = line  # .strip().split('\t')
        strand = line[0]
        tmp_start, tmp_end = get_pos(line)
        chr = get_chr(line[1:])

        # This try block is to catch assertions made about the graph. If an
        # assertion is raised it only impacts a single target for primer design
        # so complete exiting of the program is not warranted.
        try:
            # get gene annotation from bigBed (maybe deprecate) or gtf file file
            if args_gtf:
                # if the gtf doesn't have a valid gene_id attribute then use
                # the first method otherwise use the second method.
                if options['no_gene_id']:
                    gene_dict = get_weakly_connected_tx(args_gtf, strand, chr, tmp_start, tmp_end)  # hopefully filter out junk
                else:
                    gene_dict = get_from_gtf_using_gene_name(args_gtf, strand, chr, tmp_start, tmp_end)
            else:
                bed.extractBigRegion(strand, chr, tmp_start, tmp_end)
                bed.load_bed_file()
                gene_dict = bed.get_annotation()

            if options['both_flag']:
                splice_graph = SpliceGraph(annotation=gene_dict['graph'],  # use junctions from annotation
                                           chr=chr,
                                           strand=strand,
                                           read_threshold=options['read_threshold'],
                                           min_count=options['min_jct_count'])
                edge_weights_list = [sam_obj.extractSamRegion(chr, gene_dict['start'], gene_dict['end'])
                                     for sam_obj in sam_obj_list]
                edge_weights = merge_list_of_dicts(edge_weights_list)  # merge all SAM/BAM read counts to a single dictionary
                splice_graph.set_annotation_edge_weights(edge_weights)  # set edge weights supported from annotation
                splice_graph.add_all_possible_edge_weights(edge_weights)  # also use junctions from RNA-Seq
            elif options['annotation_flag']:
                splice_graph = SpliceGraph(annotation=gene_dict['graph'],  # use annotation
                                           chr=chr,
                                           strand=strand,
                                           read_threshold=options['read_threshold'],
                                           min_count=options['min_jct_count'])
                edge_weights_list = [sam_obj.extractSamRegion(chr, gene_dict['start'], gene_dict['end'])
                                     for sam_obj in sam_obj_list]
                edge_weights = merge_list_of_dicts(edge_weights_list)  # merge all SAM/BAM read counts to a single dictionary
                splice_graph.set_annotation_edge_weights(edge_weights)  # set edge weights supported from annotation
            elif options['rnaseq_flag']:
                splice_graph = SpliceGraph(annotation=None,  # set to None to not use annotation
                                           chr=chr,
                                           strand=strand,
                                           read_threshold=options['read_threshold'])
                splice_graph.set_graph_as_nodes_only(list(gene_dict['exons']))
                edge_weights_list = [sam_obj.extractSamRegion(chr, gene_dict['start'], gene_dict['end'])
                                     for sam_obj in sam_obj_list]
                edge_weights = merge_list_of_dicts(edge_weights_list)  # merge all SAM/BAM read counts to a single dictionary
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
                                               genome,
                                               name)
            output.append(tmp)
        except (AssertionError, FloatingPointError):
            t, v, trace = sys.exc_info()
            output.append([str(v)])  # just append assertion msg

    # if they specify a output file then write to it
    #if args_output:
    #    with open(args_output, 'wb') as write_handle:
    #        json.dump(output, write_handle, indent=4)

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
