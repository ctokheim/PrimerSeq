#!/usr/bin/env python
# Copyright (C) 2012-2013  Collin Tokheim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

'''
**Author:** Collin Tokheim

Splice Graph
------------

The splice_graph.py deals with gene structure as a (weighted) directed acyclic
graph (DAG) normally known as a splice graph. Like you would expect, the
:class:`~splice_graph.SpliceGraph` class represents gene structure as a splice
graph.

Transcripts Overlapping Target
------------------------------

:func:`~splice_graph.get_from_gtf_using_gene_name` returns the all transcripts
    from a gene if the gene overlaps the user's target.
:func:`~splice_graph.get_weakly_connected_tx` returns all transcripts that are
    weakly connected to the user's target. This approach is used when gene IDs
    are not valid and thus can not be used.

Flanking Exons
--------------

If read counts are used (psi < 1) then splice_graph uses the
:func:`~splice_graph.get_sufficient_psi_exons` function to find flanking exons
    with at least a user defined inclusion level. If psi == 1 then flanking
    exons are only determined by the biconnected componenets algorithm using
:func:`~splice_graph.get_flanking_biconnected_exons`.
'''

import networkx as nx
import itertools as it
import argparse
import algorithms as algs
import gtf
from utils import get_chr, get_start_pos, get_end_pos, get_pos, merge_list_of_dicts
import utils
import sys
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
        Create a nx DiGraph from list of tx in gene. FILTER_FACTOR defines a
        cutoff for using a tx of a gene. A tx must have x num of exon where x
        > MAX tx exon num / FILTER_FACTOR.
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
        Only try to find weights for already existing edges in the graph. This
        function is intended to add weight values to edges defined in the gtf
        annotation.
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
        Simple function that makes a DAG (nx.DiGraph) with only nodes and no
        edges. Meant to be used to before add_all_possible_edge_weights.
        """
        G = nx.DiGraph()
        G.add_nodes_from(exons)
        self.graph = G

    def add_all_possible_edge_weights(self, weights):  # use to have exon_forms rather than chr
        """
        Add edge/weights to graph if supported by atleast READ_THRESHOLD
        number of reads
        """
        # add novel edges if well supported
        sorted_nodes = sorted(self.graph.nodes())
        for i in range(len(sorted_nodes) - 1):
            for j in range(i + 1, len(sorted_nodes)):
                try:
                    start = sorted_nodes[i][1]   # get_start_pos(exon_forms[sorted_nodes[i]])
                    end = sorted_nodes[j][0]    # get_end_pos(exon_forms[sorted_nodes[j]])
                    if weights[(self.chr, start, end)] >= self.READ_THRESHOLD:
                        self.graph.add_edge(sorted_nodes[i], sorted_nodes[j])
                        self.graph[sorted_nodes[i]][sorted_nodes[
                            j]]['weight'] = weights[(self.chr, start, end)]
                except KeyError:
                    pass


def get_from_gtf_using_gene_name(gtf, strand, chr, start, end):
    '''
    This function finds the first gene in the gtf that completely contains the
    target interval. I should really think about checking for multiple genes
    instead of just returning.
    '''
    for gene_key in gtf[chr]:
        if gtf[chr][gene_key]['strand'] == strand and gtf[chr][gene_key]['start'] <= start and gtf[chr][gene_key]['end'] >= end:
            for ex in gtf[chr][gene_key]['exons']:
                # if start >= ex[0] and end <= ex[1]:
                if start == ex[0] and end == ex[1]:
                    gtf[chr][gene_key]['target'] = ex  # this line needed for compatability reasons
                    return gtf[chr][gene_key]
    raise utils.PrimerSeqError("Error: Did not find an appropriate gtf annotation")


def get_weakly_connected_tx(gtf, strand, chr, start, end, plus_or_minus=1000000):
    '''
    This function is meant to handle tx annotations without gene ids.
    Currently this is a function outside of the SpliceGraph class but it may
    be beneficial to later include this as a method.
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
    if not (len(weakly_con_subgraphs) > 0): raise utils.PrimerSeqError('Error: No annotations were even near your target')
    target_graph = None
    for weak_subgraph in weakly_con_subgraphs:
        for node_start, node_end in weak_subgraph.nodes():
            # if node_start <= start and node_end >= end:
            if node_start == start and node_end == end:
                target_graph = weak_subgraph
                start, end = node_start, node_end
    if target_graph is None: raise utils.PrimerSeqError('Error: Target was not contained in a tx')

    # filter tmp_tx to tx that contain atleast one node in subgraph
    filtered_tmp_tx = []
    for tx in tmp_tx:
        for exon in tx:
            if exon in target_graph.nodes():
                filtered_tmp_tx.append(tx)
                break
    if not (len(filtered_tmp_tx) > 0): utils.PrimerSeqError('Error: Your target was not contained in a tx.')

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
    '''
    Defines flanking exons as exons that cannot be skipped in
    the graph structure. Theese exons are 100% included and do not
    need estimation of inclusion level.
    '''
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
            all_paths = algs.AllPaths(sGraph, component, target,
                                      chr=sGraph.chr, strand=sGraph.strand)
            # all_paths.set_all_path_lengths()  # should no longer need this since it is done in primer.py
            all_paths.set_all_path_coordinates()

            # get sequence of upstream/target/downstream combo
            genome_chr = genome[sGraph.chr]  # chr object from pygr
            upstream_seq, target_seq, downstream_seq = genome_chr[upstream[0]:upstream[1]], genome_chr[target[0]:target[1]], genome_chr[downstream[0]:downstream[1]]
            if sGraph.strand == '-':
                upstream_seq, target_seq, downstream_seq =  \
                    -upstream_seq, -target_seq, -downstream_seq

            return [sGraph.strand, name[1:], 'NA',
                    sGraph.chr + ':' + '-'.join(map(str, upstream)), '1.0',
                    sGraph.chr + ':' + '-'.join(map(str, downstream)), '1.0',
                    all_paths, str(upstream_seq).upper(),
                    str(target_seq).upper(), str(downstream_seq).upper()]
    return ['Error: ' + name + ' was not found in a biconnected component']


def get_sufficient_psi_exons(name, target, sGraph, genome, ID, cutoff, upstream_exon, downstream_exon):
    """
    Utilizes the ExonSeek class to find flanking exons that are
    good enough to be called "constitutive".
    """
    # find appropriate flanking "constitutive" exon for primers
    exon_seek_obj = ExonSeek(target, sGraph, ID, cutoff, upstream_exon, downstream_exon)
    all_paths, upstream, downstream, component, psi_target, psi_upstream, psi_downstream = exon_seek_obj.get_info()

    # lack of successor/predecessor nodes
    if upstream is None or downstream is None:
        logging.debug("Error: %s does not have an upstream exon, downstream exon, or possibly both" % str(component))
        return ["Error: %s does not have an upstream exon, downstream exon, or possibly both" % str(component)]

    # get sequence of upstream/target/downstream combo
    genome_chr = genome[sGraph.chr]  # chr object from pygr
    upstream_seq, target_seq, downstream_seq = genome_chr[upstream[0]:upstream[1]], genome_chr[target[0]:target[1]], genome_chr[downstream[0]:downstream[1]]
    if sGraph.strand == '-':
        upstream_seq, target_seq, downstream_seq =  \
            -upstream_seq, -target_seq, -downstream_seq

    return [sGraph.strand, name[1:], psi_target,
            sGraph.chr + ':' + '-'.join(map(str, upstream)),  # upstream eg. +chr1:1000-2000
            psi_upstream,
            sGraph.chr + ':' + '-'.join(map(str, downstream)),  # downstream eg. +chr1:1000-2000
            psi_downstream,
            all_paths, upstream_seq,
            target_seq, downstream_seq]


def calculate_target_psi(target, sg_list, component):
    """
    Calculate psi for the target exon for each bam file.
    """
    logging.debug("Calculating psi for each bam file  . . .")
    psi_list = []
    for sg in sg_list:
        ap = algs.AllPaths(sg, component, target, chr=sg.chr)
        ap.trim_tx_paths()
        paths, counts = ap.estimate_counts()
        tmp_inc_count, tmp_skip_count = 0., 0.
        for i, p in enumerate(paths):
            if target in p:
                tmp_inc_count += counts[i] / (len(p) - 1)  # need to normaliz inc counts by number of jcts
            else:
                tmp_skip_count += counts[i] / (len(p) - 1)  # need to normalize skip counts by number of jcts
        psi_list.append(tmp_inc_count / (tmp_inc_count + tmp_skip_count))
    logging.debug("Finished calculating psi for each bam file.")

    return ';'.join(map(lambda x: '%.4f' % x, psi_list))  # only report to four decimal places


def construct_splice_graph(edge_weights_list, gene_dict, chr, strand, read_threshold, min_count,
                           output_type='single', both=False):
    """
    Handles construction of SpliceGraph objects
    """
    if output_type == 'single':
        # case where counts are pooled from all BAM files
        splice_graph = SpliceGraph(annotation=gene_dict['graph'],  # use junctions from annotation
                                   chr=chr,
                                   strand=strand,
                                   read_threshold=read_threshold,
                                   min_count=min_count)
        edge_weights = merge_list_of_dicts(edge_weights_list)  # merge all SAM/BAM read counts to a single dictionary
        splice_graph.set_annotation_edge_weights(edge_weights)  # set edge weights supported from annotation
        if both: splice_graph.add_all_possible_edge_weights(edge_weights)  # also use junctions from RNA-Seq
        return splice_graph
    elif output_type == 'list':
        # returns a list of splice graphs (one for each BAM file)
        single_bam_splice_graphs = []
        for eweight in edge_weights_list:
            tmp_sg = SpliceGraph(annotation=gene_dict['graph'],
                                 chr=chr,
                                 strand=strand,
                                 read_threshold=read_threshold,
                                 min_count=min_count)
            tmp_sg.set_annotation_edge_weights(eweight)
            if both: tmp_sg.add_all_possible_edge_weights(eweight)
            single_bam_splice_graphs.append(tmp_sg)
        return single_bam_splice_graphs


def main(options, args_output='tmp/debug.json'):
    """
    The gtf main function is the function designed to be called from other
    scripts. It iterates through each target exons and returns the necessary
    information for primer design.
    """
    genome, args_gtf, args_target = options['fasta'], options['gtf'], options['target']

    # the sam object interfaces with the user specified BAM/SAM file!!!
    sam_obj_list = options['rnaseq']

    # iterate through each target exon
    output = []  # output from program
    for line in args_target:  # was line in handle
        name, line = line  # bad style of reassignment
        tgt = line[0]
        strand = tgt[0]
        tmp_start, tmp_end = get_pos(tgt)
        chr = get_chr(tgt[1:])  # [1:] since strand is first character
        USER_DEFINED_FLANKING_EXONS = True if len(line) == 3 else False
        if USER_DEFINED_FLANKING_EXONS:
            up_exon = utils.get_pos(line[1])  # user's upstream exon
            down_exon = utils.get_pos(line[2])  # user's downstream exon
        else:
            up_exon = None  # user did not provide upstream exon
            down_exon = None  # user did not provide downstream exon

        # This try block is to catch assertions made about the graph. If a
        # PrimerSeqError is raised it only impacts a single target for primer
        # design so complete exiting of the program is not warranted.
        try:
            # if the gtf doesn't have a valid gene_id attribute then use
            # the first method otherwise use the second method.
            if options['no_gene_id']:
                gene_dict = get_weakly_connected_tx(args_gtf, strand, chr, tmp_start, tmp_end)  # hopefully filter out junk
            else:
                gene_dict = get_from_gtf_using_gene_name(args_gtf, strand, chr, tmp_start, tmp_end)

            # extract all edge weights only once
            edge_weights_list = [sam_obj.extractSamRegion(chr, gene_dict['start'], gene_dict['end'])
                                 for sam_obj in sam_obj_list]

            # The following options['both_flag'] determines how the splice graph is constructed.
            # The splice graph can be either constructed from annotation junctions
            # where options['both_flag']==False or RNA-Seq + annotation junctions when
            # options['both_flag']==True.

            # single pooled count data splice graph
            splice_graph = construct_splice_graph(edge_weights_list,
                                                  gene_dict,
                                                  chr,
                                                  strand,
                                                  options['read_threshold'],
                                                  options['min_jct_count'],
                                                  output_type='single',
                                                  both=options['both_flag'])
            # Second, get a splice graph for each BAM file
            single_bam_splice_graphs = construct_splice_graph(edge_weights_list,
                                                              gene_dict,
                                                              chr,
                                                              strand,
                                                              options['read_threshold'],
                                                              options['min_jct_count'],
                                                              output_type='list',
                                                              both=options['both_flag'])

            # always included case
            if options['psi'] > .9999:
                tmp = get_flanking_biconnected_exons(tgt, gene_dict['target'],
                                                     splice_graph,
                                                     genome)  # note this function ignores edge weights
            # user specified a sufficient psi value to call constitutive exons
            else:
                tmp = get_sufficient_psi_exons(tgt, gene_dict['target'],
                                               splice_graph,
                                               genome,
                                               name,
                                               options['psi'],
                                               up_exon,
                                               down_exon)  # note, this function utilizes edge wieghts

            # Error msgs are of length one, so only do psi calculations for
            # non-error msgs
            if len(tmp) > 1:
                # edit target psi value
                tmp_all_paths = tmp[-4]  # CAREFUL the index for the AllPaths object may change
                tmp[2] = calculate_target_psi(gene_dict['target'], single_bam_splice_graphs, tmp_all_paths.component)  # CAREFUL index for psi_target may change

            # append result to output list
            output.append(tmp)
        except (utils.PrimerSeqError,):
            t, v, trace = sys.exc_info()
            output.append([str(v)])  # just append assertion msg

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
