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

import networkx as nx
import numpy as np
import sys
import logging
import utils
import multinomial_em as mem

oldsettings = np.seterr(all='raise')


def get_biconnected(G):
    """
    Wrapper arround the networkx biconnected_components function. To find out
    why the biconnected components algorithm is useful for finding
    constitutive exons check the information section or wikipedia.
    """

    G_undirected = G.to_undirected()  # make sure undirected graph for biconnected components
    components = filter(lambda x: len(
        x) > 2, map(list, nx.biconnected_components(G_undirected)))  # filter out trivial dyad biconnected components

    # assert len(components) > 0, 'what nothing in it' + str(components)
    # assert components != None, 'Oddly there is a none object in the biconnected comp' + str(components)
    return components


def bellman_ford_longest_path(G, num_nodes, visited, weight='weight'):
    """
    Computes the longest path (most total weight) by only considering
    unexplained edges. That is weights of any edge already in an isoform is
    set to zero. This function tries to minimize the number of isoforms that
    could possibly be generated based on novel edges. Assumes topologically
    sorted with source node as first in topological sort. Topological sort
    version runs in O(n+m) instead of O(nm).
    """
    # initialize variables
    sorted_nodes = sorted(G.nodes())
    d = {nde: float('-inf') for nde in sorted_nodes}
    d[sorted_nodes[0]] = 0   # initialize source to have 0 distance
    p = {nde: [] for nde in sorted_nodes}
    p[sorted_nodes[0]] = [sorted_nodes[0]]  # initialize source path to be it's self

    # "edge relax"
    for tail_node in sorted_nodes:
        for head_node in G.successors(tail_node):
            # want longest path of unexplained edges, so all explained edges have zero weight
            edge_weight = G[tail_node][head_node][
                weight] if visited[tail_node][head_node] == 0 else 0

            # larger total weight case
            if d[head_node] < d[tail_node] + edge_weight:
                d[head_node] = d[tail_node] + edge_weight
                p[head_node] = p[tail_node] + [head_node]
            # same total weight case, choose edge with greater weight into head node
            elif d[head_node] == (d[tail_node] + edge_weight) and G[tail_node][head_node][weight] > G[p[head_node][-2]][head_node][weight]:
                d[head_node] = d[tail_node] + edge_weight
                p[head_node] = p[tail_node] + [head_node]

    longest_path = p[sorted_nodes[-1]]
    no_newly_visited = True
    for i in range(len(longest_path)-1):
        no_newly_visited &= (visited[longest_path[i]][longest_path[i+1]])
    #no_newly_visited = reduce(
        #lambda x, y: x and (visited[x][y] == 1), longest_path)

    return p[sorted_nodes[-1]], no_newly_visited


class AllPaths(object):
    '''
    Handle all possible paths in a biconnected component
    '''
    def __init__(self, sg, component, target, chr=None, strand=None):
        self.set_splice_graph(sg, component, target)
        self.asm_component = self.component  # save the ASM components, self.components may get trimmed if primers aren't placed on first and last exon
        self.chr = chr
        self.strand = strand

    def set_chr(self, chr):
        '''Chromosome setter'''
        self.chr = chr

    def set_strand(self, strand):
        '''Strand setter'''
        if strand == '+' or strand == '-':
            self.strand = strand
        else:
            raise ValueError('Strand should either be + or -')

    def set_splice_graph_old(self, sg, component, target):
        self.graph = sg.get_graph()
        self.tx_paths = sg.annotation
        self.original_tx_paths = sg.annotation  # tx paths all ways without trimming
        known_edges = set([(tx[i], tx[i + 1])
                           for tx in self.tx_paths
                           for i in range(len(tx) - 1)])
        self.component = component
        self.target = target
        self.sub_graph = nx.subgraph(self.graph, self.component)

        # add any possible tx that uses novel edges to list of known txs
        for tx in nx.all_simple_paths(self.sub_graph,
                                      source=self.component[0],
                                      target=self.component[-1]):
            novel = False
            for i in range(len(tx) - 1):
                if (tx[i], tx[i + 1]) not in known_edges:
                    novel = True
            if novel:
                self.tx_paths.append(tx)

        self.inc_lengths, self.skip_lengths = [], []  # call set all_path_lengths method
        self.all_path_coordinates = []  # call set_all_path_coordinates method

    def set_splice_graph(self, sg, component, target):
        """Setter for the splice graph which consists of graph/transcript attributes"""
        self.graph = sg.get_graph()
        self.tx_paths = sg.annotation
        self.original_tx_paths = sg.annotation  # tx paths all ways without trimming
        self.component = component
        self.target = target
        self.sub_graph = nx.subgraph(self.graph, self.component)

        # add novel txs
        novel_txs = self.all_paths_with_novel_junctions()
        self.tx_paths += novel_txs

    def all_paths_with_novel_junctions(self):
        """
        Create novel isoforms by finding all possible paths that include
        at least one novel junction.
        """
        iter_limit = 10000  # explicitly set error for upper limit of all paths
        tx_list = []  # store novel isforms
        known_edges = set([(tx[i], tx[i + 1])
                          for tx in self.tx_paths
                          for i in range(len(tx) - 1)])
        tmp_sub_graph = self.add_dummy_nodes_to_graph(self.sub_graph)
        tmp_src, tmp_sink = (float('-inf'), float('-inf')), (float('inf'), float('inf'))

        # only add paths that include a novel jct
        for l, tx in enumerate(nx.all_simple_paths(tmp_sub_graph,
                                                   source=tmp_src,
                                                   target=tmp_sink)):
            # set the maximum number of iterations
            if l >= iter_limit:
                raise utils.PrimerSeqError('Iteration limit reached in all paths algorithm.')

            # find all paths with novel edge
            novel = False
            # for i in range(len(tx) - 1):
            for i in range(1, len(tx) - 2):
                if (tx[i], tx[i + 1]) not in known_edges:
                    novel = True
            if novel:
                tx_list.append(tx[1:-1])
        return tx_list

    def add_dummy_nodes_to_graph(self, my_graph):
        """Add a sink and source node to a graph"""
        # setup variables
        tmp_graph = my_graph.copy()  # make sure not editing the original graph
        src_connected_nodes = [node for node in tmp_graph.nodes() if not tmp_graph.predecessors(node)]
        sink_connected_nodes = [node for node in tmp_graph.nodes() if not tmp_graph.successors(node)]

        # add edges to source and sink node
        src, sink = (float('-inf'), float('-inf')), (float('inf'), float('inf'))
        for node in src_connected_nodes:
            tmp_graph.add_edge(src, node)
        for node in sink_connected_nodes:
            tmp_graph.add_edge(node, sink)
        return tmp_graph

    def trim_tx_paths_old(self):
        '''
        Remove all exons outside the biconnected component.
        '''
        self.component = sorted(self.component, key=lambda x: (x[0], x[1]))  # make sure it is sorted

        # trim tx_paths to only contain paths within component_subgraph
        tmp = set()
        for p in self.tx_paths:
            # make sure this tx path has the biconnected component
            if self.component[0] in p and self.component[-1] in p:
                tmp.add(tuple(
                    p[p.index(self.component[0]):p.index(self.component[-1]) + 1]))  # make sure there is no redundant paths
        self.tx_paths = sorted(list(tmp), key=lambda x: (x[0], x[1]))

    def trim_tx_paths(self):
        '''
        Remove all exons outside the biconnected component.
        '''
        self.component = sorted(self.component, key=lambda x: (x[0], x[1]))  # make sure it is sorted

        # trim tx_paths to only contain paths within component_subgraph
        tmp = set()
        for p in self.tx_paths:
            # make sure this tx path has the biconnected component
            tmp_path = self._get_sub_tx(p)
            if len(tmp_path) > 1:
                tmp.add(tuple(
                    tmp_path))
                    # p[p.index(self.component[0]):p.index(self.component[-1]) + 1]))  # make sure there is no redundant paths
        self.tx_paths = sorted(list(tmp), key=lambda x: (x[0], x[1]))

    def _get_sub_tx(self, path):
        sub_path = []
        for p in path:
            if p in self.component:
                if sub_path and not self.sub_graph.has_edge(sub_path[-1], p):
                    return []
                sub_path.append(p)
            elif len(sub_path) > 1:
                return sub_path
        return sub_path

    def trim_tx_paths_using_primers(self, first_primer, second_primer, first_exon, second_exon):
        """
        Get rid of all transcripts which do not contain both the first
        and second primer. This method may keep transcripts that do not
        contain the exact user defined flanking exons.
        """
        self.component = sorted(self.component, key=lambda x: (x[0], x[1]))  # make sure it is sorted

        # trim tx_paths to only contain paths within component_subgraph
        tmp = set()
        for p in self.tx_paths:
            first_ex = utils.find_first_exon(first_primer, p)
            last_ex = utils.find_last_exon(second_primer, p)
            if first_ex is not None and last_ex is not None:
                tmp_path = p[first_ex:last_ex+1]
                if tmp_path[0][0] < first_exon[0]:
                    tmp_path[0] = (first_exon[0], tmp_path[0][1])  # don't be before user-defined exon
                if tmp_path[-1][1] > second_exon[1]:
                    tmp_path[-1] = (tmp_path[-1][0], second_exon[1])  # don't be after user-defined exon
                tmp.add(tuple(tmp_path))  # make sure no redundancies
        # self.tx_paths = sorted(list(tmp), key=lambda x: (x[0], x[1]))
        self.tx_paths = list(tmp)

    def trim_tx_paths_using_flanking_exons_and_target(self, strand, 
                                                      target_exon, up_exon, down_exon):
        tmp = set()
        for p in self.tx_paths:
            # make sure this tx path has the biconnected component
            flank_exon_flag = (up_exon in p and down_exon in p)
            target_exon_flag = target_exon in p
            if flank_exon_flag:
                if strand == '+':
                    first_index, second_index = p.index(up_exon), p.index(down_exon)
                elif strand == '-':
                    first_index, second_index = p.index(down_exon), p.index(up_exon)
                tmp.add(tuple(
                    sorted(p[first_index:second_index + 1], key=lambda x: (x[0], x[1]))))  # make sure there is no redundant paths
            elif target_exon_flag:
                tmp_p = []
                for ex in p:
                    if strand == "+":
                        if (up_exon[0] <= ex[0] <= down_exon[1]) or (up_exon[0] <= ex[1] <= down_exon[1]):
                            tmp_p.append(ex)
                    elif strand == "-":
                        if (down_exon[0] <= ex[0] <= up_exon[1]) or (down_exon[0] <= ex[1] <= up_exon[1]):
                            tmp_p.append(ex)
                tmp.add(tuple(sorted(tmp_p, key=lambda x: (x[0], x[1]))))

        # self.tx_paths = sorted(list(tmp), key=lambda x: (x[0], x[1]))
        self.tx_paths = list(tmp)

    def trim_tx_paths_using_flanking_exons(self, strand, up_exon, down_exon):
        tmp = set()
        for p in self.tx_paths:
            # make sure this tx path has the biconnected component
            if up_exon in p and down_exon in p:
                if strand == '+':
                    first_index, second_index = p.index(up_exon), p.index(down_exon)
                elif strand == '-':
                    first_index, second_index = p.index(down_exon), p.index(up_exon)
                tmp.add(tuple(
                    sorted(p[first_index:second_index + 1], key=lambda x: (x[0], x[1]))))  # make sure there is no redundant paths
        # self.tx_paths = sorted(list(tmp), key=lambda x: (x[0], x[1]))
        self.tx_paths = list(tmp)

    def trim_tx_paths_using_flanking_exons2(self, strand, up_exon, down_exon):
        """Keep TXs with appropriate flanking exons"""
        tmp = set()
        for p in self.tx_paths:
            tmp_starts, tmp_ends = zip(*p)  # get starts and ends
            if strand == "+":
                my_flag = up_exon[1] in tmp_ends and down_exon[0] in tmp_starts
            elif strand == '-':
                my_flag = up_exon[0] in tmp_starts and down_exon[1] in tmp_ends
            if my_flag:
                if strand == '+':
                    first_index = tmp_ends.index(up_exon[1])
                    second_index = tmp_starts.index(down_exon[0])
                elif strand == '-':
                    first_index = tmp_ends.index(down_exon[1])
                    second_index = tmp_starts.index(up_exon[0])
                tmp.add(tuple(
                    p[first_index:second_index + 1]))  # make sure there is no redundant paths
            my_flag = False
        self.tx_paths = sorted(list(tmp), key=lambda x: (x[0], x[1]))

    def keep_weakly_connected(self):
        '''This method filters out exons (nodes) not involved in AS events'''
        # find weakly connected subgraphs
        weakly_connected_list = nx.weakly_connected_component_subgraphs(self.sub_graph)

        # iterate to find which subgraph has the target exon
        for subgraph in weakly_connected_list:
            if self.target in subgraph.nodes():
                self.sub_graph = subgraph  # assign subgraph that actually connects to target exon

    def estimate_counts(self):
        '''
        Estimates read counts by using :func:`~algorithms.read_count_em`
        and then returns the transcript paths and read counts for those
        paths.
        '''
        # check the connectivity of the graph -- deprecated checking
        # if not nx.is_weakly_connected(self.sub_graph): raise utils.PrimerSeqError('Error: SpliceGraph should be connected')
        # make sure graph (self.sub_graph) is weakly connected
        tmp_sub_graph = self.sub_graph.copy()
        self.keep_weakly_connected()
        if not len(self.sub_graph.nodes()) > 1:
            self.sub_graph = tmp_sub_graph

        # AFE/ALE testing
        num_first_exons = len(filter(lambda x: len(self.sub_graph.predecessors(x)) == 0, self.sub_graph.nodes()))
        if num_first_exons > 1: utils.PrimerSeqError('Error: not internal AS event')
        num_last_exons = len(filter(lambda x: len(self.sub_graph.successors(x)) == 0, self.sub_graph.nodes()))
        if num_last_exons > 1: utils.PrimerSeqError('Error: not internal AS event')

        # run EM algorithm
        logging.debug('Start read count EM algorithm . . . ')
        self.count_info = mem.multinomial_em(self.tx_paths, self.sub_graph)
        logging.debug('Finished calculating counts.')

        return map(list, self.tx_paths), self.count_info

    def set_all_path_coordinates(self):
        '''
        Computes the coordinates ofr each tx
        '''
        tmp = []
        for p in self.tx_paths:
            # tmp.append(map(lambda x: (self.strand, self.chr, x[0], x[1]), self.tx_paths))
            tmp.append(map(lambda x: (self.strand, self.chr, x[0], x[1]), p))
        self.all_path_lengths = tmp

    def set_all_path_lengths(self, primer_coords):
        '''
        Computes the path length for each isoform
        '''
        # get possible lengths
        inc_length, skip_length = [], []
        for path in self.original_tx_paths:
            if self.target in path:
                inc_length.append(utils.calc_product_length(path, primer_coords))  # length of everything but target exon and flanking constitutive exons
            else:
                skip_length.append(utils.calc_product_length(path, primer_coords))  # length of everything but target exon and flanking constitutive exons
        self.inc_lengths, self.skip_lengths = list(set(inc_length)), list(set(skip_length))

    def get_shortest_path(self):
        """
        Returns the shortest isoform, this could be a skipping isoform or an
        inclusion isoform.
        """
        min_len = float('inf')
        shortest_tx = []
        for tx in self.tx_paths:
            tx_len = sum([end - start for start, end in tx[1:-1]])  # this line won't work for retained introns
            if tx_len < min_len:
                shortest_tx = tx
                min_len = tx_len
        return shortest_tx
