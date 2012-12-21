#!/usr/bin/env python
# Copyright (C) 2012  Collin Tokheim
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

oldsettings = np.seterr(all='raise')


def get_biconnected(G):
    """
    Wrapper arround the networkx biconnected_components function. To find out why the biconnected
    components algorithm is useful for finding constitutive exons check the information section or
    wikipedia.
    """

    G_undirected = G.to_undirected()  # make sure undirected graph for biconnected components
    components = filter(lambda x: len(
        x) > 2, map(list, nx.biconnected_components(G_undirected)))  # filter out trivial dyad biconnected components

    # assert len(components) > 0, 'what nothing in it' + str(components)
    # assert components != None, 'Oddly there is a none object in the biconnected comp' + str(components)
    return components


def bellman_ford_longest_path(G, num_nodes, visited, weight='weight'):
    """
    Computes the longest path (most total weight) by only considering unexplained edges. That is weights of any edge already in an isoform
    is set to zero. This function tries to minimize the number of isoforms that could possibly be generated based on novel edges.
    Assumes topologically sorted with source node as first in topological sort. Topological sort version runs in O(n+m) instead of O(nm).

    input:
        G - a networkx DAG
        num_nodes - the number of nodes in the entire graph (not just the subgraph)
        weight - the string that represents edge weight for G
    output:
        path - longest path from first node to last
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
        '''
        Chromosome setter
        '''
        self.chr = chr

    def set_strand(self, strand):
        '''
        Strand setter
        '''
        if strand == '+' or strand == '-':
            self.strand = strand
        else:
            raise ValueError('Strand should either be + or -')

    def set_splice_graph(self, sg, component, target):
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
        all_paths = list(nx.all_simple_paths(self.sub_graph, source=self.component[0], target=self.component[-1]))
        for tx in all_paths:
            novel = False
            for i in range(len(tx) - 1):
                if (tx[i], tx[i + 1]) not in known_edges:
                    novel = True
            if novel:
                self.tx_paths.append(tx)

        self.inc_lengths, self.skip_lengths = [], []  # call set all_path_lengths method
        self.all_path_coordinates = []  # call set_all_path_coordinates method

    def trim_tx_paths(self):
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

    def estimate_counts(self):
        '''
        Estimates read counts by using :func:`~algorithms.read_count_em`
        and then returns the transcript paths and read counts for those
        paths.
        '''
        # check the connectivity of the graph
        if not nx.is_weakly_connected(self.sub_graph): raise utils.PrimerSeqError('Error: SpliceGraph should be connected')
        # assert nx.is_biconnected(self.sub_graph.to_undirected()), 'Yikes! expected a biconnected component'

        # AFE/ALE testing
        num_first_exons = len(filter(lambda x: len(self.sub_graph.predecessors(x)) == 0, self.sub_graph.nodes()))
        if num_first_exons > 1: utils.PrimerSeqError('Error: not internal AS event')
        num_last_exons = len(filter(lambda x: len(self.sub_graph.successors(x)) == 0, self.sub_graph.nodes()))
        if num_last_exons > 1: utils.PrimerSeqError('ERror: not internal AS event')

        # run EM algorithm
        logging.debug('Start read count EM algorithm . . . ')
        self.count_info = read_count_em(self.tx_paths, self.sub_graph)
        logging.debug('Finished calculating counts.')

        return map(list, self.tx_paths), self.count_info

    def set_all_path_coordinates(self):
        '''
        Computes the coordinates ofr each tx
        '''
        tmp = []
        for p in self.tx_paths:
            tmp.append(map(lambda x: (self.strand, self.chr, x[0], x[1]), self.tx_paths))
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


def read_count_em(bcc_paths, sub_graph):
    '''
    Estimate multinomial probilities by using an EM algorithm by using
    junction reads.
    '''
    oldsettings = np.seterr(all='raise')  # make sure error is raised instead of numerical warning

    # useful convenience dicts
    indexToEdge = {i: e for i, e in enumerate(sub_graph.edges())}
    edgeToIndex = {e: i for i, e in enumerate(sub_graph.edges())}

    # set up count/tx info variables
    num_tx = len(bcc_paths)
    num_edges = len(sub_graph.edges())
    read_counts = np.zeros(num_edges)
    for i, val in enumerate(read_counts):
        u, v = indexToEdge[i]
        read_counts[i] = sub_graph[u][v]['weight']
    total_counts = np.sum(read_counts)

    # set up the uncommited matrix Y
    Y = np.zeros((num_edges, num_tx))
    for tx_index, path in enumerate(bcc_paths):
        for i in range(len(path) - 1):
            try:
                Y[edgeToIndex[(path[i], path[i + 1])]][tx_index] = read_counts[
                    edgeToIndex[(path[i], path[i + 1])]]
            except:
                print '*' * 10, path[i], path[i + 1], tx_index, edgeToIndex
                raise

    # set up p the probability array
    p = np.ones(num_tx) * 1 / num_tx

    # start EM
    counter, MAX_ITERS = 0, 10000  # these variables impose a max iteration restriction of 10000 on the EM algorithm
    epsilon = float('inf')
    THRESHOLD = .0001
    while epsilon > THRESHOLD and counter < MAX_ITERS:
        # E-step
        for i, row in enumerate(Y):
            try:
                if np.sum(row) != 0 and row.dot(p) != 0:
                    Y[i] = row * p / row.dot(p) * read_counts[i]
            except:
                print counter
                print row
                print p
                print row.dot(p)
                raise

        # M-step
        p_new = np.sum(Y, axis=0) / total_counts

        # convergence variable
        epsilon = np.sum(np.abs(p_new - p))

        p = p_new  # update probabilities
        p = p * (p > THRESHOLD)  # call effectively small probabilities zero to avoid numerical underflow errors
        counter += 1  # increment the iteration counter

    tx_counts = total_counts * p
    return tx_counts


def estimate_psi(exon_of_interest, paths, counts):
    """
    Uses the estimated isoform count information from read_count_em to
    estimate the exon inclusion level (psi). Note that read counts are normalized
    using the number of junctions (edges) for that isoform.
    """
    inc_count, skip_count = 0, 0
    for i, num in enumerate(counts):
        if exon_of_interest in paths[i]:
            inc_count += num / float(len(paths[i]) - 1)  # read counts / number of edges
        else:
            skip_count += num / float(len(paths[i]) - 1)  # read counts / number of edges
    psi = float(inc_count) / (inc_count + skip_count)
    return psi


def generate_isoforms(component_subgraph, s_graph):
    """
    Take in a digraph and output isoforms of its biconnected components (splice modules)

    input
        G - a directed acyclic graph
        tx_paths - a list of paths from the annotation
    output
        paths - a list of isoforms for splice modules
    """
    logging.debug('Start generating isoforms . . .')
    tx_paths, G = s_graph.annotation, s_graph.get_graph()
    #for component in get_biconnected(G):
    # define subgraph of biconnected nodes
    #component_subgraph = nx.subgraph(G, component)
    component = sorted(component_subgraph.nodes(), key=lambda x: (x[0], x[1]))  # make sure it is sorted

    # trim tx_paths to only contain paths within component_subgraph
    tmp = set()
    for p in tx_paths:
        # make sure this tx path has the biconnected component
        if component[0] in p and component[-1] in p:
            tmp.add(tuple(
                p[p.index(component[0]):p.index(component[-1]) + 1]))  # make sure there is no redundant paths
    tx_paths = sorted(list(tmp), key=lambda x: (x[0], x[1]))

    # assert statements about the connectivity of the graph
    assert nx.is_weakly_connected(component_subgraph), 'Yikes! expected weakly connected graph'
    assert nx.is_biconnected(component_subgraph.to_undirected()), 'Yikes! expected a biconnected component'

    # assert statements about AFE/ALE testing
    num_first_exons = len(filter(lambda x: len(component_subgraph.predecessors(x)) == 0, component_subgraph.nodes()))
    assert num_first_exons <= 1, 'Yikes! AFE like event is not expected'
    num_last_exons = len(filter(lambda x: len(component_subgraph.successors(x)) == 0, component_subgraph.nodes()))
    assert num_last_exons <= 1, 'Yikes! ALE like event is not expected'

    # mark all edges as unvisited
    visited_edges = {}
    for u, v in component_subgraph.edges():
        try:
            visited_edges[u][v] = 0
        except KeyError:
            visited_edges[u] = {}
            visited_edges[u][v] = 0
    num_visited_edges = 0  # no edges visited to start with

    # mark visited edges if in known annotation
    for path in tx_paths:
        for i in range(len(path) - 1):
            try:
                if not visited_edges[path[i]][path[i + 1]]:
                    visited_edges[path[i]][path[i + 1]] = 1
                    num_visited_edges += 1
            except:
                pass

    # Iterate until atleast one isoform explains each edge
    component_paths = tx_paths
    while num_visited_edges < component_subgraph.number_of_edges():
        # find maximum weight path
        path, no_new_edges = bellman_ford_longest_path(component_subgraph,
                                                       num_nodes=G.number_of_nodes(),
                                                       visited=visited_edges)
        if not no_new_edges:
            component_paths.append(path)

        # bookeeping on visited edges
        for i in range(len(path) - 1):
            if visited_edges[path[i]][path[i + 1]] == 0:
                visited_edges[path[i]][path[i + 1]] = 1
                num_visited_edges += 1
            G[path[i]][path[i + 1]]['weight'] = G[path[
                i]][path[i + 1]]['weight'] / 10.0  # down scale path

    logging.debug('Finished generating isoforms.')
    tmp = sorted(component)  # make sure nodes are sorted
    # component_paths += [filter(lambda x: x >= start and x <= stop, p) for p in tx_paths]  # add paths known from transcript annotation
    component_paths = map(list, list(
        set(map(tuple, component_paths))))  # remove redundant paths

    logging.debug('Start read count EM algorithm . . . ')
    count_info = read_count_em(component_paths, component_subgraph)
    logging.debug('Finished calculating counts.')
    #original_count = sum([component_subgraph[arc_tail][arc_head]['weight'] for arc_tail, arc_head in component_subgraph.edges()])  # total actual read counts that the subgraph has
    #module_info.append([[start, stop], list(count_info), original_count])   # module starts at first node and ends at last node
    # paths.append(component_paths)  # add paths

    return component_paths, count_info
