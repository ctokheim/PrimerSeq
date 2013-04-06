#!/usr/bin/env python
# Copyright (C) 2013  ctokheim
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

import numpy as np
import logging


def debug_graph(G):
    """Simply prints/logs info about a networkx graph"""
    logging.debug("*" * 20)
    logging.debug(G)
    logging.debug(G[u][v])
    utils.print_graph(G)
    logging.debug("*" * 20)


def debug_edges(paths, edges):
    path_edge_set = set()
    for path in paths:
        for i in xrange(len(path) - 1):
            path_edge_set.add((path[i], path[i+1]))
    print path_edge_set
    print edges


def construct_read_count_vector(graph, index_to_edge):
    """
    Construct the vector consisting of read counts
    """
    num_jcts = graph.number_of_edges()
    counts_vector = np.zeros(num_jcts)

    # handle edge weights representing exon-exon junctions
    for i in xrange(num_jcts):
        u, v = index_to_edge[i]
        try:
            counts_vector[i] = graph[u][v]['weight']
        except KeyError:
            debug_graph(graph)
            raise
    return counts_vector


def construct_uncommited_matrix(Y, counts_vec, bcc_paths, edge_to_index):
    """
    Construct a two-dimensional array Y for incorporating read count
    information in regards to edge and transcript concordance
    """
    # set the uncommited matrix to have jct cts
    for tx_index, path in enumerate(bcc_paths):
        for i in range(len(path) - 1):
            edge = (path[i], path[i + 1])
            read_count_index = edge_to_index[edge]
            Y[read_count_index][tx_index] = counts_vec[read_count_index]
    return Y


def multinomial_em(bcc_paths, sub_graph):
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
    read_counts = construct_read_count_vector(sub_graph, indexToEdge)
    total_counts = np.sum(read_counts)

    # set up the uncommited matrix Y
    num_edges = sub_graph.number_of_edges()
    Y = np.zeros((num_edges, num_tx))
    Y = construct_uncommited_matrix(Y, read_counts, bcc_paths, edgeToIndex)

    # set up p the probability array
    p = np.ones(num_tx) * 1. / num_tx

    # start EM
    counter, MAX_ITERS = 0, 10000  # these variables impose a max iteration restriction of 10000 on the EM algorithm
    epsilon = float('inf')
    THRESHOLD = .0001
    while epsilon > THRESHOLD and counter < MAX_ITERS:
        # E-step
        for i, row in enumerate(Y):
            try:
                low_count_indices = row < 1e-5
                row[low_count_indices] = 0
                # if np.sum(row) != 0 and row.dot(p) != 0:
                if row.dot(p) != 0:
                    Y[i] = row * p / row.dot(p) * read_counts[i]
            except:
                logging.debug('Row: ' + str(row))
                logging.debug('P: ' + str(p))
                logging.debug('Read count: ' + str(read_counts[i]))
                raise utils.PrimerSeqError('Numerical error in EM algorithm')

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
    if not inc_count and not skip_count:
        psi = -1  # -1 indicates a divide by zero error
    else:
        psi = float(inc_count) / (inc_count + skip_count)
    return psi
