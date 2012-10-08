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
        gene_dict[gene_id].setdefault('start', 0)
        gene_dict[gene_id]['start'] = min(gene_dict[gene_id]['start'], tx_path[0][0])  # change start if this tx has lower start position
        gene_dict[gene_id].setdefault('end', float('inf'))
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


def annotation_graph(g, FILTER_FACTOR=1000):
    """
    Create a nx DiGraph from list of tx in gene. FILTER_FACTOR defines a cutoff for
    using a tx of a gene. A tx must have x num of exon where x > MAX tx exon num / FILTER_FACTOR.
    Note by default FILTER_FACTOR is set to 1000 to prevent filtering of low exon number tx.
    """
    # filter low exon tx
    max_exons = max(map(len, g))  # figure out max num exons
    txs = filter(lambda x: len(x) > max_exons / FILTER_FACTOR, g)

    # create graph
    graph = nx.DiGraph()
    for tx in txs:
        graph.add_path(tx)
    return graph


def no_edges_graph(exons):
    """
    Simple function that makes a DAG (nx.DiGraph) with only edges and no nodes.
    Meant to be used to before add_all_possible_edge_weights.
    """
    G = nx.DiGraph()
    G.add_nodes_from(exons)
    return G


def add_annotation_edge_weights(graph, exon_forms, weights):
    """
    Only try to find weights for already existing edges in the graph.
    This function is intended to add weight values to edges defined
    in the gtf annotation.
    """
    # add edge weights to edges from annotation
    for u, v in graph.edges():
        try:
            tmpChr = get_chr(exon_forms[u])
            start = get_start_pos(exon_forms[u])
            end = get_end_pos(exon_forms[v])
            tmpWeight = weights[tmpChr][start][end]
            graph[u][v]['weight'] = tmpWeight
        except KeyError:
            graph[u][v]['weight'] = 0.01
    return graph


def add_all_possible_edge_weights(graph, chr, weights, READ_THRESHOLD=1):  # use to have exon_forms rather than chr
    """
    Add edge/weights to graph if supported by atleast READ_THRESHOLD number of reads
    """
    # add novel edges if well supported
    sorted_nodes = sorted(graph.nodes())
    for i in range(len(sorted_nodes) - 1):
        for j in range(i + 1, len(sorted_nodes)):
            try:
                # tmpChr = get_chr(exon_forms[sorted_nodes[i]])
                start = sorted_nodes[i][1]   # get_start_pos(exon_forms[sorted_nodes[i]])
                end = sorted_nodes[j][0]    # get_end_pos(exon_forms[sorted_nodes[j]])
                if weights[chr][start][end] >= READ_THRESHOLD:
                    graph.add_edge(sorted_nodes[i], sorted_nodes[j])
                    graph[sorted_nodes[i]][sorted_nodes[
                        j]]['weight'] = weights[chr][start][end]
            except KeyError:
                pass

    return graph


def get_flanking_biconnected(name, target, graph, chr, strand, genome):
    # search through each biconnected component
    for component in algs.get_biconnected(graph):
        component = sorted(component, key=lambda x: (x[0], x[1]))  # ensure first component is first exon, etc
        if target in component[1:-1]:
            # define upstream/downstream flanking exon
            if strand == '+':
                upstream = component[0]
                downstream = component[-1]
            else:
                upstream = component[-1]
                downstream = component[0]

            # get possible lengths
            inc_length, skip_length = algs.all_path_lengths(graph,
                                                            component, target)

            # get sequence of upstream/target/downstream combo
            genome_chr = genome[chr]  # chr object from pygr
            upstream_seq, target_seq, downstream_seq = genome_chr[upstream[0]:upstream[1]], genome_chr[target[0]:target[1]], genome_chr[downstream[0]:downstream[1]]
            if strand == '-':
                upstream_seq, target_seq, downstream_seq =  \
                    -upstream_seq, -target_seq, -downstream_seq

            return [strand, name[1:],
                    chr + ':' + '-'.join(map(str, upstream)),
                    chr + ':' + '-'.join(map(str, downstream)),
                    inc_length, skip_length, str(upstream_seq).upper(),
                    str(target_seq).upper(), str(downstream_seq).upper()]
    return name + ' was not found in a biconnected component'


def main(options, args_output='tmp/debug.json'):
    """
    The gtf main function is the function designed to be called from other scripts. It iterates through each target
    exons and returns the necessary information for primer design.
    """
    fasta, args_gtf, args_target = options['fasta'], options['gtf'], options['target']
    genome = SequenceFileDB(fasta)  # load fasta file from file
    gene_dict, gene_lookup = gene_annotation_reader(args_gtf)  # create splice graphs based on gene annotation

    # the sam object interfaces with the user specified BAM/SAM file!!!
    sam_obj = sam.Sam(options['rnaseq'])

    # iterate through each target exon
    output = []  # output from program
    for line in args_target:  # was line in handle
        line = line.strip()
        strand = line[0]
        try:
            genes = list(set(gene_lookup[line]))
        except KeyError:
            output.append('The coordinates of ' + line + ' did not match an exon')
            continue  # skip to next if exon lookup failed
        target = get_pos(line)
        chr = gene_dict[genes[0]]['chr']

        # check how many genes the exon matches
        if len(genes) > 1:
            output.append(str(line) + 'is in multiple genes ' + str(genes))
        elif len(genes) < 1:
            output.append(str(line) + 'is not in any genes')
        # construct output if no problems with genes
        elif options['annotation_flag']:
            # get information on flanking constitutive exons
            tmp_graph = annotation_graph(gene_dict[genes[0]]['graph'])  # construct graph based on annotation
            tmp = get_flanking_biconnected(line, target,
                                           tmp_graph,
                                           chr,
                                           strand,
                                           genome)
            output.append(tmp)
        elif options['rnaseq_flag']:
            tmp_graph = no_edges_graph(list(gene_dict[genes[0]]['exons']))
            tmp_start, tmp_end = target
            edge_weights = sam_obj.extractSamRegion(chr, tmp_start, tmp_end)
            tmp_graph = add_all_possible_edge_weights(tmp_graph.copy(),
                                                      chr,
                                                      edge_weights)  # add any edge w/ weight if sufficient read support
            tmp = get_flanking_biconnected(line, target,
                                           tmp_graph,
                                           chr,
                                           strand,
                                           genome)
            output.append(tmp)

    # if they specify a output file then write to it
    if args_output:
        with open(args_output, 'wb') as write_handle:
            json.dump(output, write_handle, indent=4)

    return output


if __name__ == '__main__':
    """Running this script directly is only for debug purposes"""
    # process command line arguments
    parser = argparse.ArgumentParser(description='Get flanking constitutive exons')
    parser.add_argument('-g', '--gtf', action='store', dest='gtf', required=True,
                        help='annotation file with legitimate gene_id\'s')
    parser.add_argument('-t', '--target', action='store', dest='target', required=True,
                        help='file of list of coordinate targets')
    parser.add_argument('-f', '--fasta', action='store', dest='fasta', required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--annotaton', dest='annotation_flag', action='store_true')
    group.add_argument('--rnaseq', dest='rnaseq_flag', action='store_true')
    parser.add_argument('-o', '--output', action='store', dest='output', required=True)
    options = vars(parser.parse_args())

    options['target'] = options['target'].replace('dash', '-').split(',')  # fix bug with - as input for strand

    # call main function
    main(options, options['output'])
