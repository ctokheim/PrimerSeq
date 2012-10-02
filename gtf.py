import csv
import re
import networkx as nx
import itertools as it
import argparse
from pygr import worldbase
import json
import algorithms as algs


def gtf_reader(fileObject, delim):
    """
    A python generator to separate out gtf parsing. Main purpose to is to parse the attribute columnn
    of the gtf and to convert numerical columns to integers.
    """
    ATTRIBUTE = 8  # index number of attribute column
    for gtf_line in csv.reader(fileObject, delimiter=delim):
        gtf_line[ATTRIBUTE] = dict(
            map(lambda x: re.split('\s+', x.replace('"', '')),
                re.split('\s*;\s*', gtf_line[ATTRIBUTE].strip().strip(';'))))  # convert attrs to dict

        # try to convert to int if possible
        for i in range(len(gtf_line)):
            try:
                gtf_line[i] = int(gtf_line[i])
            except (ValueError, TypeError):
                pass
        gtf_line[3] -= 1  # make zero based
        yield gtf_line


def splice_graph_reader(file_path, FILTER_FACTOR=2):
    """
    *creates two data structures from gtf:*

        gene_dict is a dictionary with gene_id's as keys:
            * gene_dict['My_favorite_gene']['graph'] = networkx DAG (splice graph)
            * gene_dict['My_favorite_gene']['chr'] = gene chromosome (chrXX)
            * gene_dict['My_favorite_gene']['strand'] = gene strand (+, -)

        gene_lookup is a dictionary to lookup the gene that an exon coordinate belongs to:
            * gene_lookup['(strand)(chr):(start)-(end)'] = ['My_favorite_gene', ...]
    """
    # iterate through each gtf feature
    file_input = open(file_path)
    SEQNAME, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME, ATTRIBUTE = range(9)  # These indexes are defined by the GTF spec
    gene_dict, gene_lookup = {}, {}
    for tx_id, tx in it.groupby(gtf_reader(file_input, delim='\t'), lambda x: x[ATTRIBUTE]['transcript_id']):
        tx = filter(lambda x: x[FEATURE] == 'exon', tx)  # only use exon features in gtf
        gene_id = tx[0][ATTRIBUTE]['gene_id']
        strand = tx[0][STRAND]

        # sort exons
        tx_path = sorted([(exon[START], exon[END]) for exon in tx],
                         key=lambda x: (x[0], x[1]))  # needs to be sorted because gtf files might not have them in proper order

        # add info to gene_dict
        gene_dict.setdefault(gene_id, {})  # add the gene key if it doesn't exist
        gene_dict[gene_id].setdefault('chr', tx[0][SEQNAME])  # add chr if doesn't exist
        gene_dict[gene_id].setdefault('strand', strand)  # add strand if doesn't exist
        gene_dict[gene_id].setdefault('graph', []).append(tx_path)  # append the tx path

        # lookup dict for which exon is in which gene
        for exon in tx:
            pos_string = '%s%s:%d-%d' % (strand, exon[SEQNAME], exon[START], exon[END])
            gene_lookup.setdefault(pos_string, []).append(gene_id)

    # apply filtering before creating graph
    for g in gene_dict:
        # filter low exon tx
        max_exons = max(map(len, gene_dict[g]['graph']))  # figure out max num exons
        txs = filter(lambda x: len(x) > max_exons / FILTER_FACTOR, gene_dict[g]['graph'])

        # create graph
        graph = nx.DiGraph()
        for tx in txs:
            graph.add_path(tx)
        gene_dict[g]['graph'] = graph

    file_input.close()
    return gene_dict, gene_lookup


def main(args_target, args_gtf, args_output='debug.json'):
    """
    The gtf main function is the function designed to be called from other scripts. It iterates through each target
    exons and returns the necessary information for primer design.
    """
    # guess the genome based on the gtf name
    nameToPygr = {'hg19': worldbase.Bio.Seq.Genome.HUMAN.hg19(), 'GRCh37': worldbase.Bio.Seq.Genome.HUMAN.hg19(), 'mm9': worldbase.Bio.Seq.Genome.MOUSE.mm9(), 'NCBIM37': worldbase.Bio.Seq.Genome.MOUSE.mm9()}
    for key in nameToPygr:
        if key in args_gtf:
            genome = nameToPygr[key]
            break

    gene_dict, gene_lookup = splice_graph_reader(args_gtf)  # create splice graphs

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
        target = tuple(map(int, line.split(':')[1].split('-')))

        # check how many genes the exon matches
        if len(genes) > 1:
            output.append(str(line) + 'is in multiple genes ' + str(genes))
            continue  # skip if exon is in multiple genes
        elif len(genes) < 1:
            output.append(str(line) + 'is not in any genes')
            continue  # skip if not in genes

        # search through each biconnected component
        for component in algs.get_biconnected(gene_dict[genes[0]]['graph'].copy()):
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
                inc_length, skip_length = [], []
                sub_graph = nx.subgraph(gene_dict[genes[0]]['graph'], component)
                for path in nx.all_simple_paths(sub_graph.copy(), source=component[0], target=component[-1]):
                    if target in path:
                        inc_length.append(sum(map(lambda x: x[1] - x[0], path[1:-1])))  # length of everything but target exon and flanking constitutive exons
                    else:
                        skip_length.append(sum(map(lambda x: x[1] - x[0], path[1:-1])))  # length of everything but target exon and flanking constitutive exons

                # get sequence of upstream/target/downstream combo
                genome_chr = genome[gene_dict[genes[0]]['chr']]  # chr object from pygr
                upstream_seq, target_seq, downstream_seq = genome_chr[upstream[0]:upstream[1]], genome_chr[target[0]:target[1]], genome_chr[downstream[0]:downstream[1]]
                if strand == '-':
                    upstream_seq, target_seq, downstream_seq =  \
                        -upstream_seq, -target_seq, -downstream_seq

                output.append(
                    [strand, line[1:],
                        gene_dict[genes[0]]['chr'] + ':' + '-'.join(map(str, upstream)),
                        gene_dict[genes[0]]['chr'] + ':' + '-'.join(map(str, downstream)),
                        inc_length, skip_length, str(upstream_seq).upper(),
                        str(target_seq).upper(), str(downstream_seq).upper()])
                break  # a target exon can not be in multiple biconnected components
        else:
            output.append(
                line + ' was not found in a biconnected component')

    # if they specify a output file then write to it
    if args_output:
        with open(args_output, 'wb') as writeHandle:
            json.dump(output, writeHandle, indent=4)

    return output

if __name__ == '__main__':
    """Running this script directly is only for debug purposes""" 
    # process command line arguments
    parser = argparse.ArgumentParser(description='Get flanking constitutive exons')
    parser.add_argument('-g', '--gtf', action='store', dest='gtf', required=True,
                        help='annotation file with legitimate gene_id\'s')
    parser.add_argument('-t', '--target', action='store', dest='target', required=True,
                        help='file of list of coordinate targets')
    parser.add_argument('-o', '--output', action='store', dest='output', required=True)
    args = parser.parse_args()

    args.target = args.target.replace('dash', '-').split(',')  # fix bug with - as input for strand

    # call main function
    main(args.target, args.gtf, args.output)
