import splice_graph as sg
from exon_seek import ExonSeek
import algorithms as algs
import utils
import logging


def save_isforms_and_counts(line, options):
    # get information about each row
    ID, target_coordinate = line[:2]
    strand = target_coordinate[0]
    chr = utils.get_chr(target_coordinate[1:])
    tmp_start, tmp_end = utils.get_pos(target_coordinate)
    logging.debug('Saving isoform and count information for event %s . . .' % ID)

    # get information from GTF annotation
    gene_dict, gene_name = retrieve_gene_information(options,
                                                     strand, chr, tmp_start, tmp_end)

    # get edge weights
    edge_weights_list = [sam_obj.extractSamRegion(chr, gene_dict['start'], gene_dict['end'])
                         for sam_obj in options['rnaseq']]

    # construct splice graph for each BAM file
    bam_splice_graphs = sg.construct_splice_graph(edge_weights_list,
                                                  gene_dict,
                                                  chr,
                                                  strand,
                                                  options['read_threshold'],
                                                  options['min_jct_count'],
                                                  output_type='list',
                                                  both=options['both_flag'])

    for bam_ix, my_splice_graph in enumerate(bam_splice_graphs):
        # this case is meant for user-defined flanking exons
        if line[utils.PSI_UP] == '-1' and line[utils.PSI_DOWN] == '-1':
            # find path and count information
            paths, counts = user_defined_exons(my_splice_graph, line)

            # filter out single exon paths
            my_tmp = [(path, count) for path, count in zip(paths, counts) if len(path) > 1]
            paths, counts = zip(*my_tmp)
        # this case is meant for automatic choice of flanking exons
        else:
            paths, counts = primerseq_defined_exons(my_splice_graph, line, options['psi'])
        utils.save_path_info('%s.%d' % (ID, bam_ix),
                             paths, counts,
                             save_dir='tmp/indiv_isoforms/')
    logging.debug('Finished saving isoform and count information for event %s.' % ID)


def user_defined_exons(tmp_sg, line):
    chr, strand = utils.get_chr(line[utils.TARGET]), line[utils.TARGET][0]  # get chr and strand
    upstream_exon = utils.get_pos(line[utils.UPSTREAM_EXON])  # get user-defined flanking exons
    downstream_exon = utils.get_pos(line[utils.DOWNSTREAM_EXON])
    first_primer, second_primer = utils.get_primer_coordinates(line[utils.PRIMER_COORD])

    # get possible exons for primer amplification
    tmp = sorted(tmp_sg.get_graph().nodes(), key=lambda x: (x[0], x[1]))
    first_ex = utils.find_first_exon(first_primer, tmp)
    last_ex = utils.find_last_exon(second_primer, tmp)
    my_exons = tmp[first_ex:last_ex + 1]
    # if tmp_sg.strand == '+':
    #     my_exons = tmp[tmp.index(upstream_exon):tmp.index(downstream_exon) + 1]
    # else:
    #     my_exons = tmp[tmp.index(downstream_exon):tmp.index(upstream_exon) + 1]

    # Use correct tx's and estimate counts/psi
    all_paths = algs.AllPaths(tmp_sg,
                              my_exons,
                              utils.get_pos(line[utils.TARGET]),  # tuple (start, end)
                              chr=chr,
                              strand=strand)
    # all_paths.trim_tx_paths()
    all_paths.trim_tx_paths_using_primers(first_primer, second_primer)
    all_paths.set_all_path_coordinates()
    paths, counts = all_paths.estimate_counts()  # run EM algorithm
    return paths, counts


def primerseq_defined_exons(tmp_sg, line, psi_option):
    """
    Get information about counts and paths if using PrimerSeq to define the flanking exons.
    """
    # not the best use of the ExonSeek object, initially intended to find appropriate flanking exons
    # but in this case ExonSeek is used to get the transcripts and associate counts
    ID = line[utils.ID]
    tgt_pos = utils.get_pos(line[utils.TARGET])
    exon_seek_obj = ExonSeek(tgt_pos,
                             tmp_sg,
                             ID,
                             psi_option,
                             None,  # no defined upstream exon
                             None)  # no defined downstream exon
    all_paths, upstream, downstream, component, psi_target, psi_upstream, psi_downstream = exon_seek_obj.get_info()
    return exon_seek_obj.paths, exon_seek_obj.counts


def retrieve_gene_information(options, strand, chr, start, end):
    """
    Gets information from GTF annotation either by using gene name
    or by using weakly connected transcripts
    """
    # get information regarding the gene
    if options['no_gene_id']:
        # hopefully filter out junk, but only uses weakly connected
        gene_dict, gene_name = sg.get_weakly_connected_tx(options['gtf'],
                                                          strand, chr, start, end)
    else:
        # gets everything for a single gene
        gene_dict, gene_name = sg.get_from_gtf_using_gene_name(options['gtf'],
                                                               strand, chr, start, end)
    return gene_dict, gene_name
