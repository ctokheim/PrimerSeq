import splice_graph as sg
import utils


def get_isforms_and_counts(self, line, options):
    # get information about each row
    ID, target_coordinate = line[:2]
    strand = target_coordinate[0]
    chr = utils.get_chr(target_coordinate[1:])
    tmp_start, tmp_end = utils.get_pos(target_coordinate)

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

    paths_list = []
    counts_list = []
    for my_splice_graph in bam_splice_graphs:
        # this case is meant for user-defined flanking exons
        if line[10] == '-1' and line[12] == '-1':

        # this case is meant for automatic choice of flanking exons
        else:
            # not the best use of the ExonSeek object, initially intended to find appropriate flanking exons
            # but in this case ExonSeek is used to get the transcripts and associate counts
            exon_seek_obj = ExonSeek(utils.get_pos(target_coordinate), my_splice_graph, ID, options['psi'], None, None)
            all_paths, upstream, downstream, component, psi_target, psi_upstream, psi_downstream = exon_seek_obj.get_info()
            paths_list.append(exon_seek_obj.paths)
            counts_list.append(exon_seek_obj.counts)
    return paths_list, counts_list, gene_name  # return the tx paths and count information for a single AS event


def user_defined_exons(tmp_sg, line):
    chr, strand = utils.get_chr(line[utils.TARGET]), line[utils.TARGET][0]  # get chr and strand
    upstream_exon = utils.get_pos(line[utils.UPSTREAM_EXON])  # get user-defined flanking exons
    downstream_exon = utils.get_pos(line[utils.DOWNSTREAM_EXON])
    first_primer, second_primer = utils.get_primer_coordinates(utils.PRIMER_COORD)

    # get possible exons for primer amplification
    tmp = sorted(tmp_sg.get_graph().nodes(), key=lambda x: (x[0], x[1]))
    if tmp_sg.strand == '+':
        my_exons = tmp[tmp.index(upstream_exon):tmp.index(downstream_exon) + 1]
    else:
        my_exons = tmp[tmp.index(downstream_exon):tmp.index(upstream_exon) + 1]

    # Use correct tx's and estimate counts/psi
    all_paths = algs.AllPaths(tmp_sg,
                              my_exons,
                              utils.get_pos(target_coordinate),  # tuple (start, end)
                              chr=chr,
                              strand=strand)
    all_paths.trim_tx_paths()
    all_paths.set_all_path_coordinates()
    paths, counts = all_paths.estimate_counts()  # run EM algorithm
    return paths, counts

def primerseq_defined_exons():
    pass


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
