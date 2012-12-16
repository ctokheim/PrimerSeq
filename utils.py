def get_start_pos(coordinate):
    """
    get start from 'chr:start-stop'
    """
    return int(coordinate.split(":")[1].split("-")[0])


def get_end_pos(coordinate):
    """
    get stop from 'chr:start-stop'
    """
    return int(coordinate.split(":")[1].split("-")[1])


def get_pos(coordinate):
    """
    get (start, stop) from 'chr:start-stop'
    """
    return tuple(map(int, coordinate.split(":")[1].split("-")))


def get_chr(coordinate):
    """
    get chr from 'chr:start-stop'
    """
    return coordinate.split(":")[0]


def construct_coordinate(chr, start, end):
    return '%s:%s-%s' % (chr, str(start), str(end))


def merge_list_of_dicts(list_of_dicts):
    '''
    This function mereges multiple dicts contained read counts from SAM/BAM file
    into one dictionary.
    '''
    merged_dict = {}
    for tmp_dict in list_of_dicts:
        all_keys = set(merged_dict) | set(tmp_dict)
        for key in all_keys:
            merged_dict[key] = merged_dict.get(key, 0) + tmp_dict.get(key, 0)
    return merged_dict


def calc_product_length(path, primer_coord):
    """
    Calculate product length based on the primer coordinates
    """
    # calculate length between primers
    tmp_len = 0
    first_exon_primer, second_exon_primer = False, False  # flag for telling if a primer was contained within an exon
    flag = False
    for start, end in path:
        if start <= primer_coord[0][0] and end >= primer_coord[0][1]:
            tmp_len += end - primer_coord[0][1]
            flag = True
            first_exon_primer = True
        elif start <= primer_coord[1][0] and end >= primer_coord[1][1]:
            tmp_len += primer_coord[1][0] - start
            flag = False
            second_exon_primer = True
        elif flag:
            tmp_len += end - start

    if not (first_exon_primer and second_exon_primer):
        # case where either one or both primers were not located within an exon
        final_len = 0
    else:
        # add length between primers to actual length of sequnces for each primer
        first_primer_len = primer_coord[0][1] - primer_coord[0][0]
        second_primer_len = primer_coord[1][1] - primer_coord[1][0]
        final_len = tmp_len + first_primer_len + second_primer_len
    return final_len


class PrimerSeqError(Exception):
    """
    Used to uniquely catch this exception so that only one primer
    fails to be designed
    """
    pass


class InSilicoPcrUrl(object):
    '''
    Construct url for In-Silico PCR. This class constructs all of
    the GET submission parameters.
    '''
    def __init__(self, genome='', assembly='', target='',
                 forward='', reverse='', max_size=4000, perfect=15,
                 good=15, flip=0):
        # initialize variables
        self.genome = genome
        self.assembly = assembly
        self.forward = forward
        self.reverse = reverse
        self.max_size = max_size
        self.perfect = perfect
        self.good = good
        self.flip = flip
        if target == 'Genome':
            self.target = 'genome'
        elif target == 'UCSC Genes':
            self.target = self.assembly + 'Kg'

        # base url for ucsc in-silico pcr
        self.base_url = 'http://genome.ucsc.edu/cgi-bin/hgPcr'

    def get_url(self):
        get_params = '?org=%s&db=%s&wp_target=%s&wp_f=%s&wp_r=%s&wp_size=%s&wp_perfect=%s&wp_good=%s&boolshad.wp_flipReverse=%s&Submit=submit' % (
                            self.genome, self.assembly, self.target, self.forward, self.reverse, str(self.max_size), str(self.perfect), str(self.good), str(self.flip))
        return self.base_url + get_params

