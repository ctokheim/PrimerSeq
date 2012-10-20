import csv
import re


class Gtf(object):
    """
    Separates out gtf parsing from iterating over records. Perhaps a flimsy
    use of a class.
    """
    def __init__(self, gtf_line):
        self.seqname, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attribute = gtf_line  # These indexes are defined by the GTF spec
        self.attribute = dict(
            map(lambda x: re.split('\s+', x.replace('"', '')),
                re.split('\s*;\s*', self.attribute.strip().strip(';'))))  # convert attrs to dict

        self.start, self.end = int(self.start) - 1, int(self.end)

def gtf_reader(fileObject, delim):
    """
    Iterate over a file to extract 'exon' features of tx
    """
    for gtf_line in csv.reader(fileObject, delimiter=delim):
        gtf = Gtf(gtf_line)

        # only use exon features
        if gtf.feature.lower() == 'exon':
            yield Gtf(gtf_line)

