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

import csv
import re
import argparse
from itertools import tee, imap


class Gtf(object):
    """
    Separates out gtf parsing from iterating over records. Perhaps a flimsy
    use of a class.
    """
    def __init__(self, gtf_line):
        self.gtf_list = gtf_line
        self.seqname, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attribute = gtf_line  # These indexes are defined by the GTF spec
        self.attribute = dict(
            map(lambda x: re.split('\s+', x.replace('"', '')),
                re.split('\s*;\s*', self.attribute.strip().strip(';'))))  # convert attrs to dict

        self.start, self.end = int(self.start) - 1, int(self.end)

    def __str__(self):
        return '\t'.join(self.gtf_list) + '\n'


def gtf_reader(fileObject, delim):
    """
    Iterate over a file to extract 'exon' features of tx
    """
    for gtf_line in csv.reader(fileObject, delimiter=delim):
        gtf = Gtf(gtf_line)

        # only use exon features
        if gtf.feature.lower() == 'exon':
            yield Gtf(gtf_line)


def sort_gtf(file_name, output):
    gtf_list = []
    with open(file_name) as handle:
        for gtf_line in csv.reader(handle, delimiter='\t'):
            gtf = Gtf(gtf_line)
            if gtf.feature.lower() == 'exon':
                gtf_list.append(gtf)
    gtf_list.sort(key=lambda x: (x.seqname, x.attribute['gene_id'], x.attribute['transcript_id'], x.start, x.end))

    # write the contents back to a file
    with open(output, 'wb') as write_handle:
        for gtf_obj in gtf_list:
            write_handle.write(str(gtf_obj))


def is_sorted(iterable, compare):
    """Returns if iterable is sorted given the definition of
    a compare function"""
    a, b = tee(iterable)
    next(b, None)
    return all(imap(compare, a, b))


def gtf_compare(a, b):
    """compare two lines of a GTF file to see if they are
    correctly sorted as "a" before "b"."""
    tmp_list = [a, b]  # proposed sorted order
    srt_list = sorted(tmp_list, key=lambda x: (x.seqname,
                                               x.attribute['gene_id'],
                                               x.attribute['transcript_id'],
                                               x.start, x.end))  # actual sorted order

    # return flag for whether the two entries were in sorted order
    flag = (tmp_list == srt_list)
    return flag


def gtf_iter_reader(handle):
    """Iterator yielding GTF objects."""
    for gtf_line in csv.reader(handle, delimiter='\t'):
        yield Gtf(gtf_line)


def is_gtf_sorted(file_name):
    """Returns Boolean for if gtf is sorted."""
    with open(file_name) as handle:
        mygtf_reader = gtf_iter_reader(handle)
        return is_sorted(mygtf_reader, gtf_compare)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Either performs proper sorting of GTF for PrimerSeq or checks if GTF is sorted.
                                     For Sorting GTF:\npython gtf.py -i unsorted.gtf -o sorted.gtf,
                                     For checking if GTF is sorted:\npython gtf.py -c not_sure_if_sorted.gtf""")
    parser.add_argument('-i',
                        type=str,
                        default='',
                        dest='gtf',
                        action='store',
                        help='path to gtf file to sort')
    parser.add_argument('-o',
                        type=str,
                        default='',
                        dest='output',
                        action='store',
                        help='path name of properly sorted gtf')
    parser.add_argument('-c',
                        type=str,
                        default='',
                        dest='is_sorted',
                        action='store',
                        help='path to gtf file to check if sorted correctly')
    args = parser.parse_args()

    if args.gtf and args.output:
        sort_gtf(args.gtf, args.output)  # do the work of sorting
    elif args.is_sorted:
        if is_gtf_sorted(args.is_sorted):
            print '%s is correctly sorted' % (args.is_sorted)
        else:
            print '%s is not correctly sorted. please sort before use.' % (args.is_sorted)
    else:
        print 'You must enter either both the -i and -o options or just the -c option.'
