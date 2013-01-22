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

"""
author: Collin Tokheim

This script changes the tx names in knownGene annotation to the UCSC gene name displayed in the genome browser.

It needs a gtf/bed file and the kgXref text file dump from UCSC. To obtain the kgXref text file do the following:
1. Go to the [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables)
2. Select Genes and Gene Prediction tracks from the group dropdown
3. Select UCSC Genes from the track dropdown
4. Select kgXref from the table dropdown
5. Make sure the output format is 'all fields from selected table'
6. Click 'get output'
"""

import argparse
import csv
import re
import operator
import sys


class Attribute(object):
    """
    Class to handle the attribute column in the gtf
    """
    def __init__(self, gtf_line):
        ATTRIBUTE = 8  # index num of attribute in gtf file
        self.attribute = dict(
            map(lambda x: re.split('\s+', x.replace('"', ''), maxsplit=1),
                re.split('\s*;\s*', gtf_line[ATTRIBUTE].strip().strip(';'))))  # convert attrs to dict

    def __str__(self):
        return ' '.join(['%s "%s";' % (key, self.attribute[key].replace(' ', '_')) for key in self.attribute])


def modify_gtf(options, txToGene):
    # read / modify gtf
    with open(options['annotation']) as handle:
        output, failed_lookup = [], 0
        for line in csv.reader(handle, delimiter='\t'):
            attr = Attribute(line)
            try:
                attr.attribute['gene_id'] = txToGene[attr.attribute['transcript_id']]  # change gene_id to kgXref gene in Attribute obj
            except KeyError:
                failed_lookup += 1
            line[8] = str(attr)  # perform switch of gene_id
            output.append(line)  # append line with changed gene_id to the output list
    return output, failed_lookup


def modify_bed(options, txToGene):
    # read / modify bed file
    with open(options['annotation']) as handle:
        output, failed_lookup = [], 0
        for line in csv.reader(handle, delimiter='\t'):
            line[1] = int(line[1])
            try:
                line[3] = txToGene[line[3]]
            except:
                failed_lookup += 1
            output.append(line)
    return output, failed_lookup


def main(options):
    # parse kgXref text output
    with open(options['kgxref']) as handle:
        txToGene = dict(map(lambda x: [x[0], x[4].replace(' ', '_')],
                            filter(lambda x: not x[0].startswith('#'),
                            csv.reader(handle, delimiter='\t'))))

    if options['annotation'].endswith('.gtf'):
        output, failed_lookup = modify_gtf(options, txToGene)
    elif options['annotation'].endswith('.bed'):
        output, failed_lookup = modify_bed(options, txToGene)
    else:
        print 'Please input a bed or gtf file'
        sys.exit(1)
    output.sort(key=operator.itemgetter(0, 1))

    # write output gtf
    with open(options['output'], 'wb') as handle:
        for line in output:
            line[1] = str(line[1])
            handle.write('\t'.join(line) + '\n')

    print failed_lookup, 'transcript ids were not found in the kgXref text file'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='change tx names to UCSC gene name displayed in genome browser')
    parser.add_argument('-k', dest='kgxref', required=True, help='textfile of kgXref table in UCSC genome browser')
    parser.add_argument('-a', dest='annotation', required=True, help='Path to either a bed file or gtf file')
    parser.add_argument('-o', dest='output', required=True, help='output file path to gtf/bed with correct gene names')
    options = vars(parser.parse_args())

    main(options)
