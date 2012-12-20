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

'''
File: jct_counts.py
Author: Collin Tokheim
Description: This file reads a sam file and outputs junction read
counts to file. Since this is a python script, it is only advisable
to use this on a sam file that only contains reads in your region of
interest.
'''
import csv
import argparse
import sys
import re


def main(options):
    '''
    Ouptuts jct read counts from a SAM file into the specified output file.
    '''
    # define input as either a file or stdin
    if options['sam'] == 'stdin':
        file_input = sys.stdin
    else:
        file_input = open(options['sam'])

    # iterate through each read
    inc_search, skip_search, error_search = re.compile('\d+(?=M)'), re.compile('\d+(?=N)'), re.compile('[^MN0-9]')
    QNAME, FLAG, RNAME, POS, MAPQ, CIGAR = range(6)  # define some SAM columns
    weights = {}  # hold edge count information
    for line in csv.reader(file_input, delimiter='\t'):
        if line[0][0] == '@': continue  # skip line if head character
        if error_search.search(line[CIGAR]): continue  # other cigar characters found

        # get M/N cigar info
        skips = map(int, skip_search.findall(line[CIGAR]))
        if not len(skips): continue  # skip if not junction
        incs = map(int, inc_search.findall(line[CIGAR]))
        valid_anchor_lengths = map(lambda x: x >= options['anchor'], incs)
        # if not reduce(lambda x, y: x and x >= options['anchor'], incs): continue  # small anchor length

        # add 1 to weights
        start_pos = int(line[POS]) - 1  # init to start pos of mapped read
        for i in range(len(incs) - 1):
            jct_start = start_pos + incs[i]
            jct_stop = jct_start + skips[i]
            if valid_anchor_lengths[i] and valid_anchor_lengths[i+1]:
                weights.setdefault((line[RNAME], jct_start, jct_stop), 0)
                weights[(line[RNAME], jct_start, jct_stop)] += 1
            start_pos = jct_stop
    file_input.close()  # close input

    # convert dict to list so it can written in tabular form
    output = [[chr, start, stop, weights[(chr, start, stop)]] for chr, start, stop in weights]
    output.sort(key=lambda x: (x[0], x[1], x[2]))

    # define output as either a file or stdout
    if options['output'] == 'stdout':
        file_output = sys.stdout
    else:
        file_output = open(options['output'], 'wb')

    csv.writer(file_output, delimiter='\t').writerows(output)
    file_output.close()  # close output


if __name__ == '__main__':
    # List optional command line argument
    parser = argparse.ArgumentParser(description='Outputs a junction/edge count file from RNA-Seq reads. Assumes you have already quality filtered.')
    parser.add_argument('-s', '--sam', default='stdin', action='store', dest='sam', help='sam file name or `stdin` (default)')
    parser.add_argument('-o', '--output', default='stdout', action='store', dest='output', help='output file name or `stdout` (default)')
    parser.add_argument('-a', '--anchor', default=1, type=int, action='store', dest='anchor', help='only count junction reads with an anchor length GTEQ this value (Default=1)')
    options = vars(parser.parse_args())

    main(options)  # run script
