import csv
import argparse
import sys
import re


def main(options):
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
        if not reduce(lambda x, y: x and x >= options['anchor'], incs): continue  # small anchor length

        # add 1 to weights
        start_pos = int(line[POS]) - 1  # init to start pos of mapped read
        for i in range(len(incs) - 1):
            jct_start = start_pos + incs[i]
            jct_stop = jct_start + skips[i]
            weights.setdefault((line[RNAME], jct_start, jct_stop), 1)
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
