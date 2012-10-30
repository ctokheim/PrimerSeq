'''
File: primer.py
Author: Collin Tokheim
Description: primer.py can design primers from the command line.
Note that PrimerApp.py (the GUI) just uses the primer module to design the
primers.
'''
import subprocess
import os
import shutil
import gtf
import splice_graph
import csv
import argparse  # command line parsing
import itertools as it
from pygr.seqdb import SequenceFileDB
import sam

# import for logging file
import logging
import datetime
import time
import sys
import traceback


def gene_annotation_reader(file_path, FILTER_FACTOR=2):
    """
    *creates two data structures from gtf:*

        gene_dict is a dictionary with gene_id's as keys:
            * gene_dict['chr']['My_favorite_gene']['graph'] = list of list of exons
            * gene_dict['chr']['My_favorite_gene']['strand'] = gene strand (+, -)
            * gene_dict['chr']['My_favorite_gene']['chr'] = chromosome
            * gene_dict['chr']['My_favorite_gene']['start'] = start of gene
            * gene_dict['chr']['My_favorite_gene']['end'] = end of gene
            * gene_dict['chr']['My_favorite_gene']['exons'] = the set of exons (nodes)
    """
    # logging.debug('Started reading %s' % file_path)
    # iterate through each gtf feature
    file_input = open(file_path)
    gene_dict = {}
    for tx_id, tx in it.groupby(gtf.gtf_reader(file_input, delim='\t'), lambda x: x.attribute['transcript_id']):
        tx = list(tx)
        if len(tx) == 0: continue  # no 'exon' feature case

        gene_id = tx[0].attribute['gene_id']
        strand = tx[0].strand

        # sort exons
        tx_path = sorted([(exon.start, exon.end) for exon in tx],
                         key=lambda x: (x[0], x[1]))  # needs to be sorted because gtf files might not have them in proper order

        # add info to gene_dict
        gene_dict.setdefault(tx[0].seqname, {})
        gene_dict[tx[0].seqname].setdefault(gene_id, {})  # add the gene key if it doesn't exist
        gene_dict[tx[0].seqname][gene_id].setdefault('chr', tx[0].seqname)  # add chr if doesn't exist
        gene_dict[tx[0].seqname][gene_id].setdefault('strand', strand)  # add strand if doesn't exist
        gene_dict[tx[0].seqname][gene_id].setdefault('graph', []).append(tx_path)  # append the tx path
        gene_dict[tx[0].seqname][gene_id].setdefault('start', float('inf'))
        gene_dict[tx[0].seqname][gene_id]['start'] = min(gene_dict[tx[0].seqname][gene_id]['start'], tx_path[0][0])  # change start if this tx has lower start position
        gene_dict[tx[0].seqname][gene_id].setdefault('end', 0)
        gene_dict[tx[0].seqname][gene_id]['end'] = max(gene_dict[tx[0].seqname][gene_id]['end'], tx_path[-1][1])
        gene_dict[tx[0].seqname][gene_id].setdefault('exons', set())
        for ex in tx_path:
            gene_dict[tx[0].seqname][gene_id]['exons'].add(ex)  # hold a set of non-redundant exons
    file_input.close()  # close gtf file
    # logging.debug('Finished reading %s' % file_path)
    return gene_dict


def call_primer3(target_string, jobs_ID):
    """
    Does the actual call to primer3. Will raise a CalledProcessError if primer3
    exits with a non-zero exit status.
    """
    logging.debug('Calling primer3 for %s . . .' % target_string)
    os.chdir('tmp')  # make sure primer3 files occur in tmp dir
    cmd = '../primer3/src/primer3_core < ' + jobs_ID + '.conf' + ' > ' + jobs_ID + '.Primer3'
    logging.debug(cmd)  # record command in log file
    subprocess.check_call(cmd, shell=True)
    logging.debug('Finished call of primer3 for %s' % target_string)
    os.chdir('..')  # go back to top level


def read_primer3(path):
    """
    Read the output from primer3_core into a dictionary containing
    the key:value relationship for the output tags.
    """
    primer3_dict = {}
    # read primer3 file
    with open(path, 'r') as infile:
        for line in infile:
            try:
                primer3_dict[line.split('=')[0]] = line.split('=')[1].rstrip()
            except IndexError:
                print line
    return primer3_dict


def mkdir_tmp():
    '''
    Make all the necessary directories for temporary files
    '''
    if not os.path.exists('tmp'): os.mkdir('tmp')
    if not os.path.exists('tmp/sam'): os.mkdir('tmp/sam')
    if not os.path.exists('tmp/jct'): os.mkdir('tmp/jct')
    if not os.path.exists('tmp/bed'): os.mkdir('tmp/bed')
    if not os.path.exists('tmp/wig'): os.mkdir('tmp/wig')
    if not os.path.exists('tmp/isoforms'): os.mkdir('tmp/isoforms')


def primer3(options, primer3_options):
    """
    The primer.py main function uses the gtf module to find information about constitutive flanking exons for the target exons of interest.
    It then designs primers by calling primer3. Next it parses the primer3 output and outputs the final results to a file. The output file
    is then emailed to the designated address in the command line parameters.
    """
    jobs_ID = options['job_id']

    # tmp directory
    outfiles_PATH = 'tmp/'  # directory for intermediate files
    mkdir_tmp()  # make any necessary tmp directories

    # read in targets
    # with open(options['target']) as handle:
    #    target_list = map(lambda x: x.strip().split('\t'), handle.readlines())
    #    options['target'] = target_list

    # find flanking exons
    logging.debug('Calling splice_graph.main to find flanking exons')
    flanking_info = splice_graph.main(options)
    logging.debug('Finished gtf.main')

    # iterate over all target sequences
    STRAND, EXON_TARGET, PSI_TARGET, UPSTREAM_TARGET, PSI_UPSTREAM, DOWNSTREAM_TARGET, PSI_DOWNSTREAM, ALL_PATHS, UPSTREAM_Seq, TARGET_SEQ, DOWNSTREAM_SEQ = range(11)
    output_list = []
    for z in range(len(flanking_info)):
        # no flanking exon information case
        #if type('') == type(flanking_info[z]):
        if len(flanking_info[z]) == 1:
            logging.debug(flanking_info[z][0])
            output_list.append(flanking_info[z])  # write problem msg
        # has flanking exon information case
        else:
            tar = options['target'][z][1] # flanking_info[z][1]  # target interval (used for print statements)
            tar_id = options['target'][z][0]
            ####################### Primer3 Parameter Configuration###########
            P3_FILE_FLAG = '1'
            PRIMER_EXPLAIN_FLAG = '1'
            PRIMER_THERMODYNAMIC_PARAMETERS_PATH = '../primer3/src/primer3_config/'
            SEQUENCE_ID = tar  # use the 'chr:start-stop' format for the sequence ID in primer3
            #SEQUENCE_TEMPLATE = flanking_info[z][UPSTREAM_Seq] + flanking_info[z][TARGET_SEQ].lower() + flanking_info[z][DOWNSTREAM_SEQ]
            #SEQUENCE_TARGET = str(len(flanking_info[z][UPSTREAM_Seq]) + 1) + ',' + str(len(flanking_info[z][TARGET_SEQ]))
            SEQUENCE_TEMPLATE = flanking_info[z][UPSTREAM_Seq] + flanking_info[z][TARGET_SEQ].lower() + flanking_info[z][DOWNSTREAM_SEQ]
            SEQUENCE_PRIMER_PAIR_OK_REGION_LIST = '0,' + str(len(flanking_info[z][UPSTREAM_Seq])) + ',' + str(len(flanking_info[z][UPSTREAM_Seq]) + len(flanking_info[z][TARGET_SEQ])) + ',' + str(len(flanking_info[z][DOWNSTREAM_SEQ]))
            #############################################################

            ####################### Write jobs_ID.conf##################
            with open(outfiles_PATH + jobs_ID + '.conf', 'w') as outfile:
                # hard coded options
                outfile.write('SEQUENCE_ID=' + SEQUENCE_ID + '\n')
                outfile.write('SEQUENCE_TEMPLATE=' + SEQUENCE_TEMPLATE + '\n')
                #outfile.write('SEQUENCE_TARGET=' + SEQUENCE_TARGET + '\n')
                outfile.write('SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=' + SEQUENCE_PRIMER_PAIR_OK_REGION_LIST + '\n')
                outfile.write('P3_FILE_FLAG=' + P3_FILE_FLAG + '\n')
                outfile.write('PRIMER_EXPLAIN_FLAG=' + PRIMER_EXPLAIN_FLAG + '\n')
                outfile.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=' + PRIMER_THERMODYNAMIC_PARAMETERS_PATH + '\n')  # make sure primer3 finds the config files

                # options from primer3.cfg
                for o in primer3_options:
                    outfile.write(o)
                outfile.write('=' + '\n')  # primer3 likes a '=' at the end of sequence params
                logging.debug('Wrote the input file (%s) for primer3' % (outfiles_PATH + jobs_ID + '.conf'))

            ###################### Primer3 #####################################
            if os.path.exists(outfiles_PATH + jobs_ID + '.Primer3'):
                os.remove(outfiles_PATH + jobs_ID + '.Primer3')  # delete old files

            call_primer3(tar, jobs_ID)  # command line call to Primer3!

            #################### Parse '.Primer3' ################################
            primer3_dict = read_primer3(outfiles_PATH + jobs_ID + '.Primer3')

            # checks if no output
            if(primer3_dict.keys().count('PRIMER_LEFT_0_SEQUENCE') == 0):
                logging.debug('No primer3 results for %s' % tar)
                output_list.append(['No Primer3 results for ' + tar])
                continue
            # there is output case
            else:
                logging.debug('There are primer3 results for %s' % SEQUENCE_ID)
                # get info about product sizes
                target_exon_len = len(flanking_info[z][TARGET_SEQ])
                Primer3_PRIMER_PRODUCT_SIZE = int(primer3_dict['PRIMER_PAIR_0_PRODUCT_SIZE']) - target_exon_len
                skipping_size_list = [str(int(j) + Primer3_PRIMER_PRODUCT_SIZE) for j in flanking_info[z][ALL_PATHS].skip_lengths]
                inclusion_size_list = [str(int(k) + Primer3_PRIMER_PRODUCT_SIZE) for k in flanking_info[z][ALL_PATHS].inc_lengths]
                skipping_size = ';'.join(skipping_size_list)
                inclusion_size = ';'.join(inclusion_size_list)

                # append results to output_list
                tmp = [tar_id, flanking_info[z][EXON_TARGET], flanking_info[z][PSI_TARGET], primer3_dict['PRIMER_LEFT_0_SEQUENCE'], primer3_dict['PRIMER_RIGHT_0_SEQUENCE'],
                       str((float(primer3_dict['PRIMER_LEFT_0_TM']) + float(primer3_dict['PRIMER_RIGHT_0_TM'])) / 2), skipping_size, inclusion_size,
                       flanking_info[z][UPSTREAM_TARGET], flanking_info[z][PSI_UPSTREAM], flanking_info[z][DOWNSTREAM_TARGET], flanking_info[z][PSI_DOWNSTREAM]]
                output_list.append(tmp)

    # write output information
    with open(outfiles_PATH + jobs_ID + '.txt', 'wb') as outputfile_tab:
        # define csv header
        header = ['ID', 'Target_Exon_coordinates', 'PSI_TARGET', 'Left_primer', 'Right_primer', 'Average_Tm',
                  'Skipping_production_size', 'Inclusion_production_size', 'Upstream_coordinates', 'PSI_UPSTREAM', 'Downstream_coordinates', 'PSI_DOWNSTREAM']
        output_list = [header] + output_list  # pre-pend header to output file
        csv.writer(outputfile_tab, delimiter='\t').writerows(output_list)  # output primer design to a tab delimited file
    shutil.copy(outfiles_PATH + jobs_ID + '.txt', options['output'])  # copy file to output destination


def main(options):
    """Reads in primer3 options, starts logging, and then calls the primer3 function to run primer3"""

    ##### Getting Start Time ######
    logging.debug('Start the program with [%s]', ' '.join(sys.argv))
    startTime = time.time()

    # read in primer3 options
    logging.debug('Reading in primer3.cfg . . .')
    primer3_options = []
    with open('primer3.cfg') as handle:
        for line in handle:
            # skip comment lines and lines with no set values
            if not line.startswith('#') and len(line.strip().split("=")[1]) > 0:
                primer3_options.append(line)
    logging.debug('Finished reading primer3.cfg.')

    # the primer3 function runs the primer3_core executable
    try:
        primer3(options, primer3_options)
    except:
        t, v, trace = sys.exc_info()
        logging.debug('ERROR! For more information read the following lines')
        logging.debug('Type: ' + str(t))
        logging.debug('Value: ' + str(v))
        logging.debug('Traceback:\n' + traceback.format_exc())
        raise

    # delete temporary sam files (may eventually delete more tmp files)
    logging.debug('Deleting tmp files')
    if not options['keep_temp']:
        tmp_sam = 'tmp/sam'
        for f in os.listdir(tmp_sam):
            os.remove(os.path.join(tmp_sam, f))

    ### record end of running primer3 ###
    logging.debug("Program ended")
    currentTime = time.time()
    runningTime = currentTime - startTime  # in seconds
    logging.debug("Program ran for %.2d:%.2d:%.2d" % (runningTime / 3600, (runningTime % 3600) / 60, runningTime % 60))


class ValidateRnaseq(argparse.Action):
    """
    Makes sure the -r parameter input is a SAM/BAM file
    """
    def __call__(self, parser, namespace, values, option_string=None):
        # if error print help and exit
        if not values.endswith('.sam') and not values.endswith('.bam'):
            parser.print_help()
            parser.exit(status=1,
                message='\n' + '#' * 40 + '\nThe -r parameter expects a SAM or BAM file' + '#' * 40 + '\n')
        # set value if no error,
        # simply assign the string as is
        setattr(namespace, self.dest, values)  # set the value

class ValidateCutoff(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # if error print help and exit
        if 0 <= values <= 1:
            setattr(namespace, self.dest, values)  # set the value
        else:
            parser.print_help()
            parser.exit(status=1,
                        message='\n' + '#' * 40 +
                        '\nInclusion levels can only be between 0 and 1\n' +
                        '#' * 40 + '\n')


if __name__ == '__main__':
    # command line arguments
    parser = argparse.ArgumentParser(description='Command line interface for designing primers')
    group_one = parser.add_mutually_exclusive_group(required=True)
    group_one.add_argument('-b', dest='big_bed', action='store', help='big bed file that defines the possible exons in a gene')
    group_one.add_argument('-g', dest='gtf', action='store', help='gtf file that defines the possible exons in a gene')
    parser.add_argument('--no-gene-id', dest='no_gene_id', action='store_true', help='Use this flag if your gtf does not have a valid gene_id')
    parser.add_argument('-f', required=True, dest='fasta', action='store', help='path to fasta file')
    parser.add_argument('-r', required=True, dest='rnaseq', action=ValidateRnaseq, help='path to SAM/BAM file(s) ("," delimited)')
    parser.add_argument('-t', required=True, dest='target', action='store', help='path to txt file with <strand><chr>:<start>-<end> for each target on separate lines.')
    group_two = parser.add_mutually_exclusive_group(required=True)
    group_two.add_argument('--annotaton', dest='annotation_flag', action='store_true', help='only use junctions supported from annotation')
    group_two.add_argument('--rnaseq', dest='rnaseq_flag', action='store_true', help='only use junctions supported from RNA-Seq')
    group_two.add_argument('--both', dest='both_flag', action='store_true', help='use junctions from both RNA-Seq and annotation')
    parser.add_argument('--psi', dest='psi', action=ValidateCutoff, default=1.0, type=float, help='Define inclusion level sufficient to define constitutive exon. Valid: 0<psi<1.')
    parser.add_argument('--read-threshold', dest='read_threshold', default=5, action='store', type=int, help='Define the minimum number of read support necessary to call a junction from RNA-Seq')
    parser.add_argument('--keep-temp', dest='keep_temp', action='store_true', help='Keep temporary files in your tmp directory')
    parser.add_argument('-o', required=True, dest='output', action='store', help='Output directory')
    options = vars(parser.parse_args())  # make it a dictionary

    # define job_id by the name of the target file
    tmp = options['target'].split('/\\')[-1].split('.')
    options['job_id'] = ('.'.join(tmp[:-1]) if len(tmp) > 1 else tmp[0]) + '.output'
    with open(options['target']) as handle:
        options['target'] = map(lambda x: x.strip().split('\t'), handle.readlines())

    # define logging file before using logging.debug
    if not os.path.exists('log'): os.mkdir('log')  # make directory to put log files
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(message)s',
                        filename='log/log.primer.' + str(datetime.datetime.now()),
                        filemode='w')

    ### Start loading the user's files ###
    # gtf file must be pre-loaded since there is no random access
    if options['gtf']:
        print 'Loading GTF . . .'
        print 'May take ~1 min.'
        options['gtf'] = gene_annotation_reader(options['gtf'])

    print 'Loading fasta . . .'
    options['fasta'] = SequenceFileDB(options['fasta'])  # get fasta object using pygr right away

    # the sam object interfaces with the user specified BAM/SAM file!!!
    print 'Loading Bam Files . . .'
    options['rnaseq'] = [sam.Sam(data) for data in options['rnaseq'].split(',')]
    print 'Done loading all files.'
    ### END loading files ###

    # call main function
    main(options)
