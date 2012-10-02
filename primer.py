import subprocess
import os
import shutil
import gtf
import csv
import argparse

def primer3(options, primer3_options):
    """
    The primer.py main function uses the gtf module to find information about constitutive flanking exons for the target exons of interest.
    It then designs primers by calling primer3. Next it parses the primer3 output and outputs the final results to a file. The output file
    is then emailed to the designated address in the command line parameters. NOTE: the email command is hard coded for the windows server
    computer.
    """
    jobs_ID = args.job_ID
    Exon_name = map(lambda x: x.lstrip().rstrip(), args.Exon_name.split(","))
    flanking_info = gtf.main(full_target_info, pathToGtf)
    outfiles_PATH = 'tmp'   # directory for intermediate files

    # iterate over all target sequences
    STRAND, EXON_TARGET, UPSTREAM_TARGET, DOWNSTREAM_TARGET, INC_LIST, SKIP_LIST, UPSTREAM_Seq, TARGET_SEQ, DOWNSTREAM_SEQ = range(9)
    output_list = []
    for z in range(len(flanking_info)):
        # no flanking exon information case
        if type('') == type(flanking_info[z]):
            output_list.append(flanking_info[z] + '\n')  # write problem msg
        # has flanking exon information case
        else:
            ####################### Primer3 Parameter Configuration###########
            PRIMER_FILE_FLAG = '1'
            PRIMER_EXPLAIN_FLAG = '1'
            PRIMER_SEQUENCE_ID = jobs_ID
            SEQUENCE = flanking_info[z][UPSTREAM_Seq] + flanking_info[z][TARGET_SEQ].lower() + flanking_info[z][DOWNSTREAM_SEQ]
            TARGET = str(len(flanking_info[z][UPSTREAM_Seq]) + 1) + ',' + str(len(flanking_info[z][TARGET_SEQ]))
            #############################################################

            ####################### Write jobs_ID.conf##################
            outfile = open(outfiles_PATH + jobs_ID + '.conf', 'w')

            # hard coded options
            outfile.write('PRIMER_SEQUENCE_ID=' + PRIMER_SEQUENCE_ID + '\n')
            outfile.write('SEQUENCE=' + SEQUENCE + '\n')
            outfile.write('TARGET=' + TARGET + '\n')
            outfile.write('PRIMER_FILE_FLAG=' + PRIMER_FILE_FLAG + '\n')
            outfile.write('PRIMER_EXPLAIN_FLAG=' + PRIMER_EXPLAIN_FLAG + '\n')

            # options from primer3.cfg
            for o in primer3_options:
                outfile.write(o)

            outfile.write('=' + '\n')  # primer3 likes a '=' at the end of sequence params
            outfile.close()

            ###################### Primer3 #####################################
            os.remove(outfiles_PATH + jobs_ID + '.Primer3')  # delete old files

            cmd = 'primer3_core.exe < ' + outfiles_PATH + jobs_ID + '.conf' + ' > ' + outfiles_PATH + jobs_ID + '.Primer3'
            subprocess.check_call(cmd, shell=True)

            #################### Parse '.Primer3' ################################
            # read primer3 file
            infile = open(outfiles_PATH + jobs_ID + '.Primer3', 'r')
            w_Primer3 = {}
            for line in infile:
                try:
                    w_Primer3[line.split('=')[0]] = line.split('=')[1].rstrip()
                except IndexError:
                    print line
            infile.close()

            # checks if no output
            if(w_Primer3.keys().count('PRIMER_LEFT_SEQUENCE') == 0):
                output_list.append('Failed: No_Primer3_results for ' + Exon_name[z] + '\n')
                infile.close()
                continue
            # there is output case
            else:
                # get info about product sizes
                target_exon_len = len(flanking_info[z][TARGET_SEQ])
                Primer3_PRIMER_PRODUCT_SIZE = int(w_Primer3['PRIMER_PRODUCT_SIZE']) - target_exon_len
                skipping_size_list = [str(int(j) + Primer3_PRIMER_PRODUCT_SIZE) for j in flanking_info[z][SKIP_LIST]]
                inclusion_size_list = [str(int(k) + Primer3_PRIMER_PRODUCT_SIZE) for k in flanking_info[z][INC_LIST]]
                skipping_size = ';'.join(skipping_size_list)
                inclusion_size = ';'.join(inclusion_size_list)

                # append results to output_list
                tmp = [Exon_name[z], flanking_info[z][EXON_TARGET], w_Primer3['PRIMER_LEFT_SEQUENCE'], w_Primer3['PRIMER_RIGHT_SEQUENCE'],
                       str((float(w_Primer3['PRIMER_LEFT_TM']) + float(w_Primer3['PRIMER_RIGHT_TM'])) / 2), skipping_size, inclusion_size,
                       flanking_info[z][UPSTREAM_TARGET], flanking_info[z][DOWNSTREAM_TARGET]]
                output_list.append(tmp)
                infile.close()

    # write output information
    with open(outfiles_PATH + jobs_ID + '.txt', 'wb') as outputfile_tab:
        # define csv header
        header = ['Exon_name', 'Target_Exon_coordinates', 'Left_primer', 'Right_primer', 'Average_Tm',
                  'Skipping_production_size', 'Inclusion_production_size', 'Upstream_coordinates', 'Downstream_coordinates']
        output_list = [header] + output_list  # pre-pend header to output file
        csv.writer(outputfile_tab, delimiter='\t').writerows(output_list)  # output primer design to a tab delimited file


def main(options):
    """Reads in primer3 options and then calls run function to run primer3"""
    # read in primer3 options
    primer3_options = []
    with open('primer3.cfg') as handle:
        for line in handle:
            # skip comment lines and lines with no set values
            if not line.startswith('#') and len(line.strip().split("=")[1]) > 0:
                primer3_options.append(line)

    # the primer3 function runs the primer3_core executable
    primer3(options, primer3_options)

if __name__ == '__main__':
    # command line arguments
    parser = argparse.ArgumentParser(description='Command line interface program for designing primers')
    parser.add_argument('-g', required=True, dest='gtf', action='store', help='gtf file that defines the possible exons in a gene')
    parser.add_argument('-f', required=True, dest='fasta', action='store', help='path to fasta file')
    parser.add_argument('-r', required=True, dest='rnaseq', action='store', help='path to RNA-Seq data')
    parser.add_argument('-t', required=True, dest='target', action='store', help='path to txt file with <strand><chr>:<start>-<end> for each target.')
    options = vars(parser.parse_args())

    # call main function
    main(options)
