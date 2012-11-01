'''
File: sam.py
Author: Collin Tokheim
Description: Python wrapper around ExtractSamRegion.jar and Convert2SortedBam.jar
'''

import subprocess
import jct_counts
import os
import ConfigParser

# for logging purposes
import logging
import traceback
import sys
import time

# define directories
cfg = ConfigParser.ConfigParser()
cfg.read('PrimerSeq.cfg')
cfg_options = dict(cfg.items('directory'))
TMP_DIR = cfg_options['tmp'] + '/sam/'  # directory to dump intermediate sam files
BIN_DIR = 'bin/'  # directory of jar files
JCT_DIR = cfg_options['tmp'] + '/jct/'  # directory for .jct files


class Sam(object):
    """
    Designed to be the interface with the SAM/BAM files. It automatically converts
    to sorted bam and allows for region extraction.

    NOTE: This class is a wrapper arround Convert2SortedBam.jar and ExtractSamRegion.jar
    thus it requires the JRE.
    """
    def __init__(self, sam_path):
        # complain about not ending with .sam/.bam
        if not sam_path.endswith('.sam') and not sam_path.endswith('.bam'):
            raise ValueError('RNA-Seq input should be in SAM or BAM format')

        # skip if named .sorted.bam
        if sam_path.endswith('.sorted.bam'):
            self.path = sam_path
        # skip if created .sorted.bam before
        elif os.path.exists(sam_path[:-4] + '.sorted.bam'):
            self.path = sam_path[:-4] + '.sorted.bam'
        # call Convert2SortedBam.jar if no .sorted.bam
        else:
            self.path = self.convert2SortedBam(sam_path)

    def convert2SortedBam(self, myPath):
        """
        Calls Convert2SortedBam.jar to index/sort the bam file to allow random access
        to specific regions in the BAM file
        """
        try:
            logging.debug('Converting %s to a sorted BAM file . . .' % myPath)
            beginTime = time.time()
            sorted_bam_path = myPath[:-4] + '.sorted.bam'
            cmd = 'java -jar -Xmx1024m %sConvert2SortedBam.jar %s %s' % (BIN_DIR, myPath, sorted_bam_path)
            subprocess.check_call(cmd, shell=True)
            endTime = time.time()
            runTime = endTime - beginTime
            logging.debug('Finished converting %s to %s in %.2d:%.2d:%.2d' % (
                myPath, sorted_bam_path, runTime / 3600, (runTime % 3600) / 60, runTime % 60))
            return sorted_bam_path
        except subprocess.CalledProcessError:
            t, v, trace = sys.exc_info()  # exception information
            logging.debug('Error! Call to Convert2SortedBam.jar with non-zero exit status')
            logging.debug('Type: ' + str(t))
            logging.debug('Value: ' + str(v))
            logging.debug('Traceback:\n' + traceback.format_exc())
            raise

    def __jct_to_dict(self, file_name, jct_dict={}):
        """Reads jct counts created from jctCounts.py"""
        with open(file_name) as handle:
            for line in handle:
                lineSplit = line.split('\t')
                chr, start, stop, count = lineSplit  # rename line columns something meaningful
                start, stop, count = int(start), int(stop), int(count)  # convert to int

                # add counts to jct_dict
                #jct_dict.setdefault(chr, {})
                #jct_dict[chr].setdefault(start, {})
                #jct_dict[chr][start].setdefault(stop, 0)
                #jct_dict[chr][start][stop] += count
                jct_dict[(chr, start, stop)] = count  # this simplifies the key comparted to a multi-level dict

        return jct_dict

    def __get_sam_jct(self, samFile):
        """
        Read in SAM junctions
        """
        logging.debug('Reading junctions from %s' % samFile)
        jctOutputFile = JCT_DIR + samFile.split('/')[-1][:-4] + '.jct'  # use .jct file for jct count files

        # specify options for jct_counts.py main function
        options = {}
        options['sam'] = samFile
        options['output'] = jctOutputFile
        options['anchor'] = 8

        # get reads
        jct_counts.main(options)  # note, it creates a file at jctOuputFile
        jctDict = self.__jct_to_dict(jctOutputFile)  # creates a dictionary object for jcts
        return jctDict

    def extractSamRegion(self, chr, start, end):
        """
        Retrieve jct counts from a specified region in the BAM file (self.path)
        """
        try:
            start += 1  # extraction is done in 1-based coordinates
            logging.debug('Extracting reads for %s:%d-%d' % (chr, start, end))
            tmp_sam_path = '%s%s_%d_%d.sam' % (TMP_DIR, chr, start, end)  # path to tmp sam file with region specific reads
            cmd = 'java -jar -Xmx1024m %sExtractSamRegion.jar %s %s %s %d %d' % (
                BIN_DIR, self.path, tmp_sam_path, chr, start, end)
            subprocess.check_call(cmd, shell=True)
            logging.debug('Finished getting sam reads. Parsing jcts . . .')
            junctionDict = self.__get_sam_jct(tmp_sam_path)
            logging.debug('Finished reading jcts')
            return junctionDict
        except subprocess.CalledProcessError:
            t, v, trace = sys.exc_info()  # exception information
            logging.debug('Error! Extracting sam reads using ExtractSamRegion.jar failed')
            logging.debug('Type: ' + str(t))
            logging.debug('Value: ' + str(v))
            logging.debug('Traceback:\n' + traceback.format_exc())
            raise



