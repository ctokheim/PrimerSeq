import logging
import sys
import subprocess
import traceback
import os
import ConfigParser


class BaseBedWig(object):
    """
    This class is meant to serve as the base class for the Bed
    class (bed.py) and the Wig class (wig.py). Both call the same
    ExtractBigRegion.jar file.
    """
    def __init__(self, bbfile, ext='bed'):
        # expected file paths
        cfg = ConfigParser.ConfigParser()
        cfg.read('PrimerSeq.cfg')
        cfg_options = dict(cfg.items('directory'))
        self.TMP_DIR = cfg_options['tmp']
        self.BIN_DIR = 'bin'
        self.ext = ext  # the file extension (either "bed" or "wig")
        if not os.path.exists(self.TMP_DIR + '/' + self.ext):
            os.mkdir(self.TMP_DIR + '/' + self.ext)  # mkdir in case it doesn't exist

        self.bbfile = bbfile
        self.current_file = None
        self.verbatim_content = None

        # the annotation attribute is what is of interest to other parts of my
        # program. That is, if everything ran correctly the user just wants to
        # get the self.annotation value at the end.
        self.annotation = None

        # useful variables about the target position
        self.strand = None
        self.chr = None
        self.start = None
        self.end = None

    def extractBigRegion(self, chr, start, end, strand='+'):
        try:
            logging.debug('Extracting %s lines overlaping %s:%d-%d' % (self.ext, chr, start, end))
            cmd = 'java -jar -Xmx512m "%s/ExtractBigRegion.jar" "%s" "%s/%s/%s_%d_%d.%s" %s %d %d true' % (
                self.BIN_DIR, self.bbfile, self.TMP_DIR, self.ext, chr, start, end, self.ext, chr, start, end)
            logging.debug('CMD is [%s]' % cmd)
            subprocess.check_call(cmd, shell=True)  # call to ExtractBigRegion.jar
            self.strand, self.chr, self.start, self.end = strand, chr, start, end  # hold on to target information just in case
            self.current_file = self.TMP_DIR + '/%s/%s_%d_%d.%s' % (self.ext, self.chr, self.start, self.end, self.ext)  # path to current bed file the class is working on
            logging.debug('Finished extracting %s lines for %s:%d-%d' % (self.ext, chr, start, end))
        except subprocess.CalledProcessError:
            t, v, trace = sys.exc_info()
            logging.debug('ERROR! Call to ExtractBigRegion.jar with non-zero exit status')
            logging.debug('Type: ' + str(t))
            logging.debug('Value: ' + str(v))
            logging.debug('Traceback:\n' + traceback.format_exc())
            raise

    def get_annotation(self):
        return self.annotation
