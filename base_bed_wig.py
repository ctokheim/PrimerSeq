import logging
import sys
import subprocess
import traceback


class BaseBedWig(object):
    """
    This class is meant to serve as the base class for the Bed
    class (bed.py) and the Wig class (wig.py). Both call the same
    ExtractBigRegion.jar file.
    """
    def __init__(self, bbfile):
        # expected file paths
        self.TMP_DIR = 'tmp'
        self.BIN_DIR = 'bin'

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

    def extractBigRegion(self, strand, chr, start, end):
        try:
            logging.debug('Extracting bed lines overlaping %s:%d-%d' % (chr, start, end))
            cmd = 'java -jar %s/ExtractBigRegion.jar %s %s/bed/%s_%d_%d.bed %s %d %d true' % (
                self.BIN_DIR, self.bbfile, self.TMP_DIR, chr, start, end, chr, start, end)
            subprocess.check_call(cmd, shell=True)
            self.current_file = self.TMP_DIR + '/bed/%s_%d_%d.bed'  # path to current bed file the class is working on
            self.strand, self.chr, self.start, self.end = strand, chr, start, end  # hold on to target information just in case
            logging.debug('Finished extracting bed lines for %s:%d-%d' % (chr, start, end))
        except:
            t, v, trace = sys.exc_info()
            logging.debug('ERROR! Call to ExtractBigRegion.jar with non-zero exit status')
            logging.debug('Type: ' + str(t))
            logging.debug('Value: ' + str(v))
            logging.debug('Traceback:\n' + traceback.format_exc())
            raise

    def get_annotation(self):
        return self.annotation
