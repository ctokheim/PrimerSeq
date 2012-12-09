from base_bed_wig import BaseBedWig
import csv
from collections import OrderedDict


class Wig(BaseBedWig):
    """
    Reads in depth information from big wig file
    """
    def load_wig_file(self):
        # depth_dict = {}
        depth_dict = OrderedDict()
        self.bins = []  # defines the histogram bins for plotting
        with open(self.current_file) as handle:
            for line in csv.reader(handle, delimiter='\t'):
                # self.verbatim_content.append(line)
                # for i in range(int(line[1]), int(line[2])):
                interval_start, interval_end = int(line[1]), int(line[2])
                depth_dict[(interval_start, interval_end)] = int(line[3])
                self.bins += [interval_start, interval_end]
        self.annotation = depth_dict
        self.bins = sorted(list(set(self.bins)))
