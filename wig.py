from base_bed_wig import BaseBedWig
import csv


class Wig(BaseBedWig):
    """
    Reads in depth information from big wig file
    """
    def load_wig_file(self):
        depth_dict = {}
        self.verbatim_content = []
        print self.current_file
        with open(self.current_file) as handle:
            for line in csv.reader(handle, delimiter='\t'):
                self.verbatim_content.append(line)
                for i in range(int(line[1]), int(line[2])):
                    depth_dict[i] = int(line[3])
        self.annotation = depth_dict
