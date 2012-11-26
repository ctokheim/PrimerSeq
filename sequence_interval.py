import csv
import subprocess


class SequenceInterval(object):

    def __init__(self, exons, psi_list, reverse=False):
        self.exons = exons
        self.psi_list = psi_list
        self.reverse = reverse
        self.tmp_input_file = "input.txt"
        self.path_to_sequence_interval_maps = "bin"
        self.create_input_file()  # make interval for sequence_interval_maps
        self.run_sequence_interval_maps()
        self.result_intervals, self.result_psi = self.read_sequence_interval_maps_output()

    def create_input_file(self):
        """
        Create input file for sequence_interval_maps.cpp
        """
        with open(self.tmp_input_file, 'wb') as handle:
            my_writer = csv.writer(handle, delimiter=",")
            tmp_data = []
            for i, (e_start, e_end) in enumerate(self.exons):
                tmp_data.append([e_start, e_end, self.psi_list[i]])
            my_writer.writerows(tmp_data)

    def run_sequence_interval_maps(self):
        """
        Run command line call to sequence_interval_maps.cpp binary
        """
        cmd = "%s/sequence_interval_maps" % self.path_to_sequence_interval_maps
        subprocess.check_call(cmd, shell=True)

    def read_sequence_interval_maps_output(self):
        """
        Read the output from sequence_interval_maps.cpp
        """
        with open("output.txt") as handle:
            result = list(csv.reader(handle, delimiter="\t"))
            result.sort(key=lambda x: (x[0], x[1]), reverse=self.reverse)  # make sure intervals are sorted by position
            iv = [(int(row[0]), int(row[1])) for row in result]
            psi = [float(row[2]) for row in result]
        return iv, psi

    def get_results(self):
        """

        """
        return self.result_intervals, self.result_psi
