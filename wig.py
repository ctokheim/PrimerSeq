#!/usr/bin/env python
# Copyright (C) 2012-2013  Collin Tokheim
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

from base_bed_wig import BaseBedWig
import csv
from collections import OrderedDict
import logging


class Wig(BaseBedWig):
    """
    Reads in depth information from big wig file
    """
    def load_wig_file(self):
        logging.debug('Extracting wig file from ' + self.current_file + ' . . .')
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
        logging.debug('Finished extracting wig file.')
