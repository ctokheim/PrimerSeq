#!/usr/bin/env python
# Copyright (C) 2012  Collin Tokheim
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

import csv
from base_bed_wig import BaseBedWig


class Bed(BaseBedWig):

    def load_bed_file(self):
        genes, exons = set(), set()
        self.verbatim_content = []
        self.annotation = {}
        self.annotation['start'], self.annotation['end'] = float('inf'), 0
        self.annotation['chr'], self.annotation['strand'] = self.chr, self.strand
        with open(self.current_file) as handle:
            self.annotation['graph'] = []
            for line in csv.reader(handle, delimiter='\t'):
                # ignore lines that do not match the target strand or chr
                if line[0] == self.chr and line[5] == self.strand:
                    self.verbatim_content.append(line)  # store all content just in case
                    genes.add(line[3])
                    chr, tx_start = line[0], int(line[1])
                    exon_lengths = map(int, line[-2].split(',')[:-1])
                    exon_offsets = map(int, line[-1].split(',')[:-1])
                    exon_start = [e + tx_start for e in exon_offsets]
                    exon_end = [e + exon_lengths[i] for i, e in enumerate(exon_start)]
                    self.annotation['start'] = min(self.annotation['start'], exon_start[0])
                    self.annotation['end'] = max(self.annotation['end'], exon_end[-1])
                    tx_exons = zip(exon_start, exon_end)
                    self.annotation['graph'].append(tx_exons)
                    exons |= set(tx_exons)

            self.annotation['exons'] = sorted(list(exons), key=lambda x: (x[0], x[1]))  # record all unique exons

            # make assertions that only one gene is assumed to match target
            assert len(genes) != 0, 'Your coordinates did not match anything in the annotation'
            assert len(genes) == 1, 'There is multiple genes overlapping the target. Perhaps the "name" column in the bed file is actually transcript names rather than gene names'

            self.annotation['gene_name'] = list(genes)[0]  # assign gene name after knowing its unique

            # find the exon which the target is contained
            is_contained = False
            for strt, end in self.annotation['exons']:
                if self.start >= strt and self.end <= end:
                    is_contained = True
                    self.annotation['target'] = (strt, end)
            assert is_contained is True, "Expected target to be completely within an exon"
