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

import wx
import wx.lib.mixins.listctrl as listmix
from bisect import bisect  # packaged used by TextEditMixin
import json


# define column order of output
ID = 0
TARGET = 1
PRIMER_COORD = 2
PSI_TGT = 3
FOR_PRIMER = 4
REV_PRIMER = 5
AVG_TM = 6
SKIP_SIZE = 7
INC_SIZE = 8
UPSTREAM_EXON = 9
PSI_UP = 10
DOWNSTREAM_EXON = 11
PSI_DOWN = 12
ASM_REGION = 13
GENE = 14


def get_start_pos(coordinate):
    """
    get start from 'chr:start-stop'
    """
    return int(coordinate.split(":")[1].split("-")[0])


def get_end_pos(coordinate):
    """
    get stop from 'chr:start-stop'
    """
    return int(coordinate.split(":")[1].split("-")[1])


def get_pos(coordinate):
    """
    get (start, stop) from 'chr:start-stop'
    """
    return tuple(map(int, coordinate.split(":")[1].split("-")))


def get_chr(coordinate):
    """
    get chr from 'chr:start-stop'
    """
    return coordinate.split(":")[0]


def get_primer_coordinates(primer_coord):
    """Turns chr:start1-end1;chr:start2-end2 into (start1, end1), (start2, end2)"""
    first, second = primer_coord.split(";")
    return get_pos(first), get_pos(second)


def get_path_range(paths):
    """Returns the genomic range of several paths"""
    tmp = []
    for p in paths:
        for exon in p:
            tmp.append(exon[0])
            tmp.append(exon[1])
    return min(tmp), max(tmp)


def find_first_exon(pos, exon_list):
    """
    Finds the first exon the position is contained within.
    Assumes exon_list is already sorted.
    """
    for i, (ex_start, ex_end) in enumerate(exon_list):
        if pos[0] >= ex_start and pos[1] <= ex_end:
            return i


def find_last_exon(pos, exon_list):
    """
    Finds the last exon the position is contained within.
    Assumes exon_list is already sorted.
    """
    for i in range(len(exon_list)-1, -1, -1):
        ex_start, ex_end = exon_list[i]
        if pos[0] >= ex_start and pos[1] <= ex_end:
            return i


def construct_coordinate(chr, start, end):
    return '%s:%s-%s' % (chr, str(start), str(end))


def merge_list_of_dicts(list_of_dicts):
    '''
    This function mereges multiple dicts contained read counts from SAM/BAM file
    into one dictionary.
    '''
    merged_dict = {}
    for tmp_dict in list_of_dicts:
        all_keys = set(merged_dict) | set(tmp_dict)
        for key in all_keys:
            merged_dict[key] = merged_dict.get(key, 0) + tmp_dict.get(key, 0)
    return merged_dict


def calc_product_length(path, primer_coord):
    """
    Calculate product length based on the primer coordinates
    """
    # calculate length between primers
    tmp_len = 0
    first_exon_primer, second_exon_primer = False, False  # flag for telling if a primer was contained within an exon
    flag = False
    for start, end in path:
        if start <= primer_coord[0][0] and end >= primer_coord[0][1] and start <= primer_coord[1][0] and end >= primer_coord[1][1]:
            tmp_len = primer_coord[1][0] - primer_coord[0][1]
            first_exon_primer = True
            second_exon_primer = True
        elif start <= primer_coord[0][0] and end >= primer_coord[0][1]:
            tmp_len += end - primer_coord[0][1]
            flag = True
            first_exon_primer = True
        elif start <= primer_coord[1][0] and end >= primer_coord[1][1]:
            tmp_len += primer_coord[1][0] - start
            flag = False
            second_exon_primer = True
        elif flag:
            tmp_len += end - start

    if not (first_exon_primer and second_exon_primer):
        # case where either one or both primers were not located within an exon
        final_len = 0
    else:
        # add length between primers to actual length of sequnces for each primer
        first_primer_len = primer_coord[0][1] - primer_coord[0][0]
        second_primer_len = primer_coord[1][1] - primer_coord[1][0]
        final_len = tmp_len + first_primer_len + second_primer_len
    return final_len


class PrimerSeqError(Exception):
    """
    Used to uniquely catch this exception so that only one primer
    fails to be designed
    """
    pass


class InSilicoPcrUrl(object):
    '''
    Construct url for In-Silico PCR. This class constructs all of
    the GET submission parameters.
    '''
    def __init__(self, genome='', assembly='', target='',
                 forward='', reverse='', max_size=4000, perfect=15,
                 good=15, flip=0):
        # initialize variables
        self.genome = genome
        self.assembly = assembly
        self.forward = forward
        self.reverse = reverse
        self.max_size = max_size
        self.perfect = perfect
        self.good = good
        self.flip = flip
        if target == 'Genome':
            self.target = 'genome'
        elif target == 'UCSC Genes':
            # unfortunately the ucsc genes have diff version numbers
            # for different assemblies. this means the GET parameter input
            # for wp_target may change based on the species/assembly.
            assemblyToTxSuffix = {
                'hg19': 'v14',
                'hg18': 'Nov08',
                'mm10': 'v13'
            }
            suffix = assemblyToTxSuffix[self.assembly] if self.assembly in assemblyToTxSuffix else ''
            self.target = self.assembly + 'Kg' + suffix

        # base url for ucsc in-silico pcr
        self.base_url = 'http://genome.ucsc.edu/cgi-bin/hgPcr'

    def get_url(self):
        get_params = '?org=%s&db=%s&wp_target=%s&wp_f=%s&wp_r=%s&wp_size=%s&wp_perfect=%s&wp_good=%s&boolshad.wp_flipReverse=%s&Submit=submit' % (
                            self.genome, self.assembly, self.target, self.forward, self.reverse, str(self.max_size), str(self.perfect), str(self.good), str(self.flip))
        return self.base_url + get_params


class MyListCtrl(listmix.ListCtrlAutoWidthMixin, wx.ListCtrl, listmix.TextEditMixin):
    def __init__(self, p, my_id, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, p, my_id,
                             pos, (-1, 125), style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftClick)
        self.Bind(wx.EVT_LEFT_DCLICK, self.OnLeftClick)

    def OnLeftClick(self, evt=None):
        col = self.get_col(evt)

        if not col:
            evt.Skip()
            return

        colText = self.GetColumn(col).GetText()
        if colText == "BigWig":
            # open file dialog for BigWig file selection
            dlg = wx.FileDialog(self, message='Choose your BigWig file',
                    wildcard='BigWig file (*.bw)|*.bw|BigWig file(*.bigwig)|*.bigwig')  # open file dialog

            # if they press ok
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()  # get the new filenames from the dialog
                myitem = self.GetItem(self.curRow, col)  # get item of interest
                myitem.SetText(filename)
                self.SetItem(myitem)
                dlg.Destroy()  # best to do this sooner

            # dlg = wx.MessageDialog(self, 'Row: %d, Col: %d\n%s' % (self.curRow, col, self.GetItem(self.curRow, col).GetText()), style=wx.OK)
            # dlg.ShowModal()
        elif colText == "BAM":
            # make BAM column uneditable
            return
        else:
            super(MyListCtrl, self).OnLeftDown(evt)

    def get_col(self, evt):
        x, y = evt.GetPosition()
        row, flags = self.HitTest((x, y))
        if row != self.curRow:
            return None

        self.col_locs = [0]
        loc = 0
        for n in range(self.GetColumnCount()):
            loc = loc + self.GetColumnWidth(n)
            self.col_locs.append(loc)

        col = bisect(self.col_locs, x+self.GetScrollPos(wx.HORIZONTAL)) - 1

        return col




def save_path_info(file_basename, p, cts, save_dir='tmp/isoforms/'):
    '''
    Save information about isoforms and their read counts into a json file.
    '''
    with open(save_dir + str(file_basename) + '.json', 'w') as handle:
        json.dump({'path': p, 'counts': list(cts)}, handle, indent=4)  # output path information to tmp file


def get_isoforms_and_counts(file_path):
    """Get path and count info after using save_path_info"""
    with open(file_path) as handle:
        data = json.load(handle)
    return data['path'], data['counts']

def get_seq_from_list(chr_seq, strand, pos_list):
    """
    chr_seq: pygr chromosome object (can slice to get sequence)
    strand: '+' or '-'
    pos_list: [(start1, end1), (start2, end2), ...]
    Returns the sequence of each position as a concatenated string that
    goes upstream -> downstream.
    """
    sorted_pos = sorted(pos_list, key=lambda x: (x[0], x[1]))
    if strand == '-':
        sorted_pos = reversed(sorted_pos)

    seq = ''
    for start, end in sorted_pos:
        if strand == '+':
            seq += str(chr_seq[start:end]).lower()
        else:
            seq += str(-chr_seq[start:end]).lower()
    return seq
