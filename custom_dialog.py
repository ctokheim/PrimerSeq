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

# Note to self, pubsub import must immediately follow wx
import wx
from wx.lib.pubsub import setuparg1
from wx.lib.pubsub import pub

import custom_thread as ct
import csv
import os
import utils
import sys
import primer
import logging
import json
import draw
import depth_plot
import subprocess
import gtf
import webbrowser
import ConfigParser
import re
import splice_graph as sg
from exon_seek import ExonSeek
from save_plots_html import SavePlotsHTML
import shutil
import algorithms as algs

import datetime  # import for marking html output with date
import traceback  # debugging import


class CustomDialog(wx.Dialog):
    def __init__(self, parent, id, title, text=''):
        wx.Dialog.__init__(self, parent, id, title, size=(300, 100))

        self.parent = parent
        self.text = wx.StaticText(self, -1, text)
        self.empty_text = wx.StaticText(self, -1, '')

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.empty_text, 0, wx.ALIGN_CENTER)
        sizer.Add(self.text, 0, wx.ALIGN_CENTER)
        sizer.Add(self.empty_text, 0, wx.ALIGN_CENTER)

        self.SetSizer(sizer)
        self.Show()

    def Update(self, val, update_text=''):
        if val == 100:
            self.Destroy()
        else:
            pass

    def check_dialog(self):
        pass


class PlotDialog(wx.Dialog):
    def __init__(self, parent, id, title, output_file, text=''):
        wx.Dialog.__init__(self, parent, id, title, size=(300, 100), style=wx.DEFAULT_DIALOG_STYLE)

        self.output_file = output_file

        self.parent = parent
        self.text = wx.StaticText(self, -1, text)

        self.bigwig_label = wx.StaticText(self, -1, "BigWig(s):")
        self.choose_bigwig_button = wx.Button(self, -1, "Choose . . .")
        self.bigwig = []
        self.panel_3 = wx.Panel(self, -1)
        # self.bigwig_choice_label = wx.StaticText(self, -1, "None")
        self.bigwig_choice_label = wx.TextCtrl(self, -1, "None", style=wx.TE_READONLY)
        bigwig_sizer = wx.GridSizer(1, 3, 0, 0)
        bigwig_sizer.Add(self.bigwig_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        bigwig_sizer.Add(self.choose_bigwig_button, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 0)
        bigwig_sizer.Add(self.bigwig_choice_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)

        # read in valid primer output
        with open(self.output_file) as handle:
            self.results = filter(lambda x: len(x) > 1,  # if there is no tabs then it represents an error msg in the output
                                  csv.reader(handle, delimiter='\t'))[1:]
            select_results = [', '.join([r[0], r[-1], r[1]]) for r in self.results]

        # target selection widgets
        target_sizer = wx.GridSizer(1, 2, 0, 0)
        self.target_label = wx.StaticText(self, -1, "Select Target:")
        self.target_combo_box = wx.ComboBox(self, -1, choices=select_results, style=wx.CB_DROPDOWN | wx.TE_READONLY)
        self.target_combo_box.SetMinSize((145, 27))
        target_sizer.Add(self.target_label, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        target_sizer.Add(self.target_combo_box, 0, wx.ALIGN_LEFT, 0)

        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.plot_button = wx.Button(self, -1, 'Plot')
        self.cancel_button = wx.Button(self, -1, 'Cancel')
        button_sizer.Add(self.plot_button, 0, wx.ALIGN_RIGHT)
        button_sizer.Add(self.cancel_button, 0, wx.ALIGN_LEFT)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(bigwig_sizer, 0, wx.EXPAND, 10)
        sizer.Add(target_sizer, 0, wx.EXPAND)
        sizer.Add(button_sizer, 0, wx.ALIGN_CENTER)
        sizer.SetMinSize((300, 100))

        self.Bind(wx.EVT_BUTTON, self.choose_bigwig_event, self.choose_bigwig_button)
        self.Bind(wx.EVT_BUTTON, self.plot_button_event, self.plot_button)
        self.Bind(wx.EVT_BUTTON, self.cancel_button_event, self.cancel_button)
        self.SetSizer(sizer)
        self.Show()

        pub.subscribe(self.plot_update, "plot_update")
        pub.subscribe(self.on_plot_error, "plot_error")

    def on_plot_error(self, msg):
        # self.parent.update_after_error((None,))
        pass  # prevent error on trying to call method

    def cancel_button_event(self, event):
        self.Destroy()
        event.Skip()

    def choose_bigwig_event(self, event):
        dlg = wx.FileDialog(self, message='Choose your BigWig files', defaultDir=os.getcwd(),
                            wildcard='BigWig files (*.bw)|*.bw|BigWig files (*.bigWig)|*.bigWig', style=wx.FD_MULTIPLE)  # open file dialog
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            filenames = dlg.GetPaths()  # get the new filenames from the dialog
            filenames_without_path = dlg.GetFilenames()  # only grab the actual filenames and none of the path information
            dlg.Destroy()  # best to do this sooner

            self.bigwig = filenames
            # self.bigwig_choice_label.SetLabel(', '.join(filenames_without_path))
            self.bigwig_choice_label.SetValue(', '.join(filenames_without_path))
        else:
            dlg.Destroy()

    def plot_button_event(self, event):
        if not self.target_combo_box.GetValue():
            dlg = wx.MessageDialog(self, 'Please select a BigWig file and the target exon\nyou want to plot.', style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            return

        self.target_id = str(self.target_combo_box.GetValue().split(',')[0])
        self.target_of_interest = str(self.target_combo_box.GetValue().split(', ')[1])

        # get the line from the file that matches the user selection
        for row in self.results:
            if row[0] == self.target_id:
                row_of_interest = row

        # find where the plot should span
        start, end = utils.get_pos(row_of_interest[-2])
        chr = utils.get_chr(row_of_interest[-2])
        plot_domain = utils.construct_coordinate(chr, start, end)
        gene_name = row_of_interest[-1]  # currently gene name is last column but likely will change
        self.target_pos = utils.get_pos(row_of_interest[1])

        self.plot_button.SetLabel('Ploting . . .')
        self.plot_button.Disable()

        # draw isoforms
        plot_thread = ct.PlotThread(target=self.generate_plots,
                                    args=(self.target_id,
                                          self.target_pos,
                                          plot_domain,
                                          self.bigwig,
                                          self.output_file,
                                          gene_name))

    def generate_plots(self, tgt_id, target_pos, plt_domain, bigwig, out_file, gene_name):
        # generate isoform drawing
        opts = {'json': primer.config_options['tmp'] + '/isoforms/' + tgt_id + '.json',
                'output': primer.config_options['tmp'] + '/' + 'draw/' + tgt_id + '.png',
                'target_exon': target_pos,
                'scale': 1,
                'primer_file': out_file,
                'id': tgt_id}
        self.draw_isoforms(opts)

        # generate read depth plot
        opts = {'bigwig': ','.join(bigwig),
                'position': plt_domain,
                'gene': gene_name,
                'size': 2.,
                'step': 1,
                'output': primer.config_options['tmp'] + '/depth_plot/' + tgt_id + '.png'}
        self.depth_plot(opts)

    def draw_isoforms(self, opts):
        '''
        Draw isoforms by using draw.py
        '''
        logging.debug('Drawing isoforms %s . . .' % str(opts))
        # load json file that has information isoforms and their counts
        with open(opts['json']) as handle:
            my_json = json.load(handle)

        coord = draw.read_primer_file(self.output_file, opts['id'])
        draw.main(my_json['path'], opts['target_exon'], my_json['counts'], coord, opts)
        logging.debug('Finished drawing isoforms.')

    def depth_plot(self, opts):
        '''
        Create a read depth plot by using depth_plot.py
        '''
        logging.debug('Creating read depth plot %s . . .' % str(opts))
        depth_plot.read_depth_plot(opts)
        logging.debug('Finished creating read depth plot.')

    def plot_update(self, msg):
        self.plot_button.SetLabel('Plot')
        self.plot_button.Enable()
        DisplayPlotDialog(self, -1, 'Primer Results for ' + self.target_of_interest,
                          ['tmp/depth_plot/' + self.target_id + '.png',
                           'tmp/draw/' + self.target_id + '.png'])


class SortGtfDialog(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(300, 100), style=wx.DEFAULT_DIALOG_STYLE)

        self.parent = parent

        self.gtf_label = wx.StaticText(self, -1, "GTF:")
        self.choose_gtf_button = wx.Button(self, -1, "Choose . . .")
        self.panel_3 = wx.Panel(self, -1)
        self.gtf_choice_label = wx.StaticText(self, -1, "None")
        gtf_sizer = wx.GridSizer(1, 3, 0, 0)
        gtf_sizer.Add(self.gtf_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        gtf_sizer.Add(self.choose_gtf_button, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 0)
        gtf_sizer.Add(self.gtf_choice_label, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)
        self.output_gtf_label = wx.StaticText(self, -1, "Sorted GTF:")
        self.choose_output_gtf_button = wx.Button(self, -1, "Choose . . .")
        self.panel_3 = wx.Panel(self, -1)
        self.output_gtf_choice_label = wx.StaticText(self, -1, "None")
        output_gtf_sizer = wx.GridSizer(1, 3, 0, 0)
        output_gtf_sizer.Add(self.output_gtf_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        output_gtf_sizer.Add(self.choose_output_gtf_button, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 0)
        output_gtf_sizer.Add(self.output_gtf_choice_label, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)

        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sort_button = wx.Button(self, -1, 'Sort')
        self.cancel_button = wx.Button(self, -1, 'Cancel')
        button_sizer.Add(self.sort_button, 0, wx.ALIGN_RIGHT)
        button_sizer.Add(self.cancel_button, 0, wx.ALIGN_LEFT)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(gtf_sizer, 0, wx.EXPAND, 10)
        sizer.Add(output_gtf_sizer, 0, wx.EXPAND)
        sizer.Add(button_sizer, 0, wx.ALIGN_CENTER)
        sizer.SetMinSize((300, 100))

        self.Bind(wx.EVT_BUTTON, self.choose_gtf_event, self.choose_gtf_button)
        self.Bind(wx.EVT_BUTTON, self.choose_output_gtf_event, self.choose_output_gtf_button)
        self.Bind(wx.EVT_BUTTON, self.sort_button_event, self.sort_button)
        self.Bind(wx.EVT_BUTTON, self.cancel_button_event, self.cancel_button)
        self.SetSizer(sizer)
        self.Show()

        pub.subscribe(self.sort_update, "sort_update")
        pub.subscribe(self.sort_error, "update_after_error")

    def sort_error(self, msg):
        self.parent.update_after_error((None,))

    def cancel_button_event(self, event):
        self.Destroy()
        event.Skip()

    def choose_output_gtf_event(self, event):
        dlg = wx.FileDialog(self, message='Choose your GTF file to be sorted', defaultDir=os.getcwd(),
                            wildcard='GTF file (*.gtf)|*.gtf')  # open file dialog
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()  # get the new filenames from the dialog
            filename_without_path = dlg.GetFilename()  # only grab the actual filenames and none of the path information
            dlg.Destroy()  # best to do this sooner

            self.output_gtf = filename
            self.output_gtf_choice_label.SetLabel(filename_without_path)
        else:
            dlg.Destroy()

    def choose_gtf_event(self, event):
        dlg = wx.FileDialog(self, message='Choose your GTF file to be sorted', defaultDir=os.getcwd(),
                            wildcard='GTF file (*.gtf)|*.gtf')  # open file dialog
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()  # get the new filenames from the dialog
            filename_without_path = dlg.GetFilename()  # only grab the actual filenames and none of the path information
            dlg.Destroy()  # best to do this sooner

            self.gtf = filename
            self.gtf_choice_label.SetLabel(filename_without_path)
        else:
            dlg.Destroy()

    def sort_gtf(self, infile, outfile):
        """Sort a GTF file using SortGtf.jar"""
        my_config = ConfigParser.ConfigParser()
        my_config.read('PrimerSeq.cfg')
        config_options = dict(my_config.items('memory'))
        cmd = 'java -jar -Xmx%sm "bin/SortGtf.jar" "%s" "%s"' % (config_options['sort'], infile, outfile)
        logging.debug('Sort GTF cmd: ' + cmd)
        subprocess.check_call(cmd, shell=True)

    def sort_button_event(self, event):
        self.sort_button.SetLabel('Sorting . . .')
        self.sort_button.Disable()

        # draw isoforms
        sort_thread = ct.UpdateThread(target=self.sort_gtf,
                                      args=(self.gtf, self.output_gtf),
                                      update='sort_update')

    def sort_update(self, msg):
        self.sort_button.SetLabel('Sort')
        self.sort_button.Enable()


class AddGeneIdsDialog(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(300, 100), style=wx.DEFAULT_DIALOG_STYLE)

        self.parent = parent

        self.gtf_label = wx.StaticText(self, -1, "GTF:")
        self.choose_gtf_button = wx.Button(self, -1, "Choose . . .")
        self.panel_3 = wx.Panel(self, -1)
        self.gtf_choice_label = wx.StaticText(self, -1, "None")
        gtf_sizer = wx.GridSizer(1, 3, 0, 0)
        gtf_sizer.Add(self.gtf_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        gtf_sizer.Add(self.choose_gtf_button, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 0)
        gtf_sizer.Add(self.gtf_choice_label, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)
        self.kgxref_label = wx.StaticText(self, -1, "kgXref:")
        self.choose_kgxref_button = wx.Button(self, -1, "Choose . . .")
        self.panel_3 = wx.Panel(self, -1)
        self.kgxref_choice_label = wx.StaticText(self, -1, "None")
        kgxref_sizer = wx.GridSizer(1, 3, 0, 0)
        kgxref_sizer.Add(self.kgxref_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        kgxref_sizer.Add(self.choose_kgxref_button, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 0)
        kgxref_sizer.Add(self.kgxref_choice_label, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)
        self.output_gtf_label = wx.StaticText(self, -1, "GTF W/ Genes:")
        self.choose_output_gtf_button = wx.Button(self, -1, "Choose . . .")
        self.panel_3 = wx.Panel(self, -1)
        self.output_gtf_choice_label = wx.StaticText(self, -1, "None")
        output_gtf_sizer = wx.GridSizer(1, 3, 0, 0)
        output_gtf_sizer.Add(self.output_gtf_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        output_gtf_sizer.Add(self.choose_output_gtf_button, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 0)
        output_gtf_sizer.Add(self.output_gtf_choice_label, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)

        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.add_genes_button = wx.Button(self, -1, 'Change Gene IDs')
        self.cancel_button = wx.Button(self, -1, 'Cancel')
        button_sizer.Add(self.add_genes_button, 0, wx.ALIGN_RIGHT)
        button_sizer.Add(self.cancel_button, 0, wx.ALIGN_LEFT)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(gtf_sizer, 0, wx.EXPAND, 10)
        sizer.Add(kgxref_sizer, 0, wx.EXPAND, 10)
        sizer.Add(output_gtf_sizer, 0, wx.EXPAND)
        sizer.Add(button_sizer, 0, wx.ALIGN_CENTER)
        sizer.SetMinSize((300, 100))

        self.Bind(wx.EVT_BUTTON, self.choose_gtf_event, self.choose_gtf_button)
        self.Bind(wx.EVT_BUTTON, self.choose_kgxref_event, self.choose_kgxref_button)
        self.Bind(wx.EVT_BUTTON, self.choose_output_gtf_event, self.choose_output_gtf_button)
        self.Bind(wx.EVT_BUTTON, self.add_genes_button_event, self.add_genes_button)
        self.Bind(wx.EVT_BUTTON, self.cancel_button_event, self.cancel_button)
        self.SetSizerAndFit(sizer)
        # self.SetSizer(sizer)
        self.Show()

        pub.subscribe(self.add_gene_ids_update, "add_update")

    def cancel_button_event(self, event):
        self.Destroy()
        event.Skip()

    def choose_output_gtf_event(self, event):
        dlg = wx.FileDialog(self, message='Choose your GTF file to be sorted', defaultDir=os.getcwd(),
                            wildcard='GTF file (*.gtf)|*.gtf')  # open file dialog
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()  # get the new filenames from the dialog
            filename_without_path = dlg.GetFilename()  # only grab the actual filenames and none of the path information
            dlg.Destroy()  # best to do this sooner

            self.output_gtf = filename
            self.output_gtf_choice_label.SetLabel(filename_without_path)
        else:
            dlg.Destroy()

    def choose_gtf_event(self, event):
        dlg = wx.FileDialog(self, message='Choose your GTF file to be sorted', defaultDir=os.getcwd(),
                            wildcard='GTF file (*.gtf)|*.gtf')  # open file dialog
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()  # get the new filenames from the dialog
            filename_without_path = dlg.GetFilename()  # only grab the actual filenames and none of the path information
            dlg.Destroy()  # best to do this sooner

            self.gtf = filename
            self.gtf_choice_label.SetLabel(filename_without_path)
        else:
            dlg.Destroy()

    def choose_kgxref_event(self, event):
        dlg = wx.FileDialog(self, message='Choose kgxref txt file', defaultDir=os.getcwd(),
                            wildcard='txt file (*.txt)|*.txt')  # open file dialog
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()  # get the new filenames from the dialog
            filename_without_path = dlg.GetFilename()  # only grab the actual filenames and none of the path information
            dlg.Destroy()  # best to do this sooner

            self.kgxref = filename
            self.kgxref_choice_label.SetLabel(filename_without_path)
        else:
            dlg.Destroy()

    def add_genes_button_event(self, event):
        self.add_genes_button.SetLabel('Adding . . .')
        self.add_genes_button.Disable()

        opts = {'annotation': self.gtf,
                'kgxref': self.kgxref,
                'output': self.output_gtf}

        # draw isoforms
        gene_thread = ct.UpdateThread(target=gn.main,
                                      args=(opts,),
                                      update='add_update')

    def add_gene_ids_update(self, msg):
        self.add_genes_button.SetLabel('Change Gene IDs')
        self.add_genes_button.Enable()


class InSilicoPcrDialog(wx.Dialog):
    def __init__(self, parent, id, title, output_file):
        wx.Dialog.__init__(self, parent, id, title,
                           size=(300, 175))
                           #style=wx.DEFAULT_DIALOG_STYLE)

        self.parent = parent
        self.output_file = output_file

        # read in valid primer output
        with open(self.output_file) as handle:
            self.results = filter(lambda x: len(x) > 1,  # if there is no tabs then it represents an error msg in the output
                                  csv.reader(handle, delimiter='\t'))[1:]
            select_results = [', '.join(r[:2]) for r in self.results]

        self.genome_label = wx.StaticText(self, -1, "Genome:  ")
        self.panel_3 = wx.Panel(self, -1)
        self.genome_text_field = wx.TextCtrl(self, -1, "Human")
        genome_sizer = wx.GridSizer(1, 2, 0, 0)
        genome_sizer.Add(self.genome_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        genome_sizer.Add(self.genome_text_field, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)
        self.assembly_label = wx.StaticText(self, -1, "Assembly:  ")
        self.assembly_text_field = wx.TextCtrl(self, -1, "hg19")
        assembly_sizer = wx.GridSizer(1, 2, 0, 0)
        assembly_sizer.Add(self.assembly_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        assembly_sizer.Add(self.assembly_text_field, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)
        self.type_label = wx.StaticText(self, -1, "Select Type:  ")
        self.type_combo_box = wx.ComboBox(self, -1, 'UCSC Genes', choices=['UCSC Genes', 'Genome'], style=wx.CB_DROPDOWN | wx.TE_READONLY)
        self.type_combo_box.SetMinSize((145, 27))
        type_sizer = wx.GridSizer(1, 2, 0, 0)
        type_sizer.Add(self.type_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        type_sizer.Add(self.type_combo_box, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)
        self.max_prod_size_label = wx.StaticText(self, -1, "Max Product size:  ")
        self.max_prod_size_text_field = wx.TextCtrl(self, -1, "4000")
        max_prod_size_sizer = wx.GridSizer(1, 2, 0, 0)
        max_prod_size_sizer.Add(self.max_prod_size_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        max_prod_size_sizer.Add(self.max_prod_size_text_field, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)
        self.target_label = wx.StaticText(self, -1, "Select Target:  ")
        self.target_combo_box = wx.ComboBox(self, -1, select_results[0], choices=select_results, style=wx.CB_DROPDOWN | wx.TE_READONLY)
        self.target_combo_box.SetMinSize((145, 27))
        target_sizer = wx.GridSizer(1, 2, 0, 0)
        target_sizer.Add(self.target_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        target_sizer.Add(self.target_combo_box, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 0)

        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.run_button = wx.Button(self, -1, 'Run In-Silico PCR')
        self.cancel_button = wx.Button(self, -1, 'Cancel')
        button_sizer.Add(self.run_button, 0, wx.ALIGN_RIGHT)
        button_sizer.Add(self.cancel_button, 0, wx.ALIGN_LEFT)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(genome_sizer, 0, wx.EXPAND, 10)
        sizer.Add(assembly_sizer, 0, wx.EXPAND)
        sizer.Add(type_sizer, 0, wx.EXPAND)
        sizer.Add(max_prod_size_sizer, 0, wx.EXPAND)
        sizer.Add(target_sizer, 0, wx.EXPAND)
        sizer.Add(button_sizer, 0, wx.ALIGN_CENTER)
        sizer.SetMinSize((300, 100))

        self.Bind(wx.EVT_BUTTON, self.on_run, self.run_button)
        self.Bind(wx.EVT_BUTTON, self.on_cancel, self.cancel_button)
        self.SetSizer(sizer)
        self.Show()

    def on_cancel(self, event):
        self.Destroy()
        event.Skip()

    def on_run(self, event):
        # check user input, alert user if missing data
        if not self.target_combo_box.GetValue() or not self.type_combo_box.GetValue or not self.genome_text_field.GetValue() or not self.max_prod_size_text_field.GetValue() or not self.assembly_text_field.GetValue():
            dlg = wx.MessageDialog(self, 'Please fill in all input fields.', style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            return

        # find result which matches user selection
        user_selected_id = str(self.target_combo_box.GetValue().split(',')[0])
        for row in self.results:
            if row[0] == user_selected_id:
                upstream_seq = row[4]
                downstream_seq = row[5]

        # construct url
        ucsc_url = utils.InSilicoPcrUrl(genome=str(self.genome_text_field.GetValue()),
                                        assembly=str(self.assembly_text_field.GetValue()),
                                        forward=upstream_seq,
                                        reverse=downstream_seq,
                                        target=str(self.type_combo_box.GetValue()),
                                        max_size=int(self.max_prod_size_text_field.GetValue()))
        webbrowser.open(ucsc_url.get_url())  # open url in webbrowser


class DisplayPlotDialog(wx.Dialog):
    def __init__(self, parent, id, title, img_files):
        # call super constructor
        wx.Dialog.__init__(self, parent, id, title,)
                           # style=wx.DEFAULT_DIALOG_STYLE ^ wx.RESIZE_BORDER)

        # containers for imgs
        depth_png = wx.Image(img_files[0], wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        draw_png = wx.Image(img_files[1], wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        self.draw_bitmap = wx.StaticBitmap(self, -1, draw_png, (10, 5), (draw_png.GetWidth(), draw_png.GetHeight()))
        self.depth_bitmap = wx.StaticBitmap(self, -1, depth_png, (10, 5), (depth_png.GetWidth(), depth_png.GetHeight()))

        self.parent = parent

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.depth_bitmap, 0, wx.ALIGN_CENTER)
        sizer.Add(self.draw_bitmap, 0, wx.ALIGN_CENTER)

        self.SetSizerAndFit(sizer)
        self.Show()


class SavePlotDialog(wx.Dialog):
    def __init__(self, parent, id, title, opts, text=''):
        wx.Dialog.__init__(self, parent, id, title, size=(400, 270), style=wx.DEFAULT_DIALOG_STYLE)

        self.options = opts
        self.output_file = opts['output']
        self.output_directory = None  # user needs to specify an output directory

        self.parent = parent
        self.text = wx.StaticText(self, -1, text)

        # self.data_label = wx.StaticText(self, -1, "Title,BAM,BigWig:")
        # self.data_label.SetFont(wx.Font(14, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        # self.data_text_field = wx.TextCtrl(self, -1, ''.join(map(lambda x: x + ',\n', ['Title for data,' + s.path for s in self.options['rnaseq']])), style=wx.TE_MULTILINE | wx.HSCROLL)
        # self.data_text_field.SetToolTip(wx.ToolTip("eg. Title,mySample.bam,mySample.bw"))
        # self.data_text_field.SetMinSize((396, 60))
        # self.data_text_field.SetFocus()  # set focus on text field

        tID = wx.NewId()
        self.list = utils.MyListCtrl(self, tID,
                             style=wx.LC_REPORT
                             | wx.LC_SORT_ASCENDING
                             | wx.LC_EDIT_LABELS
                             | wx.BORDER_NONE)
        # define the columns
        self.list.InsertColumn(0, 'Title')
        self.list.InsertColumn(1, 'BAM')
        self.list.InsertColumn(2, 'BigWig')
        self.set_list(self.options['rnaseq'])  # populate the list ctrl with data

        # read in valid primer output
        with open(self.output_file) as handle:
            self.total_results = list(csv.reader(handle, delimiter='\t'))[1:]  # contains designed primers and failed cases
            self.results = filter(lambda x: len(x) > 1,  self.total_results)  # if there is no tabs then it represents an error msg in the output

        # widgets for handling choice of output directory
        grid_sizer = wx.GridSizer(1, 3, 0, 0)
        self.directory_label = wx.StaticText(self, -1, "Directory:")
        self.directory_label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        self.choose_directory_button = wx.Button(self, -1, "Choose . . .")
        self.choose_directory_button.SetToolTip(wx.ToolTip('Choose your output directory'))
        self.directory_choice_label = wx.TextCtrl(self, -1, "None", style=wx.TE_READONLY)
        self.directory_choice_label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        grid_sizer.Add(self.directory_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)
        grid_sizer.Add(self.choose_directory_button, 0, wx.ALIGN_CENTER, 10)
        grid_sizer.Add(self.directory_choice_label, 0, wx.EXPAND, 10)

        # create genome widgets
        grid_sizer_genome = wx.GridSizer(1, 2, 0, 0)
        self.genome_label = wx.StaticText(self, -1, "Genome:  ")
        self.genome_label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        self.genome_choice_label = wx.TextCtrl(self, -1, "Human")
        self.genome_choice_label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        grid_sizer_genome.Add(self.genome_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)
        grid_sizer_genome.Add(self.genome_choice_label, 0, wx.EXPAND, 10)

        # create genome assembly widgets
        grid_sizer_assembly = wx.GridSizer(1, 2, 0, 0)
        self.assembly_label = wx.StaticText(self, -1, "Assembly:  ")
        self.assembly_label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        self.assembly_choice_label = wx.TextCtrl(self, -1, "hg19")
        self.assembly_choice_label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        grid_sizer_assembly.Add(self.assembly_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)
        grid_sizer_assembly.Add(self.assembly_choice_label, 0, wx.EXPAND, 10)

        # run and cancel buttons
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.save_plot_button = wx.Button(self, -1, 'Generate Report')
        self.cancel_button = wx.Button(self, -1, 'Cancel')
        button_sizer.Add(self.save_plot_button, 0, wx.ALIGN_RIGHT)
        button_sizer.Add(self.cancel_button, 0, wx.ALIGN_LEFT)

        # set sizer information
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.AddMany([(self.list, 0, wx.EXPAND, 10),
                       # (self.data_label, 0, wx.ALIGN_CENTER, 10),  # add label
                       # ((10, 10), 0),  # add spacer
                       # (self.data_text_field, 0, wx.EXPAND, 10),  # add text box
                       ((10, 10), 0),  # add spacer
                       (grid_sizer_genome, 0, wx.EXPAND, 10),
                       ((10, 10), 0),  # add spacer
                       (grid_sizer_assembly, 0, wx.EXPAND, 10),
                       ((10, 10), 0),  # add spacer
                       (grid_sizer, 0, wx.EXPAND, 10),  # add output directory choice
                       ((10, 10), 0),  # add spacer
                       (button_sizer, 0, wx.ALIGN_CENTER)])  # add button
        sizer.SetMinSize((500, 100))

        self.Bind(wx.EVT_BUTTON, self.on_save_plot, self.save_plot_button)
        self.Bind(wx.EVT_BUTTON, self.cancel_button_event, self.cancel_button)
        self.Bind(wx.EVT_BUTTON, self.on_directory_choice, self.choose_directory_button)
        self.SetSizer(sizer)
        self.Show()

        pub.subscribe(self.on_finish, "plotting_finished")

    def set_list(self, bam_list):
        """Add bam information to listctrl"""
        for i, bam in enumerate(bam_list):
            index = self.list.InsertStringItem(sys.maxint, bam.path)
            self.list.SetStringItem(index, 0, 'Title for Data')
            self.list.SetStringItem(index, 1, bam.path)
            self.list.SetStringItem(index, 2, '')
            self.list.SetItemData(index, i)

    def on_finish(self, msg):
        try:
            # testing showed this message dialog threw an exception sometimes on windows vista
            dlg = wx.MessageDialog(self, 'Finished! A webbrowser may open or open a new tab in an existing browser.', style=wx.OK)
            dlg.ShowModal()
        except:
            # wierd behaviour has some times cause the dialog to fail to display
            logging.debug('Traceback:\n' + traceback.format_exc())
        finally:
            # always reset buttons and open browser
            self.save_plot_button.SetLabel('Generate Report')
            self.save_plot_button.Enable()
            webbrowser.open(os.path.join(self.output_directory, 'index.html'))

    def on_save_plot(self, event):
        if not self.output_directory:
            dlg = wx.MessageDialog(self, 'Please enter an output directory.', style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            return  # exit if user didn't specify an output

        # get information from listctrl and make sure user specified data
        counts = self.list.GetItemCount()
        missing_bigwig_flag = False  # flag for absent bigwig file
        titles, bigwigs = [], []
        for row in xrange(counts):
            title = self.list.GetItem(itemId=row, col=0).GetText()
            bigwig = self.list.GetItem(itemId=row, col=2).GetText()
            titles.append(title)
            bigwigs.append(bigwig)
            if not bigwig:
                missing_bigwig_flag = True
        if missing_bigwig_flag:
            dlg = wx.MessageDialog(self, 'One or more BigWig files were not specified. If you intended to'
                                   ' not plot read depth then press OK. If you want to plot read depth then'
                                   ' press CANCEL and then type the file paths for the BigWig files corresponding'
                                   ' to the BAM file.', 'Confirm BigWig', wx.OK | wx.CANCEL | wx.ICON_QUESTION)
            result = dlg.ShowModal()
            dlg.Destroy()
            if result == wx.ID_CANCEL:
                return  # exit if user wanted to specify bigwig files

        self.generate_index_html()  # create the index.html web page

        # lines = map(lambda y: re.split('\s*,+\s*', y), map(str, filter(lambda x: x != '', re.split('\s*\n+\s*', self.data_text_field.GetValue()))))
        # bigwig_files = map(lambda x: x[2], lines)
        # titles = map(lambda x: x[0], lines)
        html_thread = ct.HtmlReportThread(self.generate_plots, args=(self.options, bigwigs, self.output_directory, titles))

        # disable button so they cannot press it twice
        self.save_plot_button.SetLabel('Generating . . .')
        self.save_plot_button.Disable()

    def generate_index_html(self):
        # start creating a index.html web page
        index_html = SavePlotsHTML()
        index_html.add_heading('Alternative Splicing Events')
        index_html.add_text('%d/%d primer designs are successful! ' % (len(self.results), len(self.total_results)))

        # add text file to directory and link from html to it
        shutil.copy(self.output_file, os.path.join(self.output_directory, 'output.txt'))  # copy the text file
        index_html.add_text('You can view the entire results ')
        index_html.add_link('output.txt', 'here')  # add link to text file with results
        index_html.add_line_break()
        index_html.add_line_break()

        # add links to AS events with designed primers
        for line in self.results:
            index_html.add_text(line[0] + ', ')
            index_html.add_link(line[0] + '.html',
                                line[1])  # link to each AS event
            index_html.add_text(', <i>' + line[-1] + '</i>')  # italicized gene name

            # construct url
            index_html.add_text(' -- ')
            ucsc_url = utils.InSilicoPcrUrl(genome=str(self.genome_choice_label.GetValue()),
                                            assembly=str(self.assembly_choice_label.GetValue()),
                                            forward=line[4],
                                            reverse=line[5],
                                            target='UCSC Genes',
                                            max_size=4000)
            index_html.add_link(ucsc_url.get_url(), '<i>In-Silico</i> PCR')  # open url in webbrowser
            index_html.add_line_break()  # make each link on separate line

        # Create time stamp
        index_html.add_line_break()
        index_html.add_text('Timestamp: ' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))

        # write list of links to the index.html file
        with open(os.path.join(self.output_directory, 'index.html'), 'w') as handle:
            handle.write(str(index_html))

    def cancel_button_event(self, event):
        self.Destroy()
        event.Skip()

    def on_directory_choice(self, event):
        dlg = wx.DirDialog(self, message='Choose an output directory')
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            self.output_directory = dlg.GetPath()  # get the new filenames from the dialog
            # self.directory_choice_label.SetLabel(self.output_directory)
            self.directory_choice_label.SetValue(self.output_directory)
            dlg.Destroy()  # best to do this sooner
        else:
            dlg.Destroy()

    def generate_plots(self, options, bigwigs, out_dir, titles):
        try:
            handle = open(options['output'], 'r')
            handle.readline()
            for line in csv.reader(handle, delimiter='\t'):
                if len(line) <= 1: continue  # skip cases where no good output
                ID = line[0]
                start, end = utils.get_pos(line[-2])  # ASM region column
                chr = utils.get_chr(line[-2])  # ASM region column
                tgt_pos = utils.get_pos(line[1])
                plot_domain = utils.construct_coordinate(chr, start, end)
                # only error msgs and blank lines do not have tabs
                if len(line) > 1:
                    tx_paths, counts, gene = self.get_isforms_and_counts(line, options)
                    my_html = SavePlotsHTML()
                    for index in range(len(tx_paths)):
                        path, count = tx_paths[index], counts[index]
                        self.create_plots(ID, index, plot_domain, path, tgt_pos, count, [bigwigs[index]], gene, options['output'], out_dir)
                        my_html.add_heading(titles[index])
                        my_html.add_img('%s.%d.depth.png' % (ID, index))
                        my_html.add_img('%s.%d.isoforms.png' % (ID, index))
                        my_html.add_line_break()
                    with open(os.path.join(out_dir, ID + '.html'), 'wb') as html_writer:
                        html_writer.write(str(my_html))
            handle.close()
            shutil.copy('style.css', out_dir)  # copy css to folder
        except:
            print traceback.format_exc()

    def get_isforms_and_counts(self, line, options):
        # get information about each row
        ID, target_coordinate = line[:2]
        strand = target_coordinate[0]
        chr = utils.get_chr(target_coordinate[1:])
        tmp_start, tmp_end = utils.get_pos(target_coordinate)

        # get information regarding the gene
        if options['no_gene_id']:
            gene_dict, gene_name = sg.get_weakly_connected_tx(options['gtf'], strand, chr, tmp_start, tmp_end)  # hopefully filter out junk
        else:
            gene_dict, gene_name = sg.get_from_gtf_using_gene_name(options['gtf'], strand, chr, tmp_start, tmp_end)

        # get edge weights
        edge_weights_list = [sam_obj.extractSamRegion(chr, gene_dict['start'], gene_dict['end'])
                             for sam_obj in options['rnaseq']]

        # construct splice graph for each BAM file
        bam_splice_graphs = sg.construct_splice_graph(edge_weights_list,
                                                      gene_dict,
                                                      chr,
                                                      strand,
                                                      options['read_threshold'],
                                                      options['min_jct_count'],
                                                      output_type='list',
                                                      both=options['both_flag'])

        paths_list = []
        counts_list = []
        for my_splice_graph in bam_splice_graphs:
            # not the best use of the ExonSeek object, initially intended to find appropriate flanking exons
            # but in this case ExonSeek is used to get the transcripts and associate counts
            exon_seek_obj = ExonSeek(utils.get_pos(target_coordinate), my_splice_graph, ID, options['psi'], None, None)
            all_paths, upstream, downstream, component, psi_target, psi_upstream, psi_downstream = exon_seek_obj.get_info()
            paths_list.append(exon_seek_obj.paths)
            counts_list.append(exon_seek_obj.counts)
        return paths_list, counts_list, gene_name  # return the tx paths and count information for a single AS event

    def create_plots(self, tgt_id, bam_index, plt_domain, path, tgt_exon, counts, bigwig, gene, out_file, output_directory):
        """Create plots for a single BAM/BigWig pair"""
        # generate isoform drawing
        opts = {'path': path,
                'target_exon': tgt_exon,
                'counts': counts,
                'output': os.path.join(output_directory, tgt_id + '.' + str(bam_index) + '.isoforms.png'),
                'scale': 1,
                'primer_file': out_file,
                'id': tgt_id}
        self.draw_isoforms(opts)

        # generate read depth plot
        opts = {'bigwig': ','.join(bigwig),
                'position': plt_domain,
                'gene': gene,
                'size': 2.,
                'step': 1,
                'output': os.path.join(output_directory, tgt_id + '.' + str(bam_index) + '.depth.png')}
        self.depth_plot(opts)

    def draw_isoforms(self, opts):
        '''
        Draw isoforms by using draw.py
        '''
        logging.debug('Drawing isoforms %s . . .' % str(opts))
        coord = draw.read_primer_file(self.output_file, opts['id'])
        draw.main(opts['path'], opts['target_exon'], opts['counts'], coord, opts)
        logging.debug('Finished drawing isoforms.')

    def depth_plot(self, opts):
        '''
        Create a read depth plot by using depth_plot.py
        '''
        logging.debug('Creating read depth plot %s . . .' % str(opts))
        depth_plot.read_depth_plot(opts)
        logging.debug('Finished creating read depth plot.')

    def plot_update(self, msg):
        self.plot_button.SetLabel('Plot')
        self.plot_button.Enable()
        DisplayPlotDialog(self, -1, 'Primer Results for ' + self.target_of_interest,
                          ['tmp/depth_plot/' + self.target_id + '.png',
                           'tmp/draw/' + self.target_id + '.png'])
