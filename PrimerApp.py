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

# these imports are to prevent import errors when I distribute the code
import anydbm  # probably not necessary since using pyinstaller
import dbhash  # probably not necessary since using pyinstaller

# useful imports
import wx
# from wx.lib.pubsub import Publisher
from wx.lib.pubsub import setuparg1
from wx.lib.pubsub import pub
import wx.stc as Stc
import os
import glob
import subprocess
import re
import sys
import csv
import shutil
from pygr.seqdb import SequenceFileDB
import sam
import primer
import webbrowser
import custom_thread as ct
import custom_dialog as cd
import view_output as vo
import read_counts as rc
import utils

# logging imports
import traceback
import logging
import datetime


def handle_uncaught_exceptions(t, ex, tb):
    """
    Should catch all uncaught exceptions so that an alert dialog can be
    displayed to the user so that they know an error happened
    """
    # traceback.print_tb(tb)  # print traceback to stdout so I can debug
    logging.debug('Error: ' + str(t) + str(ex))
    logging.debug('Traceback:\n' + ''.join(traceback.format_list(traceback.extract_tb(tb))))
    dlg = wx.MessageDialog(None,
                           'An uncaught error occured in PrimerSeq.'
                           'Please check the log file (%s) for details.'
                           ' You may need to press File->Reset to continue.' % log_file,
                           style=wx.OK | wx.ICON_ERROR)
    dlg.ShowModal()


class PrimerFrame(wx.Frame):
    """Top level frame in my application."""

    def __init__(self, *args, **kwds):
        # begin wxGlade: PrimerFrame.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE ^ wx.MAXIMIZE_BOX
        wx.Frame.__init__(self, *args, **kwds)
        self.SetSizeHints(10, 10, maxH=500)

        # Menu Bar
        self.primer_frame_menubar = wx.MenuBar()
        wxglade_tmp_menu = wx.Menu()
        load_quick_ex_id = wx.NewId()
        wxglade_tmp_menu.Append(load_quick_ex_id, "&Load Quick Ex.", "Load a single exon example data", wx.ITEM_NORMAL)
        load_ex_id = wx.NewId()
        wxglade_tmp_menu.Append(load_ex_id, "&Load Larger Ex.", "Load a 10 exon example data set", wx.ITEM_NORMAL)
        reset_id = wx.NewId()
        wxglade_tmp_menu.Append(reset_id, "&Reset", "Reset all inputs", wx.ITEM_NORMAL)
        quit_id = wx.NewId()
        wxglade_tmp_menu.Append(quit_id, "&Quit", "Terminate the program", wx.ITEM_NORMAL)
        self.primer_frame_menubar.Append(wxglade_tmp_menu, "File")
        wxglade_tmp_menu = wx.Menu()
        self.primer_frame_menubar.Append(wxglade_tmp_menu, "Edit")
        primer3_id = wx.NewId()
        wxglade_tmp_menu.Append(primer3_id, "&Primer3", "Edit the parameters used by Primer3 to design primers", wx.ITEM_NORMAL)
        path_id = wx.NewId()
        wxglade_tmp_menu.Append(path_id, "Primer3 Location", "Edit where primer3 and the primer3 configuraton file are located", wx.ITEM_NORMAL)
        sort_id = wx.NewId()
        wxglade_tmp_menu.Append(sort_id, "&Sort GTF", "Sort a GTF file to make it a proper input for PrimerSeq", wx.ITEM_NORMAL)
        add_genes_id = wx.NewId()
        wxglade_tmp_menu.Append(add_genes_id, "Add &Genes", "Add gene IDs to UCSCs gene annotation", wx.ITEM_NORMAL)
        wxglade_tmp_menu = wx.Menu()
        help_id = wx.NewId()
        wxglade_tmp_menu.Append(help_id, "Help", "Open help information", wx.ITEM_NORMAL)
        primer3_manual_id = wx.NewId()
        wxglade_tmp_menu.Append(primer3_manual_id, "Primer3 &Doc.", "Primer3 manual", wx.ITEM_NORMAL)
        about_id = wx.NewId()
        wxglade_tmp_menu.Append(about_id, "&About", "Information regarding PrimerSeq", wx.ITEM_NORMAL)
        self.primer_frame_menubar.Append(wxglade_tmp_menu, "Help")
        wxglade_tmp_menu = wx.Menu()
        self.SetMenuBar(self.primer_frame_menubar)
        # Menu Bar end
        self.primer_frame_statusbar = self.CreateStatusBar(1, 0)
        self.primer_notebook = wx.Notebook(self, -1, style=0)
        self.primer_notebook_pane_1 = wx.Panel(self.primer_notebook, -1)
        self.fasta_label = wx.StaticText(self.primer_notebook_pane_1, -1, "FASTA:")
        self.choose_fasta_button = wx.Button(self.primer_notebook_pane_1, -1, "Choose . . .")
        self.choose_fasta_button.SetToolTip(wx.ToolTip('Select your genome sequence\nin FASTA format'))
        self.fasta_choice_label = wx.TextCtrl(self.primer_notebook_pane_1, -1, "None", style=wx.NO_BORDER | wx.TE_READONLY)
        self.fasta_choice_label.SetBackgroundColour(self.primer_notebook.GetBackgroundColour())
        self.gtf_label = wx.StaticText(self.primer_notebook_pane_1, -1, "GTF:")
        self.choose_gtf_button = wx.Button(self.primer_notebook_pane_1, -1, "Choose . . .")
        self.choose_gtf_button.SetToolTip(wx.ToolTip('Select your gene annotation\nin a sorted GTF format'))
        self.gtf_choice_label = wx.TextCtrl(self.primer_notebook_pane_1, -1, "None", style=wx.NO_BORDER | wx.TE_READONLY)
        self.gtf_choice_label.SetBackgroundColour(self.primer_notebook.GetBackgroundColour())
        self.bam_label = wx.StaticText(self.primer_notebook_pane_1, -1, "SAM/BAM(s):")
        self.choose_bam_button = wx.Button(self.primer_notebook_pane_1, -1, "Choose . . .")
        self.choose_bam_button.SetToolTip(wx.ToolTip('Select one or multiple SAM or BAM file(s).\nWhen selecting, hold ctrl to select multiple'))
        self.bam_choice_label = wx.TextCtrl(self.primer_notebook_pane_1, -1, "None", style=wx.NO_BORDER | wx.TE_READONLY)
        self.bam_choice_label.SetBackgroundColour(self.primer_notebook.GetBackgroundColour())
        self.sizer_4_staticbox = wx.StaticBox(self.primer_notebook_pane_1, -1, "Load Files")
        self.coordinates_label = wx.StaticText(self.primer_notebook_pane_1, -1, "Coordinates:")
        self.coordinates_text_field = wx.TextCtrl(self.primer_notebook_pane_1, -1, "", style=wx.TE_MULTILINE | wx.HSCROLL)
        self.coordinates_text_field.SetToolTip(wx.ToolTip("(strand)chr:start-end\n0-based start, 1-based end"))
        self.output_label = wx.StaticText(self.primer_notebook_pane_1, -1, "Output:")
        self.choose_output_button = wx.Button(self.primer_notebook_pane_1, -1, "Choose . . .")
        self.choose_output_button.SetToolTip(wx.ToolTip('Select the output file'))
        # self.panel_4 = wx.Panel(self.primer_notebook_pane_1, -1)
        # self.output_choice_label = wx.TextCtrl(self.panel_4, -1, "None", style=wx.TE_READONLY)
        self.output_choice_label = wx.TextCtrl(self.primer_notebook_pane_1, -1, "None", style=wx.NO_BORDER | wx.TE_READONLY)
        self.output_choice_label.SetBackgroundColour(self.primer_notebook.GetBackgroundColour())
        self.run_button = wx.Button(self.primer_notebook_pane_1, -1, "Run PrimerSeq")
        self.run_button.SetToolTip(wx.ToolTip('Run PrimerSeq'))
        self.primer_notebook_pane_2 = wx.Panel(self.primer_notebook, -1)
        self.psi_label = wx.StaticText(self.primer_notebook_pane_2, -1, "Min. Flanking PSI:")
        self.psi_text_field = wx.TextCtrl(self.primer_notebook_pane_2, -1, ".95")
        self.psi_text_field.SetToolTip(wx.ToolTip("Valid: 0 < PSI <= 1"))
        self.type_label = wx.StaticText(self.primer_notebook_pane_2, -1, "Splice Junction:")
        self.type_combo_box = wx.ComboBox(self.primer_notebook_pane_2, -1, 'Annotation', choices=["Annotation", "RNA-Seq + Annotation"], style=wx.CB_DROPDOWN | wx.TE_READONLY)
        self.gene_id_label = wx.StaticText(self.primer_notebook_pane_2, -1, "Gene ID:")
        self.gene_id_combo_box = wx.ComboBox(self.primer_notebook_pane_2, -1, 'Valid', choices=["Valid", "Not Valid"], style=wx.CB_DROPDOWN | wx.TE_READONLY)
        self.design_checkbox = wx.CheckBox(self.primer_notebook_pane_2, -1, "Use shortest isoform to design primers?")
        self.temp_label = wx.StaticText(self.primer_notebook_pane_2, -1, "Keep Temporary:")
        self.temp_combo_box = wx.ComboBox(self.primer_notebook_pane_2, -1, 'No', choices=["No", "Yes"], style=wx.CB_DROPDOWN | wx.TE_READONLY)
        self.temp_combo_box.SetToolTip(wx.ToolTip('Keep intermediate files'))
        self.read_threshold_label = wx.StaticText(self.primer_notebook_pane_2, -1, "Read Threshold:")
        self.read_threshold_text_field = wx.TextCtrl(self.primer_notebook_pane_2, -1, "5")
        self.read_threshold_text_field.SetToolTip(wx.ToolTip('Minimum number of splice junction reads to\ndeclare a novel junction from RNA-Seq'))
        self.anchor_length_label = wx.StaticText(self.primer_notebook_pane_2, -1, "Anchor Length:")
        self.anchor_length_text_field = wx.TextCtrl(self.primer_notebook_pane_2, -1, "8")
        self.anchor_length_text_field.SetToolTip(wx.ToolTip('Minimum number of bases on both side of a splice junction\nfor a read to be considered a junction read'))
        self.min_jct_count_label = wx.StaticText(self.primer_notebook_pane_2, -1, "Min. Jct Count:")
        self.min_jct_count_text_field = wx.TextCtrl(self.primer_notebook_pane_2, -1, "1")
        self.min_jct_count_text_field.SetToolTip(wx.ToolTip('The minimum number of counts to assign a junction that\nis found in the gene annotation (GTF).'))
        self.normal_staticbox = wx.StaticBox(self.primer_notebook_pane_2, -1, "Normal")
        self.sizer_11_staticbox = wx.StaticBox(self.primer_notebook_pane_2, -1, "Advanced")

        self.__set_properties()
        self.__do_layout()

        # bind event handlers to different events (IDs represent events)
        self.Bind(wx.EVT_MENU, self.on_help, id=help_id)  # used to specify id as -1
        self.Bind(wx.EVT_MENU, self.on_edit_primer3_path, id=path_id)  # used to specify id as -1
        self.Bind(wx.EVT_MENU, self.on_load_example, id=load_ex_id)  # used to specify id as -1
        self.Bind(wx.EVT_MENU, self.on_load_quick_example, id=load_quick_ex_id)  # used to specify id as -1
        self.Bind(wx.EVT_MENU, self.on_reset, id=reset_id)  # used to specify id as -1
        self.Bind(wx.EVT_MENU, self.on_quit, id=quit_id)  # used to specify id as -1
        self.Bind(wx.EVT_MENU, self.on_add_genes, id=add_genes_id)
        self.Bind(wx.EVT_MENU, self.on_sort_gtf, id=sort_id)
        self.Bind(wx.EVT_MENU, self.primer3_event, id=primer3_id)
        self.Bind(wx.EVT_MENU, self.on_about, id=about_id)
        self.Bind(wx.EVT_MENU, self.on_primer3_manual, id=primer3_manual_id)
        self.Bind(wx.EVT_BUTTON, self.on_choose_fasta_button, self.choose_fasta_button)
        self.Bind(wx.EVT_BUTTON, self.on_choose_gtf_button, self.choose_gtf_button)
        self.Bind(wx.EVT_BUTTON, self.on_choose_bam_button, self.choose_bam_button)
        self.Bind(wx.EVT_BUTTON, self.on_choose_output_button, self.choose_output_button)
        self.Bind(wx.EVT_BUTTON, self.on_run_button, self.run_button)
        self.Bind(Stc.EVT_STC_START_DRAG, self.on_cancel_drag, self.fasta_choice_label)
        # end wxGlade

        self.view_output_frame = None  # No primers have been designed yet

        self.my_icon = wx.EmptyIcon()
        self.my_icon.CopyFromBitmap(wx.Bitmap("PrimerSeq.ico", wx.BITMAP_TYPE_ANY))
        self.SetIcon(self.my_icon)

        self.gtf, self.bam, self.output, self.fasta = [], [], '', None
        pub.subscribe(self.update_after_dialog, "update")
        pub.subscribe(self.update_after_run, "update_after_run")
        pub.subscribe(self.update_after_error, "update_after_error")

        # check if the user has java installed
        try:
            with open(os.devnull, 'wb') as f:
                subprocess.call('java', stdout=f, stderr=f)
        except (subprocess.CalledProcessError, OSError):
            dlg = wx.MessageDialog(self, 'You need java installed on your computer.\nYou can download java from:\n\nhttp://www.oracle.com/technetwork/java/javase/downloads/java-se-jre-7-download-432155.html', style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def __set_properties(self):
        # begin wxGlade: PrimerFrame.__set_properties
        self.SetTitle("PrimerSeq")
        self.primer_frame_statusbar.SetStatusWidths([-1])
        # statusbar fields
        primer_frame_statusbar_fields = [""]  # text in status bar by default
        for i in range(len(primer_frame_statusbar_fields)):
            self.primer_frame_statusbar.SetStatusText(primer_frame_statusbar_fields[i], i)
        self.fasta_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.fasta_choice_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.gtf_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.gtf_choice_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.bam_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.bam_choice_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.coordinates_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.coordinates_text_field.SetMinSize((396, 80))
        self.output_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.output_choice_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.output_choice_label.SetMinSize((135, 27))
        user_entry_size = (135, 27)  # size of widgets in option tab
        self.psi_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.psi_text_field.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.psi_text_field.SetMinSize(user_entry_size)
        self.type_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.type_combo_box.SetMinSize(user_entry_size)
        self.type_combo_box.SetSelection(0)
        self.gene_id_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.gene_id_combo_box.SetMinSize(user_entry_size)
        self.gene_id_combo_box.SetSelection(0)
        self.temp_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.temp_combo_box.SetMinSize(user_entry_size)
        self.temp_combo_box.SetSelection(-1)
        self.read_threshold_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.read_threshold_text_field.SetMinSize(user_entry_size)
        self.anchor_length_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.anchor_length_text_field.SetMinSize(user_entry_size)
        self.min_jct_count_label.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        self.min_jct_count_text_field.SetMinSize(user_entry_size)
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: PrimerFrame.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        sizer_10 = wx.BoxSizer(wx.VERTICAL)
        sizer_with_checkbox = wx.BoxSizer(wx.VERTICAL)
        self.sizer_11_staticbox.Lower()
        self.normal_staticbox.Lower()
        sizer_11 = wx.StaticBoxSizer(self.sizer_11_staticbox, wx.HORIZONTAL)
        normal_box_sizer = wx.StaticBoxSizer(self.normal_staticbox, wx.HORIZONTAL)
        # grid_sizer_4 = wx.GridSizer(4, 2, 0, 0)
        grid_sizer_4 = wx.GridBagSizer(4, 2)
        grid_sizer_checkbox = wx.GridSizer(1, 1, 0, 0)
        # grid_sizer_3 = wx.GridSizer(3, 2, 0, 0)
        grid_sizer_3 = wx.GridBagSizer(4, 2)
        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        sizer_3 = wx.BoxSizer(wx.VERTICAL)
        sizer_8 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_2 = wx.GridSizer(1, 3, 0, 0)
        # sizer_9 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_9 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_4_staticbox.Lower()
        sizer_4 = wx.StaticBoxSizer(self.sizer_4_staticbox, wx.VERTICAL)
        grid_sizer_1 = wx.GridSizer(3, 3, 0, 0)
        grid_sizer_1.Add(self.fasta_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.choose_fasta_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.fasta_choice_label, 1, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)  # NEW LINE
        grid_sizer_1.Add(self.gtf_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.choose_gtf_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.gtf_choice_label, 1, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)  # NEW LINE
        grid_sizer_1.Add(self.bam_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.choose_bam_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.bam_choice_label, 1, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)  # NEW LINE
        sizer_4.Add(grid_sizer_1, 1, wx.EXPAND, 0)
        sizer_3.Add(sizer_4, 1, wx.EXPAND, 0)
        sizer_8.Add(self.coordinates_label, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_8.Add(self.coordinates_text_field, 0, wx.EXPAND, 0)
        grid_sizer_2.Add(self.output_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_2.Add(self.choose_output_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        # sizer_9.Add(self.output_choice_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)
        # self.panel_4.SetSizer(sizer_9)
        # grid_sizer_2.Add(self.panel_4, 1, wx.EXPAND, 0)
        grid_sizer_2.Add(self.output_choice_label, 1, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_8.Add(grid_sizer_2, 1, wx.EXPAND, 0)
        sizer_3.Add(sizer_8, 1, wx.EXPAND, 0)
        sizer_2.Add(sizer_3, 1, wx.EXPAND, 0)
        sizer_2.Add(self.run_button, 0, wx.EXPAND, 0)
        self.primer_notebook_pane_1.SetSizer(sizer_2)
        # grid_sizer_3.Add(self.psi_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_3.Add(self.psi_text_field, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_3.Add(self.type_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_3.Add(self.type_combo_box, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 1)
        # grid_sizer_3.Add(self.gene_id_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_3.Add(self.gene_id_combo_box, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_3.Add(self.psi_label, pos=(0, 0))
        grid_sizer_3.Add(self.psi_text_field, pos=(0, 1))
        grid_sizer_3.Add(self.type_label, pos=(1, 0))
        grid_sizer_3.Add(self.type_combo_box, pos=(1, 1))
        grid_sizer_3.Add(self.gene_id_label, pos=(2, 0))
        grid_sizer_3.Add(self.gene_id_combo_box, pos=(2, 1))
        grid_sizer_3.Add(self.design_checkbox, pos=(3, 0), span=(1,2))
        # sizer_10.Add(grid_sizer_3, 1, wx.EXPAND, 0)
        normal_box_sizer.Add(grid_sizer_3)
        # sizer_10.Add(grid_sizer_3)
        sizer_10.Add(normal_box_sizer, 1, wx.EXPAND, 0)
        # grid_sizer_checkbox.Add(self.design_checkbox, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_4.Add(self.design_checkbox, pos=(0, 0), span=(1,2))
        # grid_sizer_4.Add(self.temp_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_4.Add(self.temp_combo_box, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_4.Add(self.read_threshold_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_4.Add(self.read_threshold_text_field, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_4.Add(self.anchor_length_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_4.Add(self.anchor_length_text_field, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 2)
        # grid_sizer_4.Add(self.min_jct_count_label, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        # grid_sizer_4.Add(self.min_jct_count_text_field, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_4.Add(self.temp_label, pos=(0, 0))
        grid_sizer_4.Add(self.temp_combo_box, pos=(0, 1))
        grid_sizer_4.Add(self.read_threshold_label, pos=(1, 0))
        grid_sizer_4.Add(self.read_threshold_text_field, pos=(1, 1))
        grid_sizer_4.Add(self.anchor_length_label, pos=(2, 0))
        grid_sizer_4.Add(self.anchor_length_text_field, pos=(2, 1))
        grid_sizer_4.Add(self.min_jct_count_label, pos=(3, 0))
        grid_sizer_4.Add(self.min_jct_count_text_field, pos=(3, 1))
        sizer_11.Add(grid_sizer_4, 1, wx.EXPAND, 0)
        # sizer_with_checkbox.Add(grid_sizer_checkbox, 1, wx.EXPAND, 0)
        # sizer_with_checkbox.Add(sizer_11)
        sizer_10.Add(sizer_11, 1, wx.EXPAND, 0)
        # sizer_10.Add(sizer_with_checkbox)
        self.primer_notebook_pane_2.SetSizer(sizer_10)
        self.primer_notebook.AddPage(self.primer_notebook_pane_1, "Required")
        self.primer_notebook.AddPage(self.primer_notebook_pane_2, "Optional")
        sizer_1.Add(self.primer_notebook, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()
        # end wxGlade

    def on_cancel_drag(self, event):
        """Attempt to prevent text from being draggable"""
        event.SetDragText("")
        event.Skip()

    def on_edit_primer3_path(self, event):
        """Open dialog to set primer3 properties"""
        cd.Primer3PathDialog(self, -1, 'Primer3 Location')

    def on_load_quick_example(self, event):
        """Load a single example exon"""
        self.gtf, self.bam, self.fasta = [], [], None
        self.set_fasta('example/chr2.fa', 'chr2.fa', use_dlg=False)
        self.set_bam(['example/example.sorted.bam'], ['example.sorted.bam'], use_dlg=False)
        self.set_gtf('example/example.gtf', 'example.gtf', use_dlg=False)
        self.coordinates_text_field.SetValue('-chr2:216257653-216257926')

    def on_load_example(self, event):
        """Load multiple test exons. The test data can be found in the example directory."""
        self.gtf, self.bam, self.fasta = [], [], None
        self.set_fasta('example/chr2.fa', 'chr2.fa', use_dlg=False)
        self.set_bam(['example/example.sorted.bam'], ['example.sorted.bam'], use_dlg=False)
        self.set_gtf('example/example.gtf', 'example.gtf', use_dlg=False)
        self.coordinates_text_field.SetValue('-chr2:216257653-216257926\n'
                                             '+chr2:228194321-228194499\n'
                                             '-chr2:218694566-218694605\n'
                                             '-chr2:74369398-74369487\n'
                                             '-chr2:38973285-38973321\n'
                                             '+chr2:54035435-54035597\n'
                                             '+chr2:201303844-201304051\n'
                                             '+chr2:173354212-173354386\n'
                                             '+chr2:30371110-30371407\n'
                                             '-chr2:20478343-20478580\n')

    def on_reset(self, event):
        """
        Event handler for user reseting input if a fail/crash happened.
        This method also deletes the tmp directory.
        """
        try:
            if self.load_progress:
                self.load_progress.Destroy()
                self.load_progress = None
        except AttributeError:
            pass
        self.gtf = []
        self.bam = []
        self.fasta = None
        self.bam_choice_label.SetLabel('None')
        self.fasta_choice_label.SetLabel('None')
        self.gtf_choice_label.SetLabel('None')
        self.output_choice_label.SetLabel('None')
        self.coordinates_text_field.SetValue('')
        self.delete_tmp_directory()
        self.enable_load_buttons()  # enable buttons

    def delete_tmp_directory(self):
        """Remove all temporary files"""
        shutil.rmtree(primer.config_options['tmp'])

    def on_help(self, event):
        """Open documentation in default webbrowser"""
        webbrowser.open(os.path.abspath('help/index.html'))

    def on_add_genes(self, event):
        """UCSC add Gene IDs event handler"""
        cd.AddGeneIdsDialog(self, -1, 'Add Valid Gene IDs')

    def on_sort_gtf(self, event):
        """Sort gtf event handler"""
        cd.SortGtfDialog(self, -1, 'Sort GTF')

    def primer3_event(self, event):
        '''Try to open primer3.cfg in every platform so the user can edit it. '''
        filepath = primer.config_options['primer3_cfg']
        if sys.platform.startswith('darwin'):
            subprocess.call(('open', filepath))
        elif os.name == 'nt':
            os.startfile(filepath)
        elif os.name == 'posix':
            subprocess.call(('xdg-open', filepath))

    def on_primer3_manual(self, event):
        '''Try to open primer3_manual.htm in a webbrowser. '''
        primer3_path = primer.config_options['primer3']
        if primer3_path == '../primer3':
            webbrowser.open(os.path.abspath('primer3/primer3_manual.htm'))
        else:
            webbrowser.open(os.path.abspath(primer.config_options['primer3'] + '/primer3_manual.htm'))

    def update_after_dialog(self, msg):
        '''
        Updates attributes and gui components from a started Process
        or thread. This is called by PubSub.
        '''
        try:
            self.load_progress.check_dialog()
        except:
            if msg.data[0] is None:
                self.enable_load_buttons()
            elif isinstance(msg.data[0], int):
                pass
            else:
                for obj, value in msg.data:
                    if isinstance(value, str):
                        # getattr(self, obj).SetLabel(value)
                        getattr(self, obj).SetValue(value)
                    else:
                        setattr(self, obj, value)
                self.enable_load_buttons()
            return

        if msg.data[0] is None:
            self.load_progress.Update(100)
            self.enable_load_buttons()
        elif isinstance(msg.data[0], int):
            perc, text = msg.data
            self.load_progress.Update(perc, text)
        else:
            for obj, value in msg.data:
                if isinstance(value, str):
                    # getattr(self, obj).SetLabel(value)
                    getattr(self, obj).SetValue(value)
                else:
                    setattr(self, obj, value)
            self.load_progress.Update(100)
            self.enable_load_buttons()

    def update_after_run(self, msg):
        """Destory loading dialog and display primer design results"""
        self.load_progress.Destroy()
        self.enable_load_buttons()
        self.view_output_frame = vo.ViewOutputFrame(self, -1, "Primer Design Results", self.options)
        self.view_output_frame.Show()

    def update_after_error(self, msg):
        dlg = wx.MessageDialog(self, 'An uncaught error occured in PrimerSeq. Please check the log file (%s) for details. You may need to press Reset from the File menu to continue.\n\nIf you are consistently having problems please check the PrimerSeq FAQ:\nhttp://primerseq.sourceforge.net/faq.html' % log_file, style=wx.OK | wx.ICON_ERROR)
        dlg.ShowModal()
        self.on_reset(None)

    def disable_load_buttons(self):
        """Disable loading file/run buttons when action in progress"""
        self.choose_bam_button.Disable()
        self.choose_fasta_button.Disable()
        self.choose_gtf_button.Disable()
        self.run_button.Disable()

    def enable_load_buttons(self):
        """Enable loading file/run buttons when action is finished"""
        self.choose_bam_button.Enable()
        self.choose_fasta_button.Enable()
        self.choose_gtf_button.Enable()
        self.run_button.Enable()

    def process_bam(self, fnames, fnames_without_path, anc_len):
        """This method is threaded by the set_bam method"""
        tmp_bam = []  # a list of sam.Sam obj to be returned
        for i, f in enumerate(fnames):
            wx.CallAfter(pub.sendMessage, "update", (int(float(i) / len(fnames) * 100), 'Reading %s . . .' % fnames_without_path[i]))
            tmp_bam.append(sam.Sam(f, anc_len))
        return tmp_bam

    def on_quit(self, event):  # wxGlade: PrimerFrame.<event_handler>
        """Quit PrimerSeq when user presses File -> Quit"""
        self.Destroy()
        event.Skip()

    def on_choose_output_button(self, event):  # wxGlade: PrimerFrame.<event_handler>
        """Event handler for setting output text file"""
        dlg = wx.FileDialog(self, message='Choose your output file',
                            wildcard='Text file (*.txt)|*.txt')  # open file dialog

        if dlg.ShowModal() == wx.ID_OK:
            self.set_output(dlg.GetPath(), dlg.GetFilename())
            dlg.Destroy()  # best to do this sooner
        event.Skip()

    def set_output(self, path, filename):
        """Set the output file location for PrimerSeq"""
        self.output = path
        # self.output_choice_label.SetLabel(filename)
        self.output_choice_label.SetValue(filename)

    def on_choose_fasta_button(self, event):  # wxGlade: PrimerFrame.<event_handler>
        """FASTA button event handler"""
        dlg = wx.FileDialog(self, message='Choose your FASTA file', defaultDir=os.getcwd(),
                            wildcard='FASTA file (*.fa)|*.fa|FASTA file(*.fasta)|*.fasta')  # open file dialog
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()  # get the new filenames from the dialog
            filename_without_path = dlg.GetFilename()
            dlg.Destroy()  # best to do this sooner
            self.set_fasta(filename, filename_without_path)
        else:
            dlg.Destroy()  # make sure to destroy if they hit cancel
        event.Skip()

    def set_fasta(self, filename, filename_without_path, use_dlg=True):
        """
        Loads FASTA into PrimerSeq. Uses Pygr's SequenceFileDB object to index the FASTA file
        so that subsequent loads occur almost instantly while still maintaining random access.
        Indexing is preformed within a thread to prevent the GUI from locking.

        :param str filename: full path to fasta
        :param str filename_without_path: filename without full path (for display purposes)
        :param bool use_dlg: flag for using a pop up dialog
        """
        try:
            # set the fasta attribute
            if use_dlg:
                self.load_progress = cd.CustomDialog(self, -1, 'FASTA', 'Loading FASTA . . .\n\nThis will take several minutes')
                self.load_progress.Update(0)
            self.disable_load_buttons()  # disable loading other files while another is running
            self.current_process = ct.RunThread(target=SequenceFileDB,
                                                args=(str(filename),),
                                                attr='fasta', label='fasta_choice_label', label_text=str(filename_without_path))
        except:
            if use_dlg:
                self.load_progress.Destroy()
            t, v, trace = sys.exc_info()
            print('ERROR! For more information read the following lines')
            print('Type: ' + str(t))
            print('Value: ' + str(v))
            print('Traceback:\n' + traceback.format_exc())

    def on_choose_gtf_button(self, event):  # wxGlade: PrimerFrame.<event_handler>
        """GTF button handler"""
        dlg = wx.FileDialog(self, message='Choose your GTF file', defaultDir=os.getcwd(),
                            wildcard='GTF file (*.gtf)|*.gtf')  # open file dialog
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()  # get the new filenames from the dialog
            dlg.Destroy()  # best to do this sooner
            filename_without_path = dlg.GetFilename()  # only grab the actual filenames and none of the path information
            self.set_gtf(filename, filename_without_path)
        else:
            dlg.Destroy()  # make sure to destroy if they hit cancel
        event.Skip()

    def set_gtf(self, filename, filename_without_path, use_dlg=True):
        """
        Loads exon features from GTF file into a python data structure. Uses a threaded call
        to :func:`~primer.gene_annotation_reader`.

        :param str filename: full path to GTF
        :param str filename_without_path: filename without full path (for display purposes)
        :param bool use_dlg: flag for using a pop up dialog
        """
        # set the gtf attribute
        if use_dlg:
            self.load_progress = cd.CustomDialog(self, -1, 'GTF', 'Loading GTF . . .\n\nThis will take ~1 min.')
            self.load_progress.Update(0)
        self.disable_load_buttons()
        self.current_process = ct.RunThread(target=primer.gene_annotation_reader, args=(str(filename),),
                                            attr='gtf', label='gtf_choice_label', label_text=str(filename_without_path))

    def on_choose_bam_button(self, event):  # wxGlade: PrimerFrame.<event_handler>
        """Event handler for loading bam file"""
        dlg = wx.FileDialog(self, message='Choose your bam files', defaultDir=os.getcwd(),
                            wildcard='BAM files (*.bam)|*.bam|SAM files (*.sam)|*.sam', style=wx.FD_MULTIPLE)  # open file dialog
        # if they press ok
        if dlg.ShowModal() == wx.ID_OK:
            filenames = dlg.GetPaths()  # get the new filenames from the dialog
            filenames_without_path = dlg.GetFilenames()  # only grab the actual filenames and none of the path information
            dlg.Destroy()  # best to do this sooner

            self.set_bam(filenames, filenames_without_path)  # load bam file
        else:
            dlg.Destroy()  # make sure to destroy if they hit cancel
        event.Skip()

    def set_bam(self, filenames, filenames_without_path, use_dlg=True):
        """
        Ultimately uses the SAM-JDK to set the BAM file. The BAM file is not loaded into
        memory. However, if there is no .sorted.bam extension then the user's mapped reads
        will be converted to a BAM file and then is sorted (For large files this will take a long time).
        Uses a threaded call to :meth:`~PrimerFrame.process_bam`.

        :param str filename: full path to BAM
        :param str filename_without_path: filename without full path (for display purposes)
        :param bool use_dlg: flag for using a pop up dialog
        """
        # set the bam attribute
        self.bam = []  # clear bam attribute
        if use_dlg:
            self.load_progress = cd.CustomDialog(self, -1, 'BAM', 'Loading BAM/SAM . . .\n\nThis may take several minutes')
            self.load_progress.Update(0)
        self.disable_load_buttons()
        self.current_process = ct.RunThread(target=self.process_bam,
                                            args=(map(str, filenames), filenames_without_path, int(self.anchor_length_text_field.GetValue())),
                                            attr='bam', label='bam_choice_label', label_text=str(', '.join(filenames_without_path)))

    def on_run_button(self, event):  # wxGlade: PrimerFrame.<event_handler>
        """
        Event handler for the running PrimerSeq button. Gets all necessary information from
        GUI to start primer design. Then uses a threaded call to :meth:`~PrimerFrame.run_primer_design`.
        """
        # alert the user there is missing input
        if self.gtf == [] or self.fasta is None or self.output == '':
            dlg = wx.MessageDialog(self, 'Please fill in all of the required fields.', style=wx.OK)
            dlg.ShowModal()
            return
        if self.view_output_frame:
            dlg = wx.MessageDialog(self, 'Please close your existing primer design results\n'
                                   'before designing more primers.', style=wx.OK)
            dlg.ShowModal()
            return

        strandList, chrList, startList, endList = [], [], [], []  # stores all coordinate info

        # handle the coordinates in self.coordinates_text_input
        coordinates_string = self.coordinates_text_field.GetValue()  # a string
        coordinates = map(lambda y: re.split('\s*,+\s*', y.strip()),
                          map(str, filter(lambda x: x != '',
                          re.split('\s*\n+\s*', coordinates_string))))  # ['(strand)(chr):(start)-(end)', ...]
        if not coordinates:
            dlg = wx.MessageDialog(self, 'Please fill in the exon coordinates field.', style=wx.OK)
            dlg.ShowModal()
            return

        # make sure exons actually exist in GTF before preceeding
        bad_coordinate = []
        for c in coordinates:
            found_ex = self.find_exon(c)
            bad_coordinate += [f == False for f in found_ex]
        if (len(coordinates) == 1 and any(bad_coordinate)) or (not all(bad_coordinate) and any(bad_coordinate)):
            clist = [item for sublist in coordinates for item in sublist]
            first_indx = bad_coordinate.index(True)
            first_bad_coord = clist[first_indx]
            dlg = wx.MessageDialog(self, 'One or more coordinates were incorrect. First occurence of bad coordinate'
                                   ' is %s.' % first_bad_coord, style=wx.OK)
            dlg.ShowModal()
            return
        elif all(bad_coordinate):
            dlg = wx.MessageDialog(self, 'All entered coordinates do not match an exon in the gene annotation.'
                                   'You may have accidentally enter coordinates in 1-based format.', style=wx.OK)
            dlg.ShowModal()
            return


        # options for primer.py's main function
        self.options = {}
        self.options['target'] = zip(range(1, len(coordinates) + 1), coordinates)
        self.options['gtf'] = self.gtf
        self.options['fasta'] = self.fasta
        self.options['rnaseq'] = self.bam
        self.options['psi'] = float(self.psi_text_field.GetValue())
        self.options['rnaseq_flag'] = False
        self.options['annotation_flag'] = True if str(self.type_combo_box.GetValue()) == 'Annotation' else False
        self.options['both_flag'] = True if str(self.type_combo_box.GetValue()) == 'RNA-Seq + Annotation' else False
        self.options['output'] = self.output
        self.options['read_threshold'] = int(self.read_threshold_text_field.GetValue())
        self.options['keep_temp'] = False if str(self.temp_combo_box.GetValue()) == 'No' else True
        self.options['big_bed'] = None
        self.options['no_gene_id'] = False if str(self.gene_id_combo_box.GetValue()) == 'Valid' else True
        self.options['min_jct_count'] = int(self.min_jct_count_text_field.GetValue())
        self.options['anchor_length'] = int(self.anchor_length_text_field.GetValue())
        self.options['job_id'] = 'PrimerSeq'
        self.options['short_isoform'] = self.design_checkbox.IsChecked()

        # display dialog and disable buttons while designing primers
        self.load_progress = cd.CustomDialog(self, -1, 'Run PrimerSeq', 'Designing primers . . .\n\nThis dialog will close after it is done.')
        self.load_progress.Update(0)
        self.disable_load_buttons()

        # Design primers using a thread, otherwise the GUI will lock up
        self.current_process = ct.RunPrimerSeqThread(target=self.run_primer_design,  # run primer design by calling primer.py's main function
                                                     args=(self.options,))
        event.Skip()

    def find_exon(self, coord):
        """
        Check if coordinates match a known exon from GTF. This helps alert
        the user if they enter bad coordinates.

        :param list coord: a list of coordinates [-chr8:100-1000, ...]
        :return: a list of booleans of whether coordinate was found in GTF
        """
        found = [False] * len(coord)
        chr = utils.get_chr(coord[0][1:])
        strand = coord[0][0]
        pos_list = [utils.get_pos(c) for c in coord]
        tmp = [p for sublist in pos_list for p in sublist]
        start, end = min(tmp), max(tmp)
        for gene_key in self.gtf[chr]:
            if self.gtf[chr][gene_key]['strand'] == strand and self.gtf[chr][gene_key]['start'] <= start and self.gtf[chr][gene_key]['end'] >= end:
                for k, pos in enumerate(pos_list):
                    if pos in self.gtf[chr][gene_key]['exons']:
                        found[k] = True
        return found

    def run_primer_design(self, opts):
        """
        Performs that actual execution of primer design followed by
        cleaning up some of the files in the tmp directory.

        :param dict opts: user options for running PrimerSeq
        """
        # design primers
        primer.main(opts)

        # make bam location file
        with open('tmp/indiv_isoforms/key.txt', 'wb') as my_writer:
            for s in opts['rnaseq']:
                my_writer.write(s.path + '\n')

        # remove old json files
        for f in glob.glob('tmp/indiv_isoforms/*.json'):
            os.remove(f)

        # save data in json format
        with open(opts['output']) as handle:
            handle.readline()  # skip the header
            for line in csv.reader(handle, delimiter='\t'):
                # lines of length 1 are displaying error msgs
                if len(line) > 1:
                    rc.save_isforms_and_counts(line, opts)

        # remove SAM files
        for f in glob.glob(os.path.join(primer.config_options['tmp'], 'sam/*.sam')):
            os.remove(f)

    def on_about(self, event):
        """Displays the About dialog box which explains details about PrimerSeq to the user"""

        description = """PrimerSeq designs RT-PCR primers that evaluate alternative splicing events by incorporating RNA-Seq data.
        PrimerSeq is particularly advantageous for designing a large number of primers for validating alternative splicing events found in RNA-Seq data.
        PrimerSeq incorporates RNA-Seq data in the design process to weight exons by their read counts.
        Essentially, the RNA-Seq data allows primers to be placed using actually expressed transcripts.
        This could be for your particular cell line or experimental condition, rather than using annotations that incorporate transcripts that are not expressed for your data.
        Alternatively, you can design primers that are always on constitutive exons. PrimerSeq does not limit the use of gene annotations and can be used for a wide array of species.

PrimerSeq redistributes primer3 which is licensed under GPLv2, the SAM-JDK which is licensed under Apache License V2.0, MIT, and the BigWig api which
is licensed under LGPL v2.1. There is no source code modification to any of the previous work."""

        licence = """Copyright (C) 2012-2013  Collin Tokheim

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>."""

        info = wx.AboutDialogInfo()
        tmp_icon = wx.EmptyIcon()
        tmp_icon.CopyFromBitmap(wx.Bitmap("about.ico", wx.BITMAP_TYPE_ANY))
        info.SetName('PrimerSeq')
        info.SetVersion('1.1.5.beta')
        info.SetDescription(description)
        info.SetCopyright('(C) 2012-2013 Collin Tokheim')
        info.SetWebSite('http://primerseq.sf.net')
        info.SetLicence(licence)
        info.AddDeveloper('Collin Tokheim')
        info.AddDocWriter('Collin Tokheim')
        info.SetIcon(tmp_icon)
        wx.AboutBox(info)


# end of class PrimerFrame
class PrimerApp(wx.App):
    def OnInit(self):
        wx.InitAllImageHandlers()
        self.primer_frame = PrimerFrame(None, -1, "")
        self.SetTopWindow(self.primer_frame)
        self.primer_frame.Show()
        return 1

# end of class PrimerApp

if __name__ == "__main__":
    # handle all uncaught exceptions
    sys.excepthook = handle_uncaught_exceptions

    # define logging file before using logging.debug
    if not os.path.exists(primer.config_options['log']): os.mkdir(primer.config_options['log'])  # make directory to put log files
    log_file = os.path.join(primer.config_options['log'], 'log.PrimerApp.%s.txt' % str(datetime.datetime.now()).replace(':', '.'))
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(message)s',
                        filename=log_file,
                        filemode='w')

    # make sure primer3_log directory exists
    if not os.path.exists(primer.config_options['primer3_log']): os.mkdir(primer.config_options['primer3_log'])

    # start GUI
    primerApp = PrimerApp(redirect=False)
    primerApp.MainLoop()
