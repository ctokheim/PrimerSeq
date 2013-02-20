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
import sys
import csv
import custom_dialog as cd
import utils
import traceback  # debugging import


class ViewOutputFrame(wx.Frame, listmix.ColumnSorterMixin):
    def __init__(self, parent, id, string, opts):
        # wx.Dialog.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        wx.Frame.__init__(self, parent, -1, string)
        # self.output_filename = output_file_to_load
        self.options = opts
        self.output_filename = opts['output']  # separate attribute for legacy reasons

        # ToolBar at the top of the window
        toolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL | wx.NO_BORDER)
        toolbar.SetToolBitmapSize(size=(24, 24))
        plot_id = wx.NewId()
        toolbar.AddSimpleTool(plot_id, self.get_bmp(wx.ART_NEW),
                              "Create Plots", "Create new primer plots")
        insilico_pcr_id = wx.NewId()
        toolbar.AddSimpleTool(insilico_pcr_id, self.get_bmp(wx.ART_EXECUTABLE_FILE),
                              "In Silico PCR", "Test primers with in silico PCR")
        save_plot_id = wx.NewId()
        toolbar.AddSimpleTool(save_plot_id, self.get_bmp(wx.ART_FILE_SAVE),
                              "Save Plots", "Save each plot to an html file viewable in your browser")
        reset_id = wx.NewId()
        toolbar.AddSimpleTool(reset_id, self.get_bmp(wx.ART_UNDO),
                              "Reset Values", "Reset Information")
        toolbar.AddSimpleTool(wx.ID_HELP, self.get_bmp(wx.ART_HELP),
                              "Help", "Help Information")
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(wx.ID_EXIT, self.get_bmp(wx.ART_QUIT),
                              "Exit", "Exit the program")
        toolbar.Realize()
        self.SetToolBar(toolbar)

        tID = wx.NewId()
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = utils.MyListCtrl(self, tID,
                                     style=wx.LC_REPORT
                                     | wx.LC_SORT_ASCENDING
                                     | wx.LC_EDIT_LABELS
                                     | wx.BORDER_NONE)

        # self.list.SetImageList(self.il, wx.IMAGE_LIST_SMALL)
        sizer.Add(self.list, 1, wx.EXPAND)

        # define the columns
        self.list.InsertColumn(0, 'ID')
        self.list.InsertColumn(1, 'Target Coordinate')
        self.list.InsertColumn(2, 'Primer Coordinates')
        self.list.InsertColumn(3, 'Psi Target')
        self.list.InsertColumn(4, 'Upstream Primer')
        self.list.InsertColumn(5, 'Downstream Primer')
        self.list.InsertColumn(6, 'Average TM')
        self.list.InsertColumn(7, 'Skipping Prod. Size')
        self.list.InsertColumn(8, 'Inc. Prod. Sized')
        self.list.InsertColumn(9, 'Upstream Exon Coord.')
        self.list.InsertColumn(10, 'Upstream Psi')
        self.list.InsertColumn(11, 'Downstream Exon Coord.')
        self.list.InsertColumn(12, 'Downstream Psi')
        self.list.InsertColumn(13, 'ASM Region')
        self.list.InsertColumn(14, 'Gene')
        self.read_output_file(self.output_filename)  # read in contents from file

        # Now that the list exists we can init the other base class,
        # see wx/lib/mixins/listctrl.py
        listmix.ColumnSorterMixin.__init__(self, 3)
        #self.SortListItems(0, True)

        self.SetSizer(sizer)
        self.SetAutoLayout(True)

        # bind toolbar icons
        self.Bind(wx.EVT_TOOL, self.on_insilico_pcr, id=insilico_pcr_id)
        self.Bind(wx.EVT_TOOL, self.on_save_plots, id=save_plot_id)
        self.Bind(wx.EVT_TOOL, self.on_plot, id=plot_id)
        self.Bind(wx.EVT_TOOL, self.on_help, id=wx.ID_HELP)
        self.Bind(wx.EVT_TOOL, self.on_reset, id=reset_id)
        self.Bind(wx.EVT_TOOL, self.on_exit, id=wx.ID_EXIT)

    def on_reset(self, event):
        """Refreshes the content of the listctrl in case the user edits it and wants
        the original data back."""
        self.list.DeleteAllItems()
        self.read_output_file(self.output_filename)

    def on_save_plots(self, event):
        try:
            dlg = cd.SavePlotDialog(self, -1, "Generate HTML Report", self.options)
            dlg.ShowModal()
        except Exception, e:
            print traceback.format_exc()

    def on_help(self, event):
        dlg = wx.MessageDialog(self, 'Instructions:\n\nUse the tool bar to validate '
                                     'that the primer design was successful.\n\nPress '
                                     'the Create Plots button to view how your data '
                                     'supports the primer design results\n\nPress the '
                                     'In Silico Pcr button to validate primer products\n\n'
                                     'Press the save plots button to save plots as html', style=wx.OK)
        dlg.ShowModal()

    def on_exit(self, event):
        self.Destroy()

    def on_insilico_pcr(self, event):
        dlg = cd.InSilicoPcrDialog(self, -1, "In Silico PCR",
                                   self.output_filename)
        dlg.ShowModal()

    def on_plot(self, event):
        cd.PlotDialog(self, -1, 'Plot Results', self.output_filename)

    def read_output_file(self, output_file):
        # get information from results file
        with open(output_file) as handle:
            self.results = list(csv.reader(handle, delimiter='\t'))[1:]

        for i, data in enumerate(self.results):
            index = self.list.InsertStringItem(sys.maxint, data[0])
            for col in range(1, len(data)):
                self.list.SetStringItem(index, col, data[col])
            self.list.SetItemData(index, i)
            if i % 2:
                # zebra stiped color pattern
                self.list.SetItemBackgroundColour(i, "lightgray")

        self.list.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self.list.SetColumnWidth(1, wx.LIST_AUTOSIZE)
        self.list.SetColumnWidth(2, 100)

        # show how to select an item
        self.list.SetItemState(5, wx.LIST_STATE_SELECTED, wx.LIST_STATE_SELECTED)

        self.current_selection = 0

    def get_bmp(self, pic_id):
        return wx.ArtProvider.GetBitmap(pic_id, wx.ART_TOOLBAR, wx.Size(24, 24))

    def GetListCtrl(self):
        return self.list

    def GetSortImages(self):
        return (self.sm_dn, self.sm_up)


class TestViewOutputApp(wx.App):
    '''
    This is just a test app to see if this works indepently
    '''
    def OnInit(self):
        wx.InitAllImageHandlers()
        view_output = ViewOutputFrame(None, -1, "", 'output.txt')
        self.SetTopWindow(view_output)
        view_output.Show()
        return 1

if __name__ == '__main__':
    ViewOutputApp = TestViewOutputApp(0)
    ViewOutputApp.MainLoop()
